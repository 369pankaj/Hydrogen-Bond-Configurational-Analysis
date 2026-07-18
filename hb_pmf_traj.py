#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hb_pmf_traj.py -- build the P(r,theta) / P_rand(r,theta) grid for a donor-acceptor
pair directly from an MD trajectory (MDAnalysis), reproducing the Fortran front-end
of this project and the methodology of Sapir & Harries (JCTC 2017, 13, 2851).

For every donor(heavy)-acceptor(heavy) pair within `rmax`:
    r     = |D - A|                       (minimum-image)
    theta = min over D's hydrogens H of  angle( D->A , D->H )     (0 = linear H-bond)
In `--symmetric` mode (like water-water) the acceptor's hydrogens are also considered
with the reversed axis A->D, matching the Fortran W-W code.

The observed density P(r,theta) is the normalised 2-D histogram of (r,theta).  The
random reference is P_rand(r,theta) = P_rand(r) . P_rand(theta), with
    P_rand(r)     = shell-volume ratio  [ (4/3)pi(r_hi^3 - r_lo^3) ] / [ (4/3)pi rmax^3 ]
    P_rand(theta) = read from --prand-theta file (as the Fortran does) or sin(theta).
An optional g(r) (--rdf) supplies the bulk correction F_avg for the PMF column only;
it CANCELS in the dG (Eq. 2) computed by hb_strength.py, so it is not needed for
strengths.

OUTPUT is the 6-column grid (r theta P Prand PMF P*PMF) consumed by hb_strength.py.

EXAMPLES
  # water-water (symmetric), reusing the project's random-angle distribution:
  python hb_pmf_traj.py --top ../../md.gro --traj ../../md.xtc \
      --donor-heavy "resname SOL and name OW" --acceptor "resname SOL and name OW" \
      --symmetric --exclude-same-residue --prand-theta ../P_rand_WW.dat \
      --stride 200 -o WW_grid.dat

  # urea carbonyl-O (acceptor) with water-O (donor)  -> Ou-W :
  python hb_pmf_traj.py --top ../../md.gro --traj ../../md.xtc \
      --donor-heavy "resname SOL and name OW" --acceptor "resname UREA and name OU" \
      --prand-theta ../P_rand_WH.dat --stride 200 -o OuW_grid.dat

  # urea amine-N (donor) with water-O (acceptor)     -> Nu-W :
  python hb_pmf_traj.py --top ../../md.gro --traj ../../md.xtc \
      --donor-heavy "resname UREA and name N1U N2U" --acceptor "resname SOL and name OW" \
      --prand-theta ../P_rand_WW.dat --stride 200 -o NuW_grid.dat
"""
import argparse
import os
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors, capped_distance

R_KJ = 0.0083144621


def build_donor_h_map(u, donor_heavy, hcut=1.3):
    """Map each donor heavy atom (global index) -> list of covalently bound H indices
    (same molecule, within hcut Angstrom on the first frame, minimum-image).

    Hydrogen candidates are taken from the donor's own residue (``d.residue.atoms``)
    rather than by (resid, segid) matching: a .gro topology gives every molecule one
    default segid and reuses resids across species, so a (resid, segid) filter would
    pull in hydrogens of unrelated molecules. The heavy-atom -> H distance uses the
    minimum-image convention so a donor whose molecule is wrapped across a box face on
    frame 0 still finds its hydrogens instead of being silently dropped for the run."""
    box = u.dimensions
    dmap = {}
    for d in donor_heavy:
        same = d.residue.atoms
        names = same.names.astype(str)
        is_h = np.char.startswith(names, "H") | (same.types.astype(str) == "H")
        hyd = same[is_h]
        if len(hyd) == 0:
            dmap[d.index] = []
            continue
        dv = minimize_vectors(hyd.positions - d.position[None, :], box=box)
        dist = np.linalg.norm(dv, axis=1)
        dmap[d.index] = [int(idx) for idx in hyd.indices[dist < hcut]]
    return dmap


def build_h_array(group, dmap, hmax=None):
    """Pack each atom's donor-H global indices into a fixed-width int array aligned to
    ``group`` (row i = the H indices of group atom i, -1 padded). Returns (arr, hmax).
    This lets the acceptor-side angle sweep be vectorised over pairs, one H-slot at a
    time, instead of a Python loop over every pair."""
    lists = [dmap.get(int(idx), []) for idx in group.indices]
    if hmax is None:
        hmax = max((len(l) for l in lists), default=0)
    arr = np.full((len(lists), hmax), -1, dtype=np.int64)
    for i, l in enumerate(lists):
        if l:
            arr[i, :len(l)] = l[:hmax]
    return arr, hmax


def angle_deg(v1, v2):
    """Angle (deg) between rows of v1 and rows of v2."""
    c = np.einsum("ij,ij->i", v1, v2) / (
        np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-30)
    c = np.clip(c, -1.0, 1.0)
    return np.degrees(np.arccos(c))


def _neighbor_map(dpos, apos, rmax, box, n_d):
    """For each donor (row of dpos) return an int array of acceptor local indices
    within rmax, using a capped neighbour search instead of the full N_d x N_a
    distance matrix. Selection is identical to `r <= rmax` on the full matrix --
    this only avoids computing distances to acceptors that are far away, which is
    the difference between tractable and intractable on a solvated protein."""
    nm = [np.empty(0, dtype=np.int64) for _ in range(n_d)]
    pairs = capped_distance(dpos, apos, max_cutoff=rmax, box=box,
                            return_distances=False)
    if len(pairs):
        order = np.argsort(pairs[:, 0], kind="stable")
        p = pairs[order]
        uniq, starts = np.unique(p[:, 0], return_index=True)
        chunks = np.split(p[:, 1].astype(np.int64), starts[1:])
        for di_local, acc_local in zip(uniq, chunks):
            nm[di_local] = acc_local
    return nm


def accumulate(u, donor_heavy, acceptor, dmap, acc_harr, acc_hmax,
               box_fn, args, hist, edges_r, edges_t):
    """Histogram (r,theta) over the trajectory into `hist` (in place)."""
    d_idx = donor_heavy.indices
    a_idx = acceptor.indices
    d_res = donor_heavy.resindices     # global (unique) residue index per atom
    a_res = acceptor.resindices
    # Molecule identity via bonded fragments (needs a topology with bonds, e.g. a
    # .tpr). Required for polymers: a protein is many residues in ONE molecule, so
    # resindex cannot express "same chain". Only accessed if requested.
    d_frag = donor_heavy.fragindices if args.exclude_same_molecule else None
    a_frag = acceptor.fragindices if args.exclude_same_molecule else None
    n_d = len(d_idx)
    all_acc = np.arange(len(a_idx))    # candidate set in the full-matrix path

    nfr = 0
    for ts in u.trajectory[args.start:args.stop:args.stride]:
        box = box_fn(ts)
        pos = u.atoms.positions
        dpos = pos[d_idx]
        apos = pos[a_idx]
        nmap = _neighbor_map(dpos, apos, args.rmax, box, n_d) \
            if args.neighbor_list else None
        for i, di in enumerate(d_idx):
            # cand = acceptor-group positions to test for this donor (all, or the
            # rmax neighbour shortlist). Everything downstream is indexed through it.
            cand = nmap[i] if nmap is not None else all_acc
            if len(cand) == 0:
                continue
            dvec = minimize_vectors(apos[cand] - dpos[i], box=box)   # D->A
            r = np.linalg.norm(dvec, axis=1)
            sel = r <= args.rmax
            if args.exclude_same_residue:
                # Same-molecule test on the *global* residue index, so it is
                # unambiguous even when resids are reused across species (a .gro
                # topology) -- unlike a bare resid comparison, which would drop
                # valid inter-molecular pairs that share a resid number.
                sel &= (a_res[cand] != d_res[i])
            if args.exclude_same_molecule:
                # Same bonded molecule (fragment): the correct exclusion for
                # proteins/polymers, where "same residue" is too narrow.
                sel &= (a_frag[cand] != d_frag[i])
            if args.symmetric:
                # Symmetric mode considers both partners' hydrogens for a pair, so
                # (A,B) and (B,A) are physically identical; keep each unordered pair
                # once (a_idx > di) instead of enumerating -- and double-counting --
                # both orderings.
                sel &= (a_idx[cand] > di)
            else:
                sel &= (a_idx[cand] != di)
            if not sel.any():
                continue
            loc = np.where(sel)[0]      # positions within `cand`
            jpos = cand[loc]            # positions within the acceptor group
            rr = r[loc]
            dv = dvec[loc]                                        # D->A, (m,3)
            # donor-side angles: min over this donor's hydrogens of angle(D->A, D->H)
            ang_stack = []
            for h in dmap.get(int(di), []):
                dh = minimize_vectors(pos[h] - dpos[i][None, :], box=box)  # D->H
                ang_stack.append(angle_deg(np.repeat(dh, len(jpos), axis=0), dv))
            if args.symmetric:
                # acceptor-side: for each selected acceptor use its own H's, axis A->D.
                # Vectorised one H-slot at a time over all selected pairs; slots with
                # no hydrogen contribute +inf so they are ignored by the final min.
                apos_sel = apos[jpos]
                for s in range(acc_hmax):
                    hg = acc_harr[jpos, s]                        # global H idx / -1
                    valid = hg >= 0
                    row = np.full(len(jpos), np.inf)
                    if valid.any():
                        ah = minimize_vectors(pos[hg[valid]] - apos_sel[valid],
                                              box=box)             # A->H
                        row[valid] = angle_deg(ah, -dv[valid])    # angle(A->H, A->D)
                    ang_stack.append(row)
            if not ang_stack:
                continue
            theta = np.min(np.vstack(ang_stack), axis=0)
            m = (theta <= args.theta_max) & np.isfinite(theta)
            if m.any():
                h2, _, _ = np.histogram2d(rr[m], theta[m], bins=[edges_r, edges_t])
                hist += h2
        nfr += 1
    return nfr


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build P(r,theta)/P_rand grid from a trajectory "
                    "[Sapir & Harries, JCTC 2017].",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--top", required=True, help="topology (gro/tpr/pdb)")
    ap.add_argument("--traj", required=True, help="trajectory (xtc/trr/...)")
    ap.add_argument("--donor-heavy", required=True,
                    help="MDAnalysis selection for donor heavy atoms")
    ap.add_argument("--acceptor", required=True,
                    help="MDAnalysis selection for acceptor heavy atoms")
    ap.add_argument("--symmetric", action="store_true",
                    help="also use acceptor's H (axis A->D); use for like-like (W-W)")
    ap.add_argument("--exclude-same-residue", action="store_true",
                    help="skip pairs on the same residue (small molecules)")
    ap.add_argument("--exclude-same-molecule", action="store_true",
                    help="skip pairs on the same bonded molecule/fragment; use for "
                         "polymers/proteins where a molecule spans many residues "
                         "(needs a topology with bonds, e.g. .tpr)")
    ap.add_argument("--neighbor-list", action="store_true",
                    help="use a capped neighbour search instead of the full "
                         "donor x acceptor distance matrix (much faster for large "
                         "solvated systems; identical result)")
    ap.add_argument("--prand-theta", default=None,
                    help="two-column file (theta_deg  P_rand(theta)); "
                         "if omitted, sin(theta) is used")
    ap.add_argument("-T", "--temperature", type=float, default=300.0,
                    help="temperature (K) for the PMF column; match hb_strength.py -T")
    ap.add_argument("--rmax", type=float, default=4.0)
    ap.add_argument("--dr", type=float, default=0.02)
    ap.add_argument("--dtheta", type=float, default=1.0)
    ap.add_argument("--theta-max", type=float, default=100.0)
    ap.add_argument("--shell", nargs=2, type=float, metavar=("RMIN", "RMAX"),
                    default=(2.49, 3.25), help="first-shell range for the _1st_shell file")
    ap.add_argument("--rdf", default=None,
                    help="optional g(r) file (r  g) for the F_avg PMF bulk correction")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("-o", "--out", required=True, help="output grid file")
    args = ap.parse_args(argv)

    u = mda.Universe(args.top, args.traj)
    donor_heavy = u.select_atoms(args.donor_heavy)
    acceptor = u.select_atoms(args.acceptor)
    if len(donor_heavy) == 0 or len(acceptor) == 0:
        sys.exit("empty donor or acceptor selection -- check --donor-heavy/--acceptor")
    print(f"donors:  {len(donor_heavy)} atoms  ({args.donor_heavy})")
    print(f"accept.: {len(acceptor)} atoms  ({args.acceptor})")

    if args.exclude_same_molecule and not hasattr(u.atoms, "fragindices"):
        sys.exit("--exclude-same-molecule needs bond information, but the topology "
                 f"'{args.top}' has none. Use a .tpr/.psf (or a .pdb with CONECT), "
                 "or fall back to --exclude-same-residue for single-residue molecules.")

    dmap = build_donor_h_map(u, donor_heavy)
    if args.symmetric:
        dmap.update(build_donor_h_map(u, acceptor))
    nh = np.mean([len(v) for v in dmap.values()]) if dmap else 0
    print(f"mean H per donor heavy atom: {nh:.2f}")

    # padded acceptor-H index array (only used in symmetric mode)
    acc_harr, acc_hmax = build_h_array(acceptor, dmap) if args.symmetric \
        else (np.empty((len(acceptor), 0), dtype=np.int64), 0)

    edges_r = np.arange(0.0, args.rmax + args.dr, args.dr)
    edges_t = np.arange(0.0, args.theta_max + args.dtheta, args.dtheta)
    hist = np.zeros((edges_r.size - 1, edges_t.size - 1))

    box_fn = lambda ts: ts.dimensions
    nfr = accumulate(u, donor_heavy, acceptor, dmap, acc_harr, acc_hmax,
                     box_fn, args, hist, edges_r, edges_t)
    print(f"frames processed: {nfr}")

    # --- normalise P so that sum(P*dr*dtheta) = 1  (matches Fortran nf) ---
    r_c = 0.5 * (edges_r[:-1] + edges_r[1:])
    t_c = 0.5 * (edges_t[:-1] + edges_t[1:])
    P = hist / (hist.sum() * args.dr * args.dtheta + 1e-300)

    # --- P_rand(r): shell-volume density, normalised like P (sum*dr == 1) ---
    # v_r/v_tot is a per-bin probability; divide by dr to get a density so that
    # Prand integrates to 1 consistently with P (and prand_t). Without this the
    # PMF column carries a spurious -RT*ln(dr) offset when --rdf is not given.
    v_tot = (4.0 / 3.0) * np.pi * args.rmax ** 3
    v_r = (4.0 / 3.0) * np.pi * (edges_r[1:] ** 3 - edges_r[:-1] ** 3)
    prand_r = (v_r / v_tot) / args.dr

    # --- P_rand(theta) ---
    if args.prand_theta:
        d = np.loadtxt(args.prand_theta)
        prand_t = np.interp(t_c, d[:, 0], d[:, 1])
    else:
        prand_t = np.sin(np.radians(t_c))
    prand_t /= (prand_t.sum() * args.dtheta + 1e-300)
    Prand = np.outer(prand_r, prand_t)

    # --- optional F_avg bulk correction (PMF column only) ---
    F_avg = 1.0
    if args.rdf:
        g = np.loadtxt(args.rdf)
        gr = np.interp(r_c, g[:, 0], g[:, 1])
        p_sum = np.trapz(P, t_c, axis=1)          # ~ P_tilde(r)
        with np.errstate(divide="ignore", invalid="ignore"):
            F = gr * prand_r / p_sum
        F_avg = np.nanmean(F[np.isfinite(F)])
        print(f"F_avg (bulk correction) = {F_avg:.4f}")

    RT = R_KJ * args.temperature
    with np.errstate(divide="ignore", invalid="ignore"):
        PMF = -RT * np.log(P * F_avg / (Prand + 1e-300))
    # Empty bins leave a literal +inf in the PMF column, which numpy round-trips but
    # text/Fortran consumers of the grid can fail to parse. Replace those with the
    # largest finite (least favourable) PMF present; dG_PMF = min(PMF) is unaffected.
    finite = np.isfinite(PMF)
    fill = PMF[finite].max() if finite.any() else 0.0
    PMF = np.where(finite, PMF, fill)
    # P*PMF integrand: 0 in empty observed bins (lim x*ln x = 0). PMF is now finite,
    # so this is a plain product with no 0*inf.
    PxPMF = np.where(P > 0.0, P * PMF, 0.0)

    # --- write full grid and first-shell grid ---
    def write(path, rmask=None):
        with open(path, "w") as fh:
            for i, r in enumerate(r_c):
                if rmask is not None and not rmask[i]:
                    continue
                for j, th in enumerate(t_c):
                    fh.write(f"{r:22.16E} {th:22.16E} {P[i,j]:22.16E} "
                             f"{Prand[i,j]:22.16E} {PMF[i,j]:22.16E} "
                             f"{PxPMF[i,j]:22.16E}\n")
    write(args.out)
    shell_mask = (r_c >= args.shell[0]) & (r_c <= args.shell[1])
    # Insert the _1st_shell tag before the extension of the *filename* only, so a
    # path whose directory contains '.dat' (e.g. a run_v2.data/ dir) is untouched.
    base, ext = os.path.splitext(args.out)
    shell_out = f"{base}_1st_shell{ext}" if ext else args.out + "_1st_shell"
    write(shell_out, shell_mask)
    print(f"wrote {args.out}\nwrote {shell_out}")
    print("-> compute the strength with:  python hb_strength.py " + shell_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
