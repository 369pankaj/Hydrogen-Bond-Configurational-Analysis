#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hb_strength.py  --  Hydrogen-bond strength (free energy) from a P(r,theta) grid.

Implements the reversible-work / PMF methodology of

    Sapir & Harries, "Revisiting Hydrogen Bond Thermodynamics in Molecular
    Simulations", J. Chem. Theory Comput. 2017, 13, 2851-2857
    (DOI 10.1021/acs.jctc.7b00238)

and the internal note by G. Olgenblum (2023) that this project followed.

Given the joint probability density P(r,theta) of a donor-acceptor pair and the
random reference density P_rand(r,theta) = P_rand(r).P_rand(theta), the H-bond
free energy is the reversible work of forging the observed distribution from the
random one, integrated over the first solvation shell (Eq. 2 of the paper):

    dG = -RT * INT INT  zeta*P(r,th) * ln[ zeta*P(r,th) / (zeta_rand*P_rand(r,th)) ] dr dth

with the shell-normalisers

    1/zeta      = INT INT P(r,th)      dr dth   (over [rmin,rmax] x [0,pi])
    1/zeta_rand = INT INT P_rand(r,th) dr dth

A global multiplicative bulk-correction on P (the F_avg / g(r) factor applied when
building the grid) CANCELS in this expression because zeta re-normalises P over the
shell -- so dG is identical whether or not that correction was applied to the grid.

This script is the *validated* core of the pipeline: run `validate.py` to confirm it
reproduces the reference 4 M urea/water values to 4 decimals.

INPUT GRID FORMAT (whitespace-separated columns, r-major ordering i.e. theta varies
fastest):
    col1 = r (Angstrom)      col2 = theta (degrees)
    col3 = P(r,theta)        col4 = P_rand(r,theta)
    col5 = PMF(r,theta)      [optional; = -RT ln(P.F_avg/P_rand), used for dG_PMF]
This is exactly the `Free_energy_*_1st_shell.dat` layout written by the Fortran
front-end and by hb_pmf_traj.py.

USAGE
    python hb_strength.py FILE [FILE ...]                 # table for one/many grids
    python hb_strength.py Free_energy_WW_1st_shell.dat
    python hb_strength.py Free_energy_*_1st_shell.dat --csv strengths.csv
    python hb_strength.py FULL_grid.dat --first-shell 2.49 3.25   # slice shell first
    python hb_strength.py FILE --plot pmf.png            # PMF(r,theta) landscape
    python hb_strength.py FILE -T 298                    # different temperature
"""
import argparse
import os
import sys
import numpy as np

R_KJ = 0.0083144621   # kJ / (mol K)  -- matches the Fortran / MATLAB constant


# --------------------------------------------------------------------------- #
#  Grid I/O
# --------------------------------------------------------------------------- #
def load_grid(path, first_shell=None):
    """Load a P(r,theta) grid file and reshape to (n_r, n_theta).

    Returns a dict with r-axis, theta-axis and 2-D arrays P, Prand, and PMF
    (PMF is None if the file has < 5 columns).
    """
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] < 4:
        raise ValueError(f"{path}: expected >=4 columns (r theta P Prand), "
                         f"got shape {raw.shape}")
    if first_shell is not None:
        lo, hi = first_shell
        raw = raw[(raw[:, 0] >= lo) & (raw[:, 0] <= hi)]
        if raw.size == 0:
            raise ValueError(f"{path}: no rows in first shell [{lo},{hi}]")

    n_theta = np.unique(np.round(raw[:, 1], 4)).size
    if raw.shape[0] % n_theta:
        raise ValueError(f"{path}: {raw.shape[0]} rows not divisible by "
                         f"{n_theta} angles -- is the grid r-major?")
    n_r = raw.shape[0] // n_theta
    shp = (n_r, n_theta)
    # Guard against a theta-major file (same total row count, same divisibility)
    # being silently transposed. In the required r-major layout theta varies
    # fastest, so r is constant across each reshaped row and theta is constant
    # down each column. Verify that before trusting the reshape.
    R = raw[:, 0].reshape(shp)
    T = raw[:, 1].reshape(shp)
    if not (np.allclose(R, R[:, :1]) and np.allclose(T, T[:1, :])):
        raise ValueError(f"{path}: grid is not r-major (theta varying fastest) "
                         f"as required -- refusing to reshape (would transpose)")
    g = {
        "r":     raw[:, 0].reshape(shp)[:, 0],   # 1-D r axis
        "theta": raw[:, 1].reshape(shp)[0, :],   # 1-D theta axis (deg)
        "P":     raw[:, 2].reshape(shp),
        "Prand": raw[:, 3].reshape(shp),
        "PMF":   raw[:, 4].reshape(shp) if raw.shape[1] >= 5 else None,
        "path":  path,
    }
    return g


def _int2d(f, r, th):
    """Double trapezoidal integral of f(r,theta) over the grid."""
    return np.trapz(np.trapz(f, th, axis=1), r)


# --------------------------------------------------------------------------- #
#  Free-energy definitions
# --------------------------------------------------------------------------- #
def hbond_strength(grid, temperature=300.0):
    """Return a dict of H-bond free-energy metrics (kJ/mol) for one grid."""
    RT = R_KJ * temperature
    r, th, P, Prand = grid["r"], grid["theta"], grid["P"], grid["Prand"]

    zeta = 1.0 / _int2d(P, r, th)
    zeta_rand = 1.0 / _int2d(Prand, r, th)
    N2 = zeta * P            # shell-normalised observed density
    N1 = zeta_rand * Prand   # shell-normalised random density

    # Eq. 2, Sapir & Harries 2017  (the reported convention).
    # Both densities must be positive for the reversible-work integrand to be
    # finite. Empty observed bins (N2 == 0) vanish by lim_{x->0} x*ln x = 0;
    # empty random bins (N1 == 0) come from a zero-padded reference, not a real
    # divergence, so drop them too (a lower bound, like the sparse-grid case)
    # rather than returning +inf/NaN from log(x/0) or 0*log(0). Fully-populated
    # grids (the committed Fortran grids) have both > 0 everywhere, so their
    # reported values are bit-identical.
    both = (N2 > 0.0) & (N1 > 0.0)
    if not both.all():
        print(f"warning: {int((~both).sum())} bin(s) dropped "
              f"(empty observed or reference density) -- dG is a lower bound")
    N1s = np.where(N1 > 0.0, N1, 1.0)
    N2s = np.where(N2 > 0.0, N2, 1.0)
    integ_2017 = np.where(both, N2 * np.log(N2 / N1s), 0.0)
    dG_2017 = -RT * _int2d(integ_2017, r, th)
    # Modified convention, Olgenblum et al. JCTC 2020, 16, 1249-1262 (reverse KL).
    # Same both>0 mask; reverse KL diverges where N2 == 0 and is 0*(-inf)=NaN
    # where N1 == 0, so masking both keeps it finite. Fully-populated grids
    # are unaffected.
    integ_2020 = np.where(both, N1 * np.log(N1 / N2s), 0.0)
    dG_2020 = -RT * _int2d(integ_2020, r, th)

    # Eq. 1 global minimum (deepest PMF wrt bulk); needs the bulk-corrected PMF col.
    dG_pmf = float(np.min(grid["PMF"])) if grid["PMF"] is not None else np.nan

    return {"dG": dG_2017, "dG_2020": dG_2020, "dG_PMF": dG_pmf,
            "RT": RT, "n_r": r.size, "n_theta": th.size,
            "r_range": (float(r.min()), float(r.max()))}


# --------------------------------------------------------------------------- #
#  Plotting (optional)
# --------------------------------------------------------------------------- #
def plot_pmf(grid, out_png, temperature=300.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    RT = R_KJ * temperature
    if grid["PMF"] is not None:
        Z = grid["PMF"]
    else:
        Z = -RT * np.log(grid["P"] / grid["Prand"])
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    Th, Rr = np.meshgrid(grid["theta"], grid["r"])
    vlim = np.nanpercentile(np.abs(Z), 98)
    pc = ax.pcolormesh(Th, Rr, Z, cmap="jet", shading="auto",
                       vmin=-vlim, vmax=vlim)
    ax.set_xlabel(r"$\theta$ (degrees)")
    ax.set_ylabel(r"$r$ ($\AA$)")
    ax.set_title(os.path.basename(grid["path"]))
    fig.colorbar(pc, label=r"PMF$_{\rm HB}(r,\theta)$  [kJ/mol]")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"  wrote {out_png}")


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="H-bond strength (free energy) from P(r,theta) grids "
                    "[Sapir & Harries, JCTC 2017, 13, 2851].",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("files", nargs="+", help="Free_energy_*_1st_shell.dat grid file(s)")
    ap.add_argument("-T", "--temperature", type=float, default=300.0,
                    help="temperature in K")
    ap.add_argument("--first-shell", nargs=2, type=float, metavar=("RMIN", "RMAX"),
                    default=None, help="restrict a FULL grid to this r-shell first")
    ap.add_argument("--csv", default=None, help="also write results to this CSV")
    ap.add_argument("--plot", default=None,
                    help="write PMF landscape PNG (single input file only)")
    ap.add_argument("--all", action="store_true",
                    help="also print the 2020 (reverse-KL) and dG_PMF columns")
    args = ap.parse_args(argv)

    rows = []
    for path in args.files:
        grid = load_grid(path, first_shell=args.first_shell)
        res = hbond_strength(grid, temperature=args.temperature)
        res["file"] = os.path.basename(path)
        rows.append(res)
        if args.plot:
            plot_pmf(grid, args.plot, temperature=args.temperature)

    # pretty table
    name_w = max(len(r["file"]) for r in rows)
    if args.all:
        hdr = f"{'grid':{name_w}s}  {'dG_2017':>9s}  {'dG_2020':>9s}  {'dG_PMF':>9s}  r-shell"
    else:
        hdr = f"{'grid':{name_w}s}  {'dG (kJ/mol)':>11s}  r-shell"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        rr = f"[{r['r_range'][0]:.2f},{r['r_range'][1]:.2f}]"
        if args.all:
            print(f"{r['file']:{name_w}s}  {r['dG']:9.4f}  {r['dG_2020']:9.4f}  "
                  f"{r['dG_PMF']:9.4f}  {rr}")
        else:
            print(f"{r['file']:{name_w}s}  {r['dG']:11.4f}  {rr}")
    print(f"\n(T = {args.temperature:.1f} K, RT = {R_KJ*args.temperature:.5f} kJ/mol; "
          f"dG = Sapir-Harries 2017 Eq. 2)")

    if args.csv:
        import csv
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["grid", "dG_2017_kJmol", "dG_2020_kJmol", "dG_PMF_kJmol",
                        "r_min", "r_max", "T_K"])
            for r in rows:
                w.writerow([r["file"], f"{r['dG']:.4f}", f"{r['dG_2020']:.4f}",
                            f"{r['dG_PMF']:.4f}", f"{r['r_range'][0]:.3f}",
                            f"{r['r_range'][1]:.3f}", f"{args.temperature:.1f}"])
        print(f"wrote {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
