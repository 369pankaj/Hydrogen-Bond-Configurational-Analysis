#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hb_run.py -- automation front-end for the H-bond strength pipeline.

Reads a chemistry registry (molecules.yaml) and a job list (systems.yaml), then
for every system auto-generates each valid donor->acceptor species and drives
hb_pmf_traj.py + hb_strength.py. This is what makes the pipeline reusable across
binary osmolyte-water mixtures, ternary mixtures, and protein-containing systems
without hand-writing selections per pair.

Species generation rules (per system's `components` list):
  * a molecule with no donors is acceptor-only  -> only ever on the acceptor side
    (e.g. TMAO/DMSO/betaine: water donates TO them);
  * a like-like self-pair flagged `like_like_symmetric` (water) -> one SYMMETRIC
    species using both partners' H's;
  * any other same-species pair (urea-urea) -> one directional self species
    (donor atoms -> acceptor atoms), self-pairs removed;
  * cross-species pairs -> one directional species per (donor comp, acceptor comp)
    ordering, so both urea->water and water->urea appear.

Self / intramolecular exclusion is by residue for small molecules and by bonded
fragment ("molecule") for proteins (needs a .tpr/.psf); override per system with
`exclude: molecule|residue`.

USAGE
  ./hb_run.py                     # run every (uncommented) system in systems.yaml
  ./hb_run.py 4M 4M_tpr           # only these systems
  ./hb_run.py --dry-run 4M        # print the planned species; run nothing
  ./hb_run.py --stride 100 4M     # override a default for this invocation
"""
import argparse
import csv
import os
import subprocess
import sys

try:
    import yaml
except ImportError:
    sys.exit("hb_run.py needs PyYAML:  pip install pyyaml")

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------- #
#  Registry -> selections
# --------------------------------------------------------------------------- #
def donor_sel(comp, mols):
    """MDAnalysis donor-heavy selection for a component, or '' if it never donates."""
    m = mols[comp]
    if m.get("kind") == "protein":
        return m.get("donor_sel", "").strip()
    names = m.get("donors") or []
    return f"resname {m['resname']} and name {' '.join(names)}" if names else ""


def acceptor_sel(comp, mols):
    """MDAnalysis acceptor selection for a component, or '' if it never accepts."""
    m = mols[comp]
    if m.get("kind") == "protein":
        return m.get("acceptor_sel", "").strip()
    names = m.get("acceptors") or []
    return f"resname {m['resname']} and name {' '.join(names)}" if names else ""


def exclude_level(comp, mols, override):
    """How to remove intramolecular pairs for a same-species species."""
    if override:
        return override
    return "molecule" if mols[comp].get("kind") == "protein" else "residue"


class Species:
    __slots__ = ("name", "donor", "acceptor", "symmetric", "exclude")

    def __init__(self, name, donor, acceptor, symmetric, exclude):
        self.name = name
        self.donor = donor
        self.acceptor = acceptor
        self.symmetric = symmetric     # bool
        self.exclude = exclude         # 'residue' | 'molecule' | None


def generate_species(components, mols, exclude_override):
    """Enumerate the donor->acceptor species for one system's component list."""
    for cd in components:
        dsel = donor_sel(cd, mols)
        if not dsel:
            continue
        for ca in components:
            asel = acceptor_sel(ca, mols)
            if not asel:
                continue
            if cd == ca:
                lvl = exclude_level(cd, mols, exclude_override)
                if mols[cd].get("like_like_symmetric"):
                    yield Species(f"{cd}_sym", dsel, asel, True, lvl)
                else:
                    # e.g. urea->urea : donor N, acceptor O, drop same-molecule
                    yield Species(f"{cd}_self", dsel, asel, False, lvl)
            else:
                # cross species are always inter-molecular -> no self exclusion
                yield Species(f"{cd}D_{ca}A", dsel, asel, False, None)


# --------------------------------------------------------------------------- #
#  Optional first-shell auto-detection from the pair RDF
# --------------------------------------------------------------------------- #
def auto_shell(top, traj, donor, acceptor, rmax, stride):
    """First-shell [rmin, rmax] from the donor-acceptor RDF: rise of the first
    peak to the following minimum. Returns None on failure (caller falls back)."""
    try:
        import numpy as np
        import MDAnalysis as mda
        from MDAnalysis.analysis.rdf import InterRDF
        u = mda.Universe(top, traj)
        d = u.select_atoms(donor)
        a = u.select_atoms(acceptor)
        if len(d) == 0 or len(a) == 0:
            return None
        rdf = InterRDF(d, a, nbins=120, range=(1.8, rmax))
        rdf.run(step=max(stride, 1))
        r, g = rdf.results.bins, rdf.results.rdf
        peak = int(np.argmax(g))                       # first (dominant) peak
        # first local minimum after the peak
        hi = None
        for k in range(peak + 1, len(g) - 1):
            if g[k] <= g[k - 1] and g[k] <= g[k + 1]:
                hi = r[k]
                break
        if hi is None:
            return None
        # lower edge: where g rises through half its peak height before the peak
        half = 0.5 * g[peak]
        lo = r[0]
        for k in range(peak, 0, -1):
            if g[k] < half:
                lo = r[k]
                break
        return (round(float(lo), 3), round(float(hi), 3))
    except Exception as e:
        print(f"    auto-shell failed ({e.__class__.__name__}: {e}); using default")
        return None


# --------------------------------------------------------------------------- #
#  Driving the engine
# --------------------------------------------------------------------------- #
def build_cmd(sp, cfg, grid_path):
    cmd = [sys.executable, os.path.join(HERE, "hb_pmf_traj.py"),
           "--top", cfg["topology"], "--traj", cfg["trajectory"],
           "--donor-heavy", sp.donor, "--acceptor", sp.acceptor,
           "--temperature", str(cfg["temperature"]),
           "--stride", str(cfg["stride"]),
           "--rmax", str(cfg["rmax"]),
           "--shell", str(cfg["shell"][0]), str(cfg["shell"][1]),
           "-o", grid_path]
    if cfg.get("prand_theta"):
        cmd += ["--prand-theta", cfg["prand_theta"]]
    if sp.symmetric:
        cmd += ["--symmetric"]
    if sp.exclude == "residue":
        cmd += ["--exclude-same-residue"]
    elif sp.exclude == "molecule":
        cmd += ["--exclude-same-molecule"]
    if cfg.get("neighbor_list"):
        cmd += ["--neighbor-list"]
    return cmd


def run_system(name, syscfg, defaults, mols, results_dir, dry_run):
    cfg = dict(defaults)
    cfg.update({k: v for k, v in syscfg.items() if k != "components"})
    for key in ("topology", "trajectory"):
        if not cfg.get(key):
            print(f"[{name}] SKIP: missing '{key}'"); return True
        if not dry_run and not os.path.exists(cfg[key]):
            print(f"[{name}] SKIP: {key} not found: {cfg[key]}"); return True

    components = syscfg["components"]
    exclude_override = syscfg.get("exclude")
    species = list(generate_species(components, mols, exclude_override))

    print("=" * 78)
    print(f" system     : {name}")
    print(f" components : {', '.join(components)}")
    print(f" topology   : {cfg['topology']}")
    print(f" species    : {len(species)} -> {', '.join(s.name for s in species)}")
    print("=" * 78)

    outdir = os.path.join(results_dir, name)
    if not dry_run:
        os.makedirs(outdir, exist_ok=True)

    shells, ok = [], True
    for sp in species:
        # resolve 'auto' shell per species (each pair has its own geometry)
        if str(cfg.get("shell")) == "auto":
            det = None if dry_run else auto_shell(
                cfg["topology"], cfg["trajectory"], sp.donor, sp.acceptor,
                cfg["rmax"], cfg["stride"])
            spcfg = dict(cfg, shell=det or defaults.get("shell", [2.49, 3.25]))
        else:
            spcfg = cfg
        grid = os.path.join(outdir, f"{sp.name}_grid.dat")
        cmd = build_cmd(sp, spcfg, grid)

        excl = sp.exclude or "none"
        mode = "symmetric" if sp.symmetric else "directional"
        print(f">>> [{name}/{sp.name}] {mode}, exclude={excl}, "
              f"shell={spcfg['shell']}")
        print(f"    donor    = {sp.donor}")
        print(f"    acceptor = {sp.acceptor}")
        if dry_run:
            print("    CMD: " + " ".join(repr(c) if " " in c else c for c in cmd))
            continue
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"    FAILED (exit {e.returncode})"); ok = False; continue
        shells.append(os.path.join(outdir, f"{sp.name}_grid_1st_shell.dat"))

    if dry_run or not shells:
        return ok

    print(f">>> [{name}] strengths:")
    csv_path = os.path.join(outdir, f"strengths_{name}.csv")
    try:
        subprocess.run([sys.executable, os.path.join(HERE, "hb_strength.py"),
                        *shells, "-T", str(cfg["temperature"]),
                        "--all", "--csv", csv_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"    strengths FAILED (exit {e.returncode})"); ok = False
    else:
        print(f">>> [{name}] wrote {csv_path}")
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Automated H-bond strengths across binary/ternary/protein systems.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("systems", nargs="*", help="system names to run (default: all)")
    ap.add_argument("--molecules", default=os.path.join(HERE, "molecules.yaml"))
    ap.add_argument("--systems-file", default=os.path.join(HERE, "systems.yaml"))
    ap.add_argument("--results-dir", default=os.path.join(HERE, "results"))
    ap.add_argument("--dry-run", action="store_true",
                    help="print planned species + commands, execute nothing")
    ap.add_argument("--stride", type=int, help="override default stride")
    ap.add_argument("--temperature", type=float, help="override default temperature")
    args = ap.parse_args(argv)

    with open(args.molecules) as fh:
        mols = yaml.safe_load(fh)["molecules"]
    with open(args.systems_file) as fh:
        spec = yaml.safe_load(fh)
    defaults = spec.get("defaults", {})
    if args.stride is not None:
        defaults["stride"] = args.stride
    if args.temperature is not None:
        defaults["temperature"] = args.temperature
    all_systems = spec.get("systems", {})

    wanted = args.systems or list(all_systems)
    unknown = [s for s in wanted if s not in all_systems]
    if unknown:
        sys.exit(f"unknown system(s): {unknown}\nknown: {list(all_systems)}")

    failed = 0
    for name in wanted:
        if not run_system(name, all_systems[name], defaults, mols,
                          args.results_dir, args.dry_run):
            failed = 1
    print("=" * 78)
    print("DRY RUN complete." if args.dry_run else
          ("ALL SYSTEMS COMPLETED." if not failed else "*** SOME SYSTEMS FAILED ***"))
    return failed


if __name__ == "__main__":
    sys.exit(main())
