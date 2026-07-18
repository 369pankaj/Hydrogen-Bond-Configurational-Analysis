#!/usr/bin/env python3
"""
validate.py -- regression test for the H-bond strength pipeline.

Confirms that hb_strength.py reproduces the reference 4 M urea/water H-bond
free energies (Smith_KWB, 1st CPI paper) from the committed grid files, and
cross-checks the water-water value against Sapir & Harries (JCTC 2017, 13, 2851).

Run:   python validate.py            (uses grids in ../ by default)
       python validate.py /path/to/grids
Exit code 0 = all values reproduced within tolerance.
"""
import os
import sys
import numpy as np
from hb_strength import load_grid, hbond_strength

TOL = 1.0e-3  # kJ/mol

# canonical grid file  ->  (label, reported dG kJ/mol)
REFERENCE = [
    ("Free_energy_WW_1st_shell.dat",         "W-W  (water-water)",       -3.4941),
    ("Free_energy_U-W_1st_shell.dat",        "U-W  (urea-water, all)",   -1.3816),
    ("Free_energy_CO-W_nf_uw_1st_shell.dat", "Ou-W (carbonyl O-water)",  -4.1218),
    ("Free_energy_NH-W_nf_uw_1st_shell.dat", "Nu-W (amine N-water)",     -1.1155),
]
# Literature cross-check: Sapir & Harries 2017, Table 1, WW dG (Eq. 2)
PAPER_WW = {"pure water": -3.7, "0.8 M trehalose": -3.8}


def main():
    gdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")
    print(f"grids in: {gdir}\n")
    print(f"{'species':26s} {'computed':>10s} {'reported':>10s} {'|diff|':>9s}  ok")
    print("-" * 70)
    ok_all = True
    for fname, label, ref in REFERENCE:
        path = os.path.join(gdir, fname)
        if not os.path.exists(path):
            print(f"{label:26s}  MISSING FILE: {fname}")
            ok_all = False
            continue
        val = hbond_strength(load_grid(path))["dG"]
        ok = abs(val - ref) < TOL
        ok_all &= ok
        print(f"{label:26s} {val:10.4f} {ref:10.4f} {abs(val-ref):9.4f}  "
              f"{'PASS' if ok else 'FAIL'}")

    print("\nLiterature cross-check (Sapir & Harries 2017, Table 1, WW dG, Eq. 2):")
    ww = hbond_strength(load_grid(
        os.path.join(gdir, "Free_energy_WW_1st_shell.dat")))["dG"]
    for sysname, pval in PAPER_WW.items():
        print(f"  this work (4 M urea): {ww:6.3f} kJ/mol   vs  {sysname}: "
              f"{pval:+.1f} kJ/mol   (delta {ww-pval:+.2f})")

    print("\n" + ("ALL VALUES REPRODUCED (PASS)" if ok_all
                  else "*** SOME VALUES FAILED ***"))
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
