# config.sh -- systems, species and run defaults for the H-bond strength pipeline.
#
# This file is *sourced* by run_hb.sh and submit.slurm. Edit the SYSTEMS and
# SPECIES tables to add new trajectories or donor/acceptor pairs; edit the
# defaults (or override them from the environment, e.g. `STRIDE=100 ./run_hb.sh`).
#
# Nothing here is run on its own -- there are no side effects beyond defining
# variables, so it is safe to `source config.sh` from an interactive shell too.

# Absolute path to the directory holding these scripts (works no matter where
# the driver is invoked from).
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --------------------------------------------------------------------------- #
#  Run defaults  (every one is overridable from the environment)
# --------------------------------------------------------------------------- #
# Python interpreter with MDAnalysis + numpy available.
PYTHON="${PYTHON:-python}"

# Directory holding the reference random-angle files  P_rand_*.dat
# (the Fortran front-end wrote these one level up, in Analysis/hbond_energy).
PRAND_DIR="${PRAND_DIR:-$PIPELINE_DIR/..}"

# Where per-system grids, logs and CSV summaries are written.
RESULTS_DIR="${RESULTS_DIR:-$PIPELINE_DIR/results}"

# Frame stride. 1 = every frame; 200-500 is plenty for converged strengths and
# is what the docstring examples use. Smaller stride -> tighter statistics.
STRIDE="${STRIDE:-200}"

# First-shell integration window [Angstrom] over which dG (Eq. 2) is integrated.
SHELL_MIN="${SHELL_MIN:-2.49}"
SHELL_MAX="${SHELL_MAX:-3.25}"

# Temperature [K]. Used for BOTH the PMF column written by hb_pmf_traj.py and the
# dG reported by hb_strength.py, so they always agree (keep them equal).
TEMPERATURE="${TEMPERATURE:-300}"

# --------------------------------------------------------------------------- #
#  Systems:  name -> "topology|trajectory"
#  Add a line per production run you want analysed. Paths may be absolute or
#  relative to wherever you launch the driver.
# --------------------------------------------------------------------------- #
declare -gA SYSTEMS=(
  [4M]="/work/Urea_Water/Smith_KWB/4M/md.gro|/work/Urea_Water/Smith_KWB/4M/md.xtc"
  [4M_500]="/work/Urea_Water/Smith_KWB/4M/4M_500/md.gro|/work/Urea_Water/Smith_KWB/4M/4M_500/md.xtc"
  [4M_1000]="/work/Urea_Water/Smith_KWB/4M/4M_1000/md.gro|/work/Urea_Water/Smith_KWB/4M/4M_1000/md.xtc"
)
# Order systems are processed / reported in (associative arrays are unordered).
SYSTEM_ORDER=(4M 4M_500 4M_1000)

# --------------------------------------------------------------------------- #
#  Species:  name -> "donor_heavy | acceptor | prand_file | extra_flags"
#
#   donor_heavy  : MDAnalysis selection for the donor HEAVY atom(s)
#   acceptor     : MDAnalysis selection for the acceptor heavy atom(s)
#   prand_file   : random-angle reference (looked up in $PRAND_DIR)
#   extra_flags  : passed verbatim to hb_pmf_traj.py (may be empty)
#
#  WW  is symmetric (both partners' H's are used) and excludes self-pairs.
#  OuW / NuW are directional, so no --symmetric and no --exclude-same-residue.
# --------------------------------------------------------------------------- #
declare -gA SPECIES=(
  [WW]="resname SOL and name OW|resname SOL and name OW|P_rand_WW.dat|--symmetric --exclude-same-residue"
  [OuW]="resname SOL and name OW|resname UREA and name OU|P_rand_WH.dat|"
  [NuW]="resname UREA and name N1U N2U|resname SOL and name OW|P_rand_WW.dat|"
)
# Order species are run / reported in.
SPECIES_ORDER=(WW OuW NuW)
