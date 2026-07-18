#!/usr/bin/env bash
# run_hb.sh -- build P(r,theta) grids and compute H-bond strengths for one or
#              more systems, driving hb_pmf_traj.py + hb_strength.py.
#
# USAGE
#   ./run_hb.sh                 # every system in SYSTEM_ORDER (see config.sh)
#   ./run_hb.sh all             # same as above
#   ./run_hb.sh 4M              # just the 4M system
#   ./run_hb.sh 4M 4M_500       # a subset
#   STRIDE=100 ./run_hb.sh 4M   # override any default from config.sh via env
#
# Per system it writes, under $RESULTS_DIR/<system>/ :
#   <SP>_grid.dat              full 6-column grid            (per species SP)
#   <SP>_grid_1st_shell.dat    first-shell slice used for dG
#   strengths_<system>.csv     dG table (2017 / 2020 / PMF)
#   run.log                    full stdout/stderr of the run
#
# Exit status is non-zero if any system fails.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh"

die() { echo "run_hb.sh: $*" >&2; exit 1; }

run_one_system() {
  local sysname="$1"
  local spec="${SYSTEMS[$sysname]:-}"
  [[ -n "$spec" ]] || die "unknown system '$sysname' (known: ${SYSTEM_ORDER[*]})"

  local top="${spec%%|*}" traj="${spec##*|}"
  [[ -r "$top"  ]] || die "[$sysname] topology not readable: $top"
  [[ -r "$traj" ]] || die "[$sysname] trajectory not readable: $traj"

  local outdir="$RESULTS_DIR/$sysname"
  mkdir -p "$outdir"
  local log="$outdir/run.log"
  : > "$log"

  echo "========================================================================"
  echo " system   : $sysname"
  echo " topology : $top"
  echo " traj     : $traj"
  echo " stride   : $STRIDE      shell: [$SHELL_MIN, $SHELL_MAX] A      T: $TEMPERATURE K"
  echo " outdir   : $outdir"
  echo "========================================================================"

  local shells=() sp donor acc prand extra grid
  for sp in "${SPECIES_ORDER[@]}"; do
    IFS='|' read -r donor acc prand extra <<< "${SPECIES[$sp]}"
    local prand_path="$PRAND_DIR/$prand"
    [[ -r "$prand_path" ]] || die "[$sysname/$sp] P_rand file not readable: $prand_path"
    grid="$outdir/${sp}_grid.dat"

    echo ">>> [$sysname/$sp] donor='$donor'  acceptor='$acc'  prand='$prand'  ${extra:-(directional)}"
    # $extra is intentionally unquoted so multi-word flag strings split into args.
    # shellcheck disable=SC2086
    "$PYTHON" "$HERE/hb_pmf_traj.py" \
        --top "$top" --traj "$traj" \
        --donor-heavy "$donor" --acceptor "$acc" \
        --prand-theta "$prand_path" \
        --temperature "$TEMPERATURE" \
        --shell "$SHELL_MIN" "$SHELL_MAX" \
        --stride "$STRIDE" \
        $extra \
        -o "$grid" 2>&1 | tee -a "$log"

    shells+=("$outdir/${sp}_grid_1st_shell.dat")
  done

  echo ">>> [$sysname] H-bond strengths:"
  "$PYTHON" "$HERE/hb_strength.py" "${shells[@]}" \
      -T "$TEMPERATURE" --all \
      --csv "$outdir/strengths_${sysname}.csv" 2>&1 | tee -a "$log"
  echo ">>> [$sysname] done -> $outdir/strengths_${sysname}.csv"
  echo
}

main() {
  local -a targets=()
  if [[ $# -eq 0 || ( $# -eq 1 && "${1:-}" == "all" ) ]]; then
    targets=("${SYSTEM_ORDER[@]}")
  else
    targets=("$@")
  fi

  local sysname failed=0
  for sysname in "${targets[@]}"; do
    if ! run_one_system "$sysname"; then
      echo "run_hb.sh: system '$sysname' FAILED" >&2
      failed=1
    fi
  done

  echo "========================================================================"
  if [[ $failed -eq 0 ]]; then
    echo "ALL SYSTEMS COMPLETED.  Summaries:"
    for sysname in "${targets[@]}"; do
      echo "  $RESULTS_DIR/$sysname/strengths_${sysname}.csv"
    done
  else
    echo "*** ONE OR MORE SYSTEMS FAILED -- see logs above ***"
  fi
  return $failed
}

main "$@"
