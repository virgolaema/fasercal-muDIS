#!/bin/bash
# scripts/run_toy.sh — one-shot driver for the FASERcal charm-tagging toy.
#
#   1. shower_charm.py : re-shower the muon-DIS POWHEG LHE with Pythia8 and cache
#                        per-charm-hadron truth (sign, flight length, decay-muon
#                        kinematics) to the .npz  [HEAVY: needs pythia8 / LCG view]
#   2. make_report.py  : apply the toy FASERcal+spectrometer response and build
#                        the PDF  [FAST: pure numpy, iterate freely]
#
# Reuses the FASERnu POWHEG productions (production_v1) on EOS.  Paths come from
# json/ev.json (see json/settings_template.json).
#
# Usage:  ./scripts/run_toy.sh [--nevents N] [--report-only]

set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LCG_VIEW="${LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt}"
NEVENTS=20000
SEEDS="1-10"
REPORT_ONLY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --nevents) NEVENTS="$2"; shift 2 ;;
        --seeds) SEEDS="$2"; shift 2 ;;
        --report-only) REPORT_ONLY=1; shift ;;
        -h|--help) sed -n '2,14p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

set +u; source "$LCG_VIEW/setup.sh"; set -u

if [[ "$REPORT_ONLY" -eq 0 ]]; then
    echo "=== run_toy.sh: showering charm muons (production_v1, seeds $SEEDS) ==="
    python3 "$REPO/python/shower_charm.py" --nevents "$NEVENTS" --seeds "$SEEDS" \
        2>&1 | tee "$REPO/logs/shower_charm.log"
fi

echo "=== run_toy.sh: building report ==="
python3 "$REPO/python/make_report.py" 2>&1 | tee "$REPO/logs/make_report.log"
echo "=== run_toy.sh done ==="
