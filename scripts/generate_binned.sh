#!/bin/bash
# scripts/generate_binned.sh — regenerate the muon-DIS samples in FIXED-ENERGY
# BINS instead of with the flux as a beam PDF.
#
# WHY.  The original production handed the muon flux to POWHEG as the "PDF" of a
# fictitious 7 TeV muon beam.  That beam has a REMNANT which hadronises and
# contaminates the hadronic final state: energy conservation failed in 80% of
# events and E_had/nu reached 2.5 at large x, which forced the analysis to take
# the hadronic four-vector from momentum conservation rather than from the
# particles (see docs/XRECO.md).  Generating at fixed beam energy with
# `fixed_lepton_beam 1` gives x_lepton = 1, so there is no lepton-side remnant at
# all, and E_mu is exactly the bin energy.  The flux is applied afterwards as a
# per-bin weight (python/flux_bins.py).
#
# Usage:
#   ./scripts/generate_binned.sh --outdir /eos/.../production_binned \
#        [--bins 0-19] [--samples mum_proton,...] [--numevts 20000] [--ncores N]
#
# Produces  <outdir>/<sample>/bin<NN>/pwgevents-0001.lhe  plus a bins.json
# recording the energy and flux weight of every bin.

set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GEN="/afs/cern.ch/work/e/evilla/private/faser/muondisgenerator"

OUTDIR=""; BINS="0-19"; NUMEVTS=20000
SAMPLES="mum_proton,mum_neutron,mup_proton,mup_neutron"
NCORES="$(nproc 2>/dev/null || echo 8)"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)  OUTDIR="$2"; shift 2 ;;
        --bins)    BINS="$2"; shift 2 ;;
        --samples) SAMPLES="$2"; shift 2 ;;
        --numevts) NUMEVTS="$2"; shift 2 ;;
        --ncores)  NCORES="$2"; shift 2 ;;
        -h|--help) sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 1; }

# per-sample POWHEG beam/target codes (as in the generator's own scripts)
declare -A IH1=( [mum_proton]=13  [mum_neutron]=13  [mup_proton]=-13 [mup_neutron]=-13 )
declare -A IH2=( [mum_proton]=1   [mum_neutron]=2   [mup_proton]=1   [mup_neutron]=2   )
declare -A EB2=( [mum_proton]="0.938d0" [mum_neutron]="0.940d0" \
                 [mup_proton]="0.938d0" [mup_neutron]="0.940d0" )

source "$GEN/scripts/env.sh" >/dev/null 2>&1

# bin energies from the same code the reweighting uses, so they cannot drift
mapfile -t EGEN < <(cd "$REPO" && python3 -c "
import sys; sys.path.insert(0,'python')
from flux_bins import flux_weights
for e in flux_weights()['e_centre']: print('%.4f'%e)")

IFS='-' read -r B0 B1 <<< "${BINS/,/-}"
mkdir -p "$OUTDIR"
(cd "$REPO" && python3 -c "
import json,sys; sys.path.insert(0,'python')
from flux_bins import flux_weights
b=flux_weights()
json.dump({k:(v.tolist() if hasattr(v,'tolist') else v) for k,v in b.items()},
          open('$OUTDIR/bins.json','w'), indent=1)")

run_one() {
    local sample=$1 ibin=$2
    # EGEN is a bash array in the parent; xargs spawns a fresh shell, so it is
    # passed through as a space-separated string and rebuilt here.
    local -a EG=( $EGEN_STR )
    local e="${EG[$ibin]}"
    local -A IH1 IH2 EB2
    local -a _k _v
    IFS='|' read -r _ks _vs <<< "$IH1_STR"; _k=( $_ks ); _v=( $_vs )
    for i in "${!_k[@]}"; do IH1[${_k[$i]}]=${_v[$i]}; done
    IFS='|' read -r _ks _vs <<< "$IH2_STR"; _k=( $_ks ); _v=( $_vs )
    for i in "${!_k[@]}"; do IH2[${_k[$i]}]=${_v[$i]}; done
    IFS='|' read -r _ks _vs <<< "$EB2_STR"; _k=( $_ks ); _v=( $_vs )
    for i in "${!_k[@]}"; do EB2[${_k[$i]}]=${_v[$i]}; done
    local d; d="$(printf '%s/%s/bin%02d' "$OUTDIR" "$sample" "$ibin")"
    [[ -s "$d/pwgevents-0001.lhe" ]] && { echo "  skip $sample bin$ibin (done)"; return 0; }
    mkdir -p "$d"; cd "$d"
    source "$GEN/scripts/env.sh" >/dev/null 2>&1
    sed -e "s/^ih1 .*/ih1 ${IH1[$sample]}/" \
        -e "s/^ih2 .*/ih2 ${IH2[$sample]}/" \
        -e "s/^ebeam1 .*/ebeam1 ${e}d0/" \
        -e "s/^ebeam2 .*/ebeam2 ${EB2[$sample]}/" \
        -e "s/^fixed_lepton_beam 0/fixed_lepton_beam 1/" \
        -e "s/^numevts .*/numevts $NUMEVTS/" \
        -e "s/^LEPpdf /!LEPpdf /" \
        "$GEN/config/powheg_original.input" > powheg_original.input
    cp "$GEN/config/pwgseeds.dat" .
    for ig in 1 2 3; do
        sed "s/xgriditeration.*/xgriditeration $ig/; s/parallelstage.*/parallelstage 1/" \
            powheg_original.input > powheg.input
        echo 1 | pwhg_main > "s1_$ig.log" 2>&1
    done
    for st in 2 3 4; do
        sed "s/parallelstage.*/parallelstage $st/" powheg_original.input > powheg.input
        echo 1 | pwhg_main > "s$st.log" 2>&1
    done
    echo "  done $sample bin$ibin  E=$e GeV  events=$(grep -c '<event>' pwgevents-0001.lhe 2>/dev/null || echo 0)"
}
export -f run_one
export OUTDIR NUMEVTS GEN
export EGEN_STR="${EGEN[*]}"
# the per-sample maps must reach the xargs subshells too
export IH1_STR="${!IH1[*]}|${IH1[*]}" IH2_STR="${!IH2[*]}|${IH2[*]}" EB2_STR="${!EB2[*]}|${EB2[*]}"

echo "=== binned generation: bins $B0-$B1, samples $SAMPLES, $NUMEVTS evts/bin ==="
for sample in ${SAMPLES//,/ }; do
    for ibin in $(seq "$B0" "$B1"); do
        echo "$sample $ibin"
    done
done | xargs -P "$NCORES" -n 2 bash -c 'run_one "$0" "$1"'
echo "=== done -> $OUTDIR ==="
