# muondis-fasercal

Toy Monte Carlo study of **charm vs anticharm identification in FASERcal**, a 3D
granular plastic-scintillator detector — a very different regime from the
FASERnu emulsion.  It reuses the muon-DIS POWHEG+Pythia8 inputs from the
[`generatoroutputanalysis`](../generatoroutputanalysis) contamination study.

## The physics question

FASERnu (emulsion) tags charm by seeing the sub-mm decay kink of the charmed
hadron.  FASERcal cannot: its cells are ~cm-scale and it has **no magnetic
field**.  The only robust charge-sign handle left is the **muon from the
semileptonic charm decay**:

```
c  -> s W+ , W+ -> mu+   (charm)        muon charge = +1
c~ -> s W- , W- -> mu-   (anticharm)    muon charge = -1
```

so the decay-muon charge *is* the charm sign.  But signing that muon needs a
**downstream magnetised spectrometer** (FASER2-like) — FASERcal only sees it
pass through MIP-like.  This toy quantifies how well that works.

## What the toy computes

- **Decay-muon lab momentum** (and hence its ionisation range) — it is "soft"
  only within the jet; in the lab it is tens of GeV and traverses the whole
  detector, reaching the spectrometer.
- **Charm-hadron flight length** vs scintillator cell size — the vertex handle
  that FASERnu uses and FASERcal largely lacks.
- **c/c~ tagging efficiency, purity, and asymmetry dilution** `(1-2η)` as a
  function of the spectrometer acceptance and charge-confusion assumptions.

## Layout

| path | role |
|------|------|
| `python/shower_charm.py` | Pythia8 re-showering → per-charm-hadron `.npz` cache (heavy, needs LCG view) |
| `python/detector.py`     | toy FASERcal + spectrometer response (all assumptions, tunable) |
| `python/make_report.py`  | `.npz` → PDF report (fast, pure numpy) |
| `python/config.py`       | path resolution via `json/ev.json` |
| `scripts/run_toy.sh`     | shower → report driver |

## Run

```bash
cp json/settings_template.json json/ev.json   # then edit paths
./scripts/run_toy.sh --nevents 20000          # full: shower + report
./scripts/run_toy.sh --report-only            # re-run just the report from cache
```

The `.npz` caches and PDFs live on EOS (see `json/ev.json`); the repo holds only
code.  The detector numbers in `detector.py` are **transparent toy defaults, not
FASERcal specifications** — the study is about how the result depends on them.
