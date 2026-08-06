# Binned generation: why, and what it changed

*Written 2026-08-06. Supersedes the generation strategy described in
[XRECO.md](XRECO.md) §"the hadronic energy problem".*

## The old strategy and its cost

The original samples handed the muon flux to POWHEG as the "beam PDF" of a
fictitious 7 TeV muon beam (arXiv:2506.13889 Eq. 2.1). It is elegant — one run
covers the whole spectrum — but that fictitious beam has a **remnant**, and the
remnant hadronises into the event.

The damage was worse than a few spurious GeV of activity:

* energy conservation failed in **80 %** of events;
* `E_had / nu` reached **2.5** at large `x`, exactly where `nu` is small and the
  measurement matters;
* no particle-level cut removed it. Beam-remnant status codes (61–63), a 1 mrad
  forward cone, and mother-chain ancestry to beam 2 were all tried and all
  failed — the last because Pythia's colour reconnection leaves the ancestry
  graph too connected to separate the two beam sides.

That forced the hadronic four-vector to be taken from momentum conservation
instead of from the particles, which in turn made `E_had` a deterministic
function of the muon kinematics. The calorimetric methods (Jacquet–Blondel,
Sigma) were therefore being fed an input partly built from the very quantity
they exist to replace. It closed by construction, so it could not be wrong — and
for the same reason it could not be a real test.

## The new strategy

Generate in **20 fixed-energy bins**, logarithmic from 10 to 3000 GeV
(`python/flux_bins.py`), with `fixed_lepton_beam 1` so that `x_lepton = 1`: the
lepton carries the full beam momentum and there is **no lepton-side remnant at
all**. The flux is applied afterwards as a per-bin weight,

```
N_i = [ Int_{bin i} f(x) dx ] x sigma_hat(E_i)
```

which reproduces the beam-PDF normalisation by construction. `E_mu` is exact and
equal to the bin energy rather than something sampled inside POWHEG.

Driver: `scripts/generate_binned.sh`. Showering: `python/shower_binned.py`.

## Three bugs this exposed

### 1. The remnant — the one we were hunting

| | flux-as-PDF | binned |
|---|---|---|
| events with negative energy balance | 80 % | **4.5 %** |
| median energy balance | — | **+0.028 GeV** |
| `E_had(particles) / E_had(conservation)` | up to 2.5 | **0.9997** |
| same, restricted to `x > 0.2` | — | **0.9996** |

The particle sum now closes where it used to fail worst, so
`xreco.hadronic_truth()` takes the particle sum again and the momentum-
conservation workaround is retired (it remains available as
`source="conservation"` for comparison).

### 2. Charm below threshold in the old sample

Charm-hadron pair production needs `W ≳ 2 m_D ≈ 3.7 GeV`. Checking the hadronic
invariant mass of every event that contains a charm hadron:

| sample | min `W` | fraction with `W < 3.7 GeV` |
|---|---|---|
| flux-as-PDF, all charm | 1.95 GeV | 2.8 % |
| flux-as-PDF, `E_mu < 20 GeV` | 2.02 GeV | **56.8 %** |
| binned, all charm | 3.28 GeV | **< 0.05 %** |
| binned, `E_mu < 20 GeV` | 4.29 GeV | **0** |

More than half of the low-energy charm in the old sample sat in events with
**insufficient invariant mass to have produced it**. That charm came from the
fictitious beam remnant, not from the DIS vertex.

Because remnant-produced charm is symmetric in `c`/`c̄`, it *diluted* the
fragmentation asymmetry — but only where it lived, i.e. at low energy and low
`x`. The inclusive tagged asymmetry moves from `-0.061 ± 0.036` to
`-0.087 ± 0.024`, and the low-`x` bins move by up to 2.9 sigma, while the signal
region `x > 0.2` is unchanged within 0.9 sigma. See STATUS.md for the bin-by-bin
comparison.

### 3. A bug of our own: non-unit POWHEG weights

The first binned pass used a constant weight per bin. POWHEG events are not unit
weight — in a typical bin ~3 % carry a **negative** weight (NLO subtraction) and
a few carry a large positive one. `shower_dis.py` had always used
`infoPython().weight()` per event; the binned script now does too.

With this fixed, the closure against the old production is:

| | flux-as-PDF | binned | ratio |
|---|---|---|---|
| total DIS yield | 8.99e5 | 8.69e5 | **0.967** |
| `<E_mu>` weighted | 1001 GeV | 940 GeV | 0.94 |

3.4 % on the absolute normalisation, between two generation strategies that
share nothing but the flux grid and the PDF set. That is the closure test.

## Negative weights dominate the charm statistics

The effective statistics `N_eff = (sum w)^2 / sum w^2` are far below the raw
count, and the loss is concentrated in exactly the subsample we care about:

| subsample | `N_raw` | `N_eff` | retained |
|---|---|---|---|
| inclusive (binned) | 3 335 495 | 1 604 615 | 48 % |
| **contains charm** | 296 111 | 10 314 | **3.5 %** |

The mechanism: **37 % of charm events carry a negative weight**, against 4 %
of non-charm events, while `|w|` is essentially the same for both. POWHEG's
negative-weight events populate the hard-radiation region, and that is where
`g -> c cbar` happens, so charm is genuinely enriched in the negative-weight
class.

This is not a bug, but it has a practical consequence that was previously
mistaken for a puzzle: it is why the tagged sample had `N_eff = 778` for 12 019
raw events. **Generating more events buys ~14x less than the raw count suggests
for any charm observable.**

## Allocation: the remaining inefficiency

Generating the same 10 k events in every bin does not match where the rate is.
Per-bin share of the tagged-charm yield:

| bins | share of yield | share of generated events |
|---|---|---|
| 0–11 (11–266 GeV) | ~2.5 % | 60 % |
| 12–15 (353–831 GeV) | ~17 % | 20 % |
| **16–19 (1.1–2.6 TeV)** | **~80 %** | 20 % |

**This has been done.** Bins 12–19 were regenerated with 80 k events each into a
separate directory and merged by `shower_binned.py --indir <dir1> <dir2>`, which
divides the per-bin weight by the *combined* event count so the normalisation is
untouched and only `N_eff` improves. Bins 12–15 were included as well as 16–19
because the `x > 0.2` region is spread more broadly in energy than the inclusive
charm sample: bins 16–19 carry 80 % of all tagged charm but only 45 % of the
`x > 0.2` yield, while bins 12–15 carry a further 33 %.

The merged sample is **3 335 495 events**. At `x > 0.2` the tagged-charm
`N_eff` rose from 55 to 280, and the hadron-composition `N_eff` from 619 to
3 468. The `x > 0.2` asymmetry moved from `-0.402 ± 0.134` at low statistics to
**`-0.264 ± 0.060`** — i.e. the low-statistics value was a fluctuation, and the
converged number agrees with the original production's `-0.198 ± 0.036` at
0.9 sigma.

## Two things that are correct and look wrong

**`ebeam2` must equal the PDG mass exactly.** The original production used
`0.940` for the neutron (PDG `0.93957`). With the flux-as-beam-PDF setup the
beam remnant quietly absorbed the resulting 28 MeV of target momentum; with
`fixed_lepton_beam 1` there is no remnant, and Pythia fails *every* neutron event
with `ProcessContainer::constructProcess: setting mass failed`. The proton was
unaffected only by luck (`0.938 < m_p`, so Pythia clamped it to rest).

**`LesHouches:matchInOut = off` is required.** POWHEG writes fixed-target events
in a frame where the lepton carries `E = ebeam1 + m_target/2`, so the light-cone
energy of the final state does not match the beams reconstructed from the
`<init>` line to Pythia's 1e-6 tolerance. Pythia's incoming-parton reassignment
(`ProcessContainer.cc:1008`) then rejects every event — it killed 14 of the 40
neutron bins. The reassignment is a no-op wherever POWHEG already conserves
momentum exactly, i.e. everywhere: a bin that worked with it on is unchanged to
5 significant digits (151.35902 vs 151.35908 GeV mean final-state energy).

## Frames

The muon is read from `py.process` (hard process) and the hadrons from
`py.event` (after shower and decay). That is deliberate. The target is at rest in
both records to 2 MeV, and `Q^2` and `x` are invariant under the residual
longitudinal boost. The 1.6 % difference between `nu` measured in the two records
is FSR off the muon leg, which is physics: the hard-process muon is the right
truth for a PDF study, the final-state one is what a spectrometer would see. The
combination closes energy conservation to **+0.07 GeV on a 113 GeV event**.
