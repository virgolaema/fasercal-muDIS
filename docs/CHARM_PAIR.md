# Both charm hadrons, not just the tagged one

*2026-08-06. Extraction: `python/charm_pair.py`; cache
`/eos/home-e/evilla/faser/charm_pair_v1.npz` (base binned production).*

The analysis tags charm with the **leading semileptonic muon**, which silently
drops the second charm hadron. This note asks what that costs, and whether the
second hadron is a usable signature in its own right.

## There are always two

**99.56 %** of charm events contain two weakly-decaying charm hadrons (47 331 of
47 541 raw). An earlier verbal claim — that in the flavour-excitation topology
you see only one charm hadron — was **wrong**. Both charm quarks always end up in
charm hadrons; what differs between topologies is their kinematic asymmetry, not
their number.

## The second one is not soft

| | leading | sub-leading |
|---|---|---|
| median `p` | 135 GeV | **37 GeV** |
| median `theta` | 7.8 mrad | 23.5 mrad |
| median decay length `L = (p/m) c tau` | 11.2 mm | **3.0 mm** |
| `L` > 1 cm (one 3DCal voxel) | 75 % | **31 %** |
| `L` > 3 mm | 93 % | 62 % |
| `p` > 100 GeV | 83 % | 37 % |

A third of sub-leading charm hadrons fly more than one 1 cm cube before decaying.
That is a displaced vertex, not a point-like deposit.

**This is necessary, not sufficient.** There is no reconstruction simulation
here: resolving a second vertex also needs the daughter tracks (2–3 charged, few
GeV) separated inside a hadronic shower depositing tens of GeV at the same place.
Voxel clustering and shower overlap are not modelled. Conversely, 1 cm is the
*sampling pitch*, not the vertex resolution — with ~200 samplings along a track
the pointing is better than the pitch, as it is for the 0.3 mrad muon angle.

## Where the Lambda_c goes — the asymmetry mechanism, explicitly

Species fractions among charm events with a pair:

| species | leading | sub-leading |
|---|---|---|
| D0 | 54.6 % | 52.1 % |
| D+ | 30.8 % | 24.7 % |
| Ds+ | 8.3 % | 6.9 % |
| **Lambda_c** | **6.0 %** | **16.1 %** |

At `x_truth > 0.2` the sub-leading hadron is a Lambda_c **41.5 %** of the time,
with mean charm sign **+1.000** — i.e. essentially every one is a Lambda_c+,
containing the `c` and not the `cbar`.

The sum rule appears by itself at hadron level, as it must (one `c` and one
`cbar` per event):

| | A(leading) | A(sub-leading) | sum |
|---|---|---|---|
| all | −0.059 | +0.048 | **−0.011** |
| `x > 0.2` | +0.045 | −0.045 | **+0.000** |

This is an independent check that the simulation invents no parton-level
asymmetry, consistent with `int (c - cbar) dx = -1.5e-4` on the NNPDF grid.

**So the measured asymmetry comes from branching ratios, not kinematics.** At
large `x` the `c` picks up the remnant diquark to form a Lambda_c+, which is the
softer of the two; Lambda_c has the lowest semileptonic BR of any charm hadron
(4.8 %, against 6.6 % for D0 and 16.7 % for D+). The `c` side therefore loses its
muon tag more often than the `cbar` side, and the tagged sample is `cbar`-rich.

## Does dropping the second muon bias the asymmetry? No.

Recomputing `A` with three different choices of which muon to tag on:

| choice | A (all) | A (`x > 0.2`) |
|---|---|---|
| hardest muon *(what the analysis does)* | −0.1248 | −0.4013 |
| softest muon | −0.1072 | −0.4154 |
| average of both | −0.1160 | −0.4084 |
| **spread from the choice** | **0.018** | **0.014** |
| statistical error | ±0.062 | ±0.134 |

The selection effect is **3.5x smaller than the statistical error** overall and
**10x smaller** at large `x`. The asymmetry survives deliberately tagging on the
wrong muon.

Reason: only **5.7 %** of tagged events have both charm hadrons semileptonic
(0.7 % at `x > 0.2`). In 94 % of events there is no choice to make.

## DIS muon vs charm muon: safe where it matters

The analysis identifies the scattered muon by **ancestry**, i.e. with truth
information a detector does not have. A realistic criterion would be "the hardest
muon is the DIS muon". How often would that fail?

| truth `x` | failure rate | `<nu>` |
|---|---|---|
| < 0.01 | **4.4 %** | 880 GeV |
| 0.05–0.10 | 2.3 % | 178 GeV |
| 0.10–0.20 | 0.5 % | 68 GeV |
| **> 0.2** | **0.8 %** | 51 GeV |

The confusion lives at very low `x`, where `nu` is huge and the scattered muon is
soft. At `x > 0.2` the mean scattered-muon momentum is **881 GeV** against
**6.4 GeV** for the hardest charm muon — a factor **140**. The rate itself rests
on 3 raw events and is poorly measured, but the kinematic separation does not
depend on those events. **The ancestry assumption is harmless in the signal
region.**

## Do we need to see both muons?

**Leptonic analysis: marginal.** It removes the ~5 % DIS/charm mis-assignment,
which is already negligible at large `x`. Lepton-only `x` is dominated by the
muon momentum resolution (half-width 0.70) regardless.

**Calorimetric analysis: more useful.** Semileptonic neutrinos remove energy
invisibly and one-sidedly:

| | share of charm events | `<E_nu / nu>` lost |
|---|---|---|
| 0 semileptonic | 81.9 % | 2.0 % |
| **1 semileptonic** | 17.1 % | **12.0 %** |
| >= 2 semileptonic | 1.0 % | **23.2 %** |

A factor ~12 between the zero-muon and two-muon categories, so **counting the
muons tells you which correction to apply** — this is what the `subtract_nu` flag
does, but currently on average rather than per category. The gain is in the
1-muon category (17 % of charm events), not the 2-muon one (1 %). At `x > 0.2`
the loss is the same 12.3 %, so the bias does not worsen where it matters.

The correction remains statistical: you learn that energy is missing and how much
on average, not how much in that event.

The 2 % in the "0 semileptonic" row is not an error — those are semi-**electronic**
decays, which produce neutrinos without a muon.

**A third use, the most solid.** The opposite-sign dimuon is genuinely
opposite-sign **96.6 %** of the time. At 1 % of charm events it is too rare to
measure with, but it is a clean control sample against pi/K decay-in-flight
background — which remains unmodelled (see the note's caveat list).

## What this does NOT validate

The approximation is safe for the **asymmetry**. It says nothing about
**selection efficiency and purity**: the ML tagger numbers used in the note
(21 % efficiency at 82 % purity for `c -> H`) are taken from the CDR and have not
been validated against this topology or this muon multiplicity.
