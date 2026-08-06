# Both charm hadrons, not just the tagged one

*2026-08-06. Extraction: `python/charm_pair.py`; cache
`/eos/home-e/evilla/faser/charm_pair_v2.npz` — the merged binned production,
3 335 495 events (base 10 k/bin plus 80 k/bin extensions of bins 12–19).*

The analysis tags charm with the **leading semileptonic muon**, which silently
drops the second charm hadron. This note asks what that costs, and whether the
second hadron is a usable signature in its own right.

## There are always two

**99.46 %** of charm events contain two weakly-decaying charm hadrons (294 511 of
296 111 raw). An earlier verbal claim — that in the flavour-excitation topology
you see only one charm hadron — was **wrong**. Both charm quarks always end up in
charm hadrons; what differs between topologies is their kinematic asymmetry, not
their number.

## The second one is not soft

| | leading | sub-leading |
|---|---|---|
| median `p` | 242 GeV | **67 GeV** |
| median `theta` | 5.0 mrad | 17.2 mrad |
| median decay length `L = (p/m) c tau` | 20.2 mm | **5.2 mm** |
| `L` > 1 cm (one 3DCal voxel) | 76 % | **32 %** |
| `L` > 3 mm | 94 % | 63 % |
| `p` > 100 GeV | 83 % | 37 % |

A third of sub-leading charm hadrons fly more than one 1 cm cube before decaying.
That is a displaced vertex, not a point-like deposit.

**This is necessary, not sufficient.** There is no reconstruction simulation
here: resolving a second vertex also needs the daughter tracks (2–3 charged, few
GeV) separated inside a hadronic shower depositing tens of GeV at the same place.
Voxel clustering and shower overlap are not modelled. Conversely, 1 cm is the
*sampling pitch*, not the vertex resolution — with ~200 samplings along a track
the pointing is better than the pitch, as it is for the 0.3 mrad muon angle.

## Hadron composition, split by charm side

This is the table the fragmentation argument needs, and it does not exist in any
other cache: the per-hadron cache has species but no `x`, and the per-event cache
has `x` but no species.

Renormalised to D0 + D+ + Lambda_c (D_s excluded, so the literature column is
comparable; D_s is ~8 % in data).

**Inclusive over x** — `N_raw` = 532 322, `N_eff` = 19 509:

| hadron | literature | our MC, c side | our MC, cbar side | c/cbar |
|---|---|---|---|---|
| D0 | 66.4 % | 53.1 % | 62.8 % | 0.84 |
| D+ | 26.8 % | 26.9 % | 34.3 % | 0.78 |
| **Lambda_c** | 6.8 % | **20.0 %** | **2.9 %** | **6.90** |

**At `x_truth > 0.2`** — `N_raw` = 6 301, `N_eff` = 3 468:

| hadron | literature | our MC, c side | our MC, cbar side | c/cbar |
|---|---|---|---|---|
| D0 | 66.4 % | **21.0 %** | 63.0 % | 0.34 |
| D+ | 26.8 % | **10.9 %** | 36.2 % | 0.30 |
| **Lambda_c** | 6.8 % | **68.1 %** | **0.8 %** | **87.4** |

The point for a presentation: **the cbar side reproduces the literature at every
x** (63.0 % D0 at large x against 66.4 % inclusive in data), while the c side
transforms completely. Only one side moves. That is the "hadronisation against
valence versus sea" statement made quantitative — the cbar has only vacuum
qqbar pairs available at any x, whereas the c gains access to the valence diquark
as x grows and the struck quark and remnant come close in rapidity.

The inclusive world average of 6.8 % is an average over both sides and therefore
hides the effect, which is why fragmentation looks "universal" in the literature.

### Evolution with x

| x | f(Lambda_c) c side | f(Lambda_c) cbar side | Lambda_c/Lambda_cbar | N_eff |
|---|---|---|---|---|
| < 0.01 | 11.0 % | 3.5 % | 3.1 | 11 282 |
| 0.01–0.05 | 16.1 % | 2.4 % | 6.5 | 3 129 |
| 0.05–0.10 | 25.8 % | 3.6 % | 7.1 | 1 670 |
| 0.10–0.20 | 56.4 % | 0.7 % | 75.2 | 4 525 |
| 0.20–0.30 | 66.8 % | 1.0 % | 70.0 | 2 199 |
| **0.30–1.00** | **70.3 %** | 0.5 % | **148.3** | 1 269 |

Raw counts behind the extremes: 2 094 Lambda_c against 26 Lambda_cbar at
`x > 0.2`, and 770 against 5 at `x > 0.3`.

## Where the Lambda_c goes — the asymmetry mechanism, explicitly

Species fractions among charm events with a pair:

| species | leading | sub-leading |
|---|---|---|
| D0 | 55.3 % | 51.6 % |
| D+ | 30.2 % | 26.2 % |
| Ds+ | 8.4 % | 6.3 % |
| **Lambda_c** | **5.4 %** | **15.6 %** |

At `x_truth > 0.2` the sub-leading hadron is a Lambda_c **43.8 %** of the time,
with mean charm sign **+0.996** — i.e. essentially every one is a Lambda_c+,
containing the `c` and not the `cbar`.

The sum rule appears by itself at hadron level, as it must (one `c` and one
`cbar` per event):

| | A(leading) | A(sub-leading) | sum |
|---|---|---|---|
| all | −0.030 | +0.009 | **−0.021** |
| `x > 0.2` | −0.006 | +0.006 | **+0.000** |

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
| hardest muon *(what the analysis does)* | −0.0813 | −0.2637 |
| softest muon | −0.0776 | −0.2968 |
| average of both | −0.0794 | −0.2803 |
| **spread from the choice** | **0.004** | **0.033** |
| statistical error | ±0.023 | ±0.060 |

The selection effect is **6x smaller than the statistical error** overall and
**2x smaller** at large `x`. The asymmetry survives deliberately tagging on the
wrong muon.

Reason: only **5.6 %** of tagged events have both charm hadrons semileptonic
(4.4 % at `x > 0.2`, 25 raw events). In 94 % of events there is no choice to make.

## DIS muon vs charm muon: safe where it matters

The analysis identifies the scattered muon by **ancestry**, i.e. with truth
information a detector does not have. A realistic criterion would be "the hardest
muon is the DIS muon". How often would that fail?

| truth `x` | failure rate | `<nu>` |
|---|---|---|
| < 0.01 | **4.4 %** | 839 GeV |
| 0.05–0.10 | 2.4 % | 150 GeV |
| 0.10–0.20 | 0.6 % | 57 GeV |
| **> 0.2** | **1.1 %** | 54 GeV |

The confusion lives at very low `x`, where `nu` is huge and the scattered muon is
soft. At `x > 0.2` the mean scattered-muon momentum is **758 GeV** against
**6.2 GeV** for the hardest charm muon — a factor **122**. The rate itself rests
on 3 raw events out of 503 and is poorly measured, but the kinematic separation does not
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
| 0 semileptonic | 82.3 % | 2.2 % |
| **1 semileptonic** | 16.8 % | **10.9 %** |
| >= 2 semileptonic | 0.9 % | **20.2 %** |

A factor ~12 between the zero-muon and two-muon categories, so **counting the
muons tells you which correction to apply** — this is what the `subtract_nu` flag
does, but currently on average rather than per category. The gain is in the
1-muon category (17 % of charm events), not the 2-muon one (1 %). At `x > 0.2`
the loss is the same 12.0 %, so the bias does not worsen where it matters.

The correction remains statistical: you learn that energy is missing and how much
on average, not how much in that event.

The 2 % in the "0 semileptonic" row is not an error — those are semi-**electronic**
decays, which produce neutrinos without a muon.

**A third use, the most solid.** The opposite-sign dimuon is genuinely
opposite-sign **98.3 %** of the time. At 1 % of charm events it is too rare to
measure with, but it is a clean control sample against pi/K decay-in-flight
background — which remains unmodelled (see the note's caveat list).

## What this does NOT validate

The approximation is safe for the **asymmetry**. It says nothing about
**selection efficiency and purity**: the ML tagger numbers used in the note
(21 % efficiency at 82 % purity for `c -> H`) are taken from the CDR and have not
been validated against this topology or this muon multiplicity.
