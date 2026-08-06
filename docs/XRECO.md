# Bjorken-x reconstruction in FASERcal — RESULT

*2026-08-03. Answers the "showstopper" question. Report: `reports/fasercal_xreco.pdf`.*

## Bottom line

**The showstopper is real for the lepton-only method, and FASERcal's calorimeter
substantially rescues it.** At the CDR's own resolutions, the large-x charm
signal retained (efficiency × purity) is:

| method | x ≥ 0.2 | vs lepton-only | if E_in **not** measured |
|---|---:|---:|---|
| lepton-only (what the paper uses) | 0.149 | — | **unavailable** |
| Jacquet–Blondel | 0.288 | ×1.9 | 0.302 |
| **Σ (fixed-target)** | **0.421** | **×2.8** | **0.421 — unchanged** |
| double-angle | 0.491 | ×3.3 | **unavailable** |

*Only x ≥ 0.2 is quoted. Above x = 0.4 the perturbative sample has **12 raw
events**, so nothing there is statistically meaningful.*

**Σ is the method to use.** Double-angle scores marginally higher at the nominal
point but degrades steeply with the hadronic angular resolution, whereas **Σ is
flat across the whole range** —
it depends on E + p_z, which is first-order insensitive to θ_h. It is also the
only method that does not require the incoming muon energy to be measured.

### The incoming muon energy decides which methods even exist

The downstream spectrometer sees only the **outgoing** muon. If E_in is not
measured anywhere, then:

| method | dependence on E_in | consequence |
|---|---|---|
| lepton-only | ν = E_in − E′, amplified by 1/y | **impossible** |
| double-angle | x ∝ E_in (linear) | **impossible** |
| Jacquet–Blondel | only via the mild (1−y) | usable |
| **Σ** | **reconstructs it from the final state** | **works by design** |

So in that scenario **Σ is the only good method left**, and it loses nothing.

### How to measure E_in in FASERcal: you already can

FASERcal cannot measure the incoming muon momentum directly — there is no field
upstream of the target, the muon does not stop so there is no range measurement,
and at ~665 GeV multiple scattering is far too small to be used. The only two
options are an upstream spectrometer, or **energy conservation using the final
state**. The second is free:

    E_in = E'_mu (spectrometer)  +  E_had (calorimeter)  -  M

**This is numerically identical to the Σ method** (both give median 532 GeV, rel.
error IQR 0.37) — at these small angles E + p_z ≈ 2E, so the two constructions
coincide. And it works for a reason worth stating explicitly:

    nu = E_in(reco) - E'(reco) = (E' + E_had - M) - E' = E_had - M

**The 47 % spectrometer error appears in both terms and cancels exactly.** The
resulting ν resolution is **13 %**, against **410 %** for lepton-only, where E_in
and E′ are measured independently and their errors add instead of cancelling. It
is an algebraic cancellation, not a tuned improvement — which is why Σ is robust.

Note JB gets *better* without a measured E_in (0.288 → 0.302): the
Σ-reconstructed E_in, coming from the 9 % calorimeter, is more accurate than the
54 % spectrometer measurement. Even when E_in *is* measured, it is better to
reconstruct it from the final state.

> The honest hypothesis — *"FASERcal may reconstruct large-x better than FASERν
> despite worse angular resolution, because it has an energy measurement emulsion
> cannot provide"* — is **supported**.

## How much of the selected sample is actually intrinsic charm?

The cut x ≥ 0.2 is well chosen: at **truth** level the charm above it is almost
entirely IC-attributable (IC fraction = 1 − N_perturbative/N_fitted):

| x bin | IC fraction |
|---|---:|
| 0.1–0.2 | 85 % |
| 0.2–0.3 | **98 %** |
| 0.3–0.4 | 99 % |
| > 0.4 | ~100 % |
| **integrated x ≥ 0.2** | **~99 %** |

**But that is not what you select.** Migration from low x — where there is no IC —
dilutes it, and the dilution is just the purity in physics terms:

| selection | IC fraction | rel. IC yield |
|---|---:|---:|
| truth x ≥ 0.2 | 99 % | — |
| lepton-only | 82 % | 0.28 |
| Jacquet–Blondel | 38 % | 1.74 |
| Σ | 64 % | 1.00 |
| double-angle | 72 % | 1.01 |
| **Σ AND double-angle** | **95 %** | 0.58 |
| all three | 97 % | 0.43 |

So a single method leaves the selected sample **one-third ordinary charm** (Σ) to
**nearly two-thirds** (JB). Note lepton-only gives the purest *single*-method
sample (82 %) but keeps only 24 % of the signal, so in absolute terms it has
~3.5× fewer IC events than Σ.

> **Caveat.** The perturbative sample above x = 0.2 has **N_eff = 7** (966 raw
> events, one carrying 9 % of the weight). The 99 % is robust only because the
> perturbative contribution is small to begin with — a factor-2 error still
> leaves ~97 %. But the perturbative *yield itself* should not be quoted, and any
> ratio-based statement there needs a dedicated high-statistics `production_pc`.

## Combining methods — it works, and substantially

The methods fail for **different reasons**, so their errors are not all
correlated. Correlation of log(x_reco/x_true) on charm events:

| | lepton | JB | Σ | DA |
|---|---:|---:|---:|---:|
| lepton-only | 1.00 | −0.03 | 0.65 | −0.05 |
| Jacquet–Blondel | −0.03 | 1.00 | −0.01 | **0.94** |
| Σ | 0.65 | −0.01 | 1.00 | **0.08** |
| double-angle | −0.05 | 0.94 | 0.08 | 1.00 |

**Σ and double-angle are nearly independent (0.08)** — Σ is driven by the
calorimetric energy, double-angle by the angles. JB and double-angle are
redundant (0.94), both being hadronic-system-driven.

Requiring both:

| selection | efficiency | purity | eff × purity |
|---|---:|---:|---:|
| lepton-only | 23.7 % | 63.8 % | 0.151 |
| Σ | 77.1 % | 57.6 % | 0.444 |
| double-angle | 82.5 % | 61.3 % | 0.506 |
| **Σ AND double-angle** | **66.3 %** | **85.7 %** | **0.568** |
| all three | 54.9 % | 94.7 % | 0.520 |

**Σ AND double-angle is the best selection: 0.568, ×3.8 over lepton-only**, and it
lifts the IC fraction from 64 % to **95 %**. Adding JB buys purity (94.7 %) but
costs more efficiency than it gains, because JB is redundant with double-angle.

Two practical notes:

- **The combination needs E_in measured**, since double-angle does. If it is not,
  the fallback is **Σ AND Jacquet–Blondel** (their correlation is −0.01, also
  independent): efficiency 60.8 %, purity 84.5 %, product 0.514 — still ×3.4 over
  lepton-only.
- This is a simple AND of two cuts. A proper **constrained kinematic fit** using
  all measurements with their covariances (as in arXiv:2506.13889 §4.2) should do
  better still, and is the natural next step.

## Why lepton-only fails

| | median |
|---|---:|
| E_in | 665 GeV |
| E′_µ | 469 GeV |
| **ν = E_in − E′_µ** | **52 GeV** |
| **y** | **0.161** |

ν is only **7.7 %** of E_in — a small difference of two large numbers — so a
fractional error ε on each muon leg propagates to ≈ ε/y on ν. The CDR
spectrometer gives **σ_p/p = 47 %** at the median scattered-muon momentum
(iron-core, multiple-scattering dominated: 20 % floor at 20 GeV, 63 % at 1 TeV),
so **ν carries a ~290 % error**. The migration matrix shows essentially **no
correlation** between x_true and x_reco.

This is not a defect of the spectrometer — it is built for muon ID and charge,
where it excels (0.7 % misidentification). It is the wrong instrument for x.

## Why Σ works

    2 E_in = Σ_had(E + p_z) + E′(1 + cos θ_µ) − M

reconstructs the incoming energy **from the final state**, then proceeds as
lepton-only with that E_in. Three properties matter here:

1. **No incoming-muon measurement needed** — important, since the downstream
   spectrometer only sees the outgoing muon.
2. **First-order insensitive to θ_h.** p_z = |p| cos θ_h and θ_h ≈ 81 mrad, so a
   20 mrad smearing changes cos θ_h by 0.2 %. Hence the flat curve.
3. **Driven by the calorimeter (≈9 %), not the spectrometer (47 %).**

(Note the collider Σ method uses E − p_z; for a fixed target with the beam along
+z that combination is identically M and carries no information, so the **E + p_z**
form is the correct analogue. This is derived, not taken from the literature.)

The **angular** worry turned out to be a non-issue for the muon: the 3DCAL fits a
track over ~200 layers, giving σ_θ ≈ σ_hit·√(12/N)/L ≈ **0.3 mrad** against a
4.3 mrad median scattering angle. The "cm voxels vs µm emulsion" concern is
largely defused by the number of samples.

## Method definitions

Fixed-target kinematics, masses neglected, x = Q²/(2Mν):

| method | formula | needs |
|---|---|---|
| lepton-only | Q² = 2E_in E′(1−cos θ_µ), ν = E_in − E′ | both muon energies |
| Jacquet–Blondel | ν = E_had − M, Q² = p_T,had²/(1−y) | calorimeter only |
| Σ | E_in from E + p_z balance, then as lepton-only | calorimeter + muon angle |
| double-angle | E′ = E_in sin θ_h / sin(θ_µ+θ_h) | angles + E_in |

**JB needs the target-mass subtraction**: the measured hadronic energy is ν + M,
and omitting M biases x_JB low by M/ν — 9 % on the x ≥ 0.2 sample where ν ≈ 9 GeV.

## Inputs and their provenance

| quantity | value | source |
|---|---|---|
| muon momentum resolution | 20 % (20 GeV) → 63 % (1 TeV), applied as a Gaussian in **1/p** | **CDR §2.6.2** |
| hadronic energy resolution | 9 % | **Bern talk** (p_jet, NC) |
| muon angular resolution | 0.3 mrad | derived from 3DCAL geometry |
| **hadronic angular resolution** | **20 mrad** | **ASSUMPTION — not from any document.** Scanned 0–100 mrad; the conclusion is stable because Σ is insensitive to it |

## Validation

All four methods return **x_reco/x_true = 1.0000** with a perfect detector, so
the framework is verified before any smearing is applied.

## A methodological problem, and how it was solved

> **Resolved at source, 2026-08-06.** The events have been regenerated in fixed
> energy bins, which removes the beam remnant entirely. The hadronic four-vector
> is now the **particle sum** again and the workaround below is retired. It
> survives as `xreco.hadronic_truth(d, source="conservation")` for comparison.
> See [BINNED_PRODUCTION.md](BINNED_PRODUCTION.md). The rest of this section
> records the problem and why the workaround was necessary at the time.

The hadronic four-vector was computed **from momentum conservation**
(X = q + p_target ⇒ E_had = ν + M, p_T = p_out sin θ_µ, p_z = |p_in| − p_out cos θ_µ),
**not** by summing Pythia final-state particles.

The particle sum was contaminated. The muon flux was implemented as a beam PDF of
a fictitious 7 TeV beam, and that beam's remnant deposits a few GeV of hadronic
activity that is not part of the DIS vertex. Energy conservation
(E_in + M − E′ − E_had ≥ 0) failed in **80 % of events**, rising to **98 %** on the
x ≥ 0.2 charm sample, where E_had/ν reached **2.5**. Negligible at large ν
(E_had/ν = 1.02 above 500 GeV), it dominated exactly where the physics is.

Three fixes were tried and **all failed**:

| attempt | result |
|---|---|
| drop beam-remnant particles (status 61–63) | no change — the remnant *partons* are not final; their hadronic descendants are |
| drop very forward particles (θ < 1 mrad) | removes 28 GeV on average but does not move the median |
| require ancestry to reach beam 2 | *worse* (E_had/ν = 8.5) — Pythia's colour reconnection leaves the mother graph too connected to separate the two beam sides |

Conservation was exact and made the closure exact by construction, which is what
a resolution study needs.

### What the workaround was costing

It made E_had a deterministic function of the muon kinematics, so Jacquet–Blondel
and Σ were partly fed the very quantity they exist to replace. On the binned
sample the two definitions can finally be compared on the same events:

| method | particles (bias / half-width) | conservation |
|---|---|---|
| lepton-only | −0.283 / 0.702 | −0.283 / 0.702 |
| Jacquet–Blondel | −0.049 / 13.71 | −0.001 / 14.31 |
| **Σ** | +0.029 / **1.213** | +0.010 / **1.130** |
| double-angle (Σ E_in) | −0.010 / 3.225 | +0.006 / 3.267 |
| double-angle | −0.005 / 3.563 | −0.001 / 3.532 |

Lepton-only is identical in the two columns, as it must be — it never touches the
hadrons — which validates the comparison. The workaround was flattering the
calorimetric methods: it hid a 5 % bias on JB and a 3 % bias on Σ, and improved
Σ's resolution by 7 %. Those are the real fragmentation fluctuations and
semileptonic-neutrino losses, previously suppressed by construction.

**The physics conclusion is unchanged**: Σ and the double-angle method remain far
ahead of lepton-only.

## Are the calorimetric methods viable given neutral particles? — Yes.

This was the main threat to the whole argument, and it is now **modelled, not
assumed away**.

**Neutrinos from charm semileptonic decays** escape entirely. This is an
irreducible, one-sided loss and it is *specific to the signal events*. Measured
on the 25-seed sample:

| | value |
|---|---|
| charm events with a charm neutrino | **35 %** |
| when nonzero, median E_ν | **13.2 GeV = 8.4 % of ν** |
| x ≥ 0.2 charm events losing >10 % of ν | **21 %** |

Effect on the result (eff × purity at x ≥ 0.2):

| method | no ν loss | **with ν loss** | change |
|---|---:|---:|---:|
| lepton-only | 0.151 | 0.151 | 0 % (uses no E_had) |
| Jacquet–Blondel | 0.292 | 0.262 | **−10 %** |
| **Σ** | 0.433 | **0.444** | +3 % |
| double-angle | 0.506 | 0.506 | 0 % |

**It costs Jacquet–Blondel 10 % and Σ nothing.** Σ is protected because the loss
biases ν low, which pushes x_Σ *up* and hence into the selection — efficiency and
purity move in opposite directions and largely cancel.

*Caveat:* the 9 % hadronic resolution from the Bern talk is derived from full
simulation of neutrino events, which also contain charm decays, so this may
partially double-count. The effect is switchable (`subtract_nu`) for that reason.

**Neutral hadrons** (n, K_L), which respond poorly in a non-compensating
calorimeter, are only **8 %** of the visible energy on the x ≥ 0.2 charm sample
(EM 22 %, charged hadrons 63 %). Too small to overturn anything, and already
folded into the quoted σ_had.

**Still not modelled:**

- Non-compensation as an explicit e/h response difference with event-by-event π⁰
  fluctuations — subsumed into the single σ_had parameter.
- The ~9 % is quoted for the full 3DCAL+ECAL+AHCAL chain; the 3DCAL alone will be
  worse, and this study excludes the AHCAL from the fiducial volume.

## Plots

| page | content |
|---|---|
| 1 | truth-x → reco-x **migration matrices**, all four methods. Lepton-only shows almost no correlation; Σ retains a clear diagonal |
| 2 | **efficiency × purity vs σ_θh**, at x ≥ 0.2 and 0.4. Shows Σ's flatness and double-angle's degradation |

## Next

1. Fold in the neutrino energy loss for charm events — the one unmodelled effect that works against the calorimetric methods.
2. Get a real σ_θh from the FASERcal reconstruction, replacing the 20 mrad assumption.
3. Redo with the 3DCAL-only hadronic resolution if the AHCAL is excluded.
4. Propagate to a fitted-vs-perturbative discrimination significance. Note the current MC has only **12 perturbative charm events above x = 0.4**, so ratio-based statements there are statistics-limited and need more seeds.
