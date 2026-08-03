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

The hadronic four-vector is computed **from momentum conservation**
(X = q + p_target ⇒ E_had = ν + M, p_T = p_out sin θ_µ, p_z = |p_in| − p_out cos θ_µ),
**not** by summing Pythia final-state particles.

The particle sum is contaminated. The muon flux is implemented as a beam PDF of a
fictitious 7 TeV beam, and that beam's remnant deposits a few GeV of hadronic
activity that is not part of the DIS vertex. Energy conservation
(E_in + M − E′ − E_had ≥ 0) fails in **80 % of events**, rising to **98 %** on the
x ≥ 0.2 charm sample, where E_had/ν reaches **2.5**. Negligible at large ν
(E_had/ν = 1.02 above 500 GeV), it dominates exactly where the physics is.

Three fixes were tried and **all failed**:

| attempt | result |
|---|---|
| drop beam-remnant particles (status 61–63) | no change — the remnant *partons* are not final; their hadronic descendants are |
| drop very forward particles (θ < 1 mrad) | removes 28 GeV on average but does not move the median |
| require ancestry to reach beam 2 | *worse* (E_had/ν = 8.5) — Pythia's colour reconnection leaves the mother graph too connected to separate the two beam sides |

Conservation is exact and makes the closure exact by construction, which is what
a resolution study needs.

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
