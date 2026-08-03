# Bjorken-x reconstruction — status and a blocker found

*2026-08-03. Study in progress; this records what the first pass established and
the problem that stopped it.*

## The question

Intrinsic charm shows up as a **large-x excess of charm**: fitted/perturbative
= 1.8× at x = 0.2, 5.7× at 0.4, 10.9× at 0.5, 44× at 0.7. arXiv:2506.13889 §4
reconstructs x from the muon alone and finds that **σ_p = 10 % drops the x ≥ 0.2
excess from ×8 to ~×2, and σ_p = 30 % removes it entirely.**

FASERcal's spectrometer sits at **σ_p/p ≈ 47 %** for the median scattered muon
(CDR §2.6.2) — far inside the "signal gone" regime. The hypothesis to test is
that **calorimetric** methods survive where the lepton-only method dies, because
the hadronic resolution (~9 %) is five times better than the muon momentum.

## Why the lepton-only method is structurally fragile — quantified

| quantity | median |
|---|---:|
| E_in | 665 GeV |
| E′_µ | 469 GeV |
| **ν = E_in − E′_µ** | **52 GeV** |
| **y = ν/E_in** | **0.161** |
| Q² | 4.3 GeV² |

**ν is only 7.7 % of E_in.** It is a small difference of two large numbers, so a
fractional error ε on each muon leg propagates to roughly ε/y on ν:

> σ_p/p = 47 % on each leg → **~290 % error on ν** → x is destroyed.

This is the mechanism behind the paper's result, now with FASERcal's real
numbers. It is not a subtle effect.

The **angular** term, by contrast, is *not* the problem — contrary to the initial
worry that cm voxels would be fatal versus µm emulsion. The 3DCAL measures a
track with 1 cm voxels over ~200 layers, so the straight-line fit gives
σ_θ ≈ σ_hit·√(12/N)/L ≈ **0.3 mrad**, against a median scattering angle of
4.3 mrad — a 7 % effect. Multiple scattering adds ~0.06 mrad at 470 GeV. **The
granularity penalty is largely defused by the number of samples.**

## The blocker: the cached hadronic energy is contaminated

With a **perfect detector** (all resolutions set to zero), the four methods
should each return x_reco/x_true = 1. They do not:

| method | perfect-detector closure |
|---|---:|
| lepton-only | **1.000** ✓ |
| Jacquet–Blondel | 0.620 ✗ |
| Σ | 0.520 ✗ |
| double-angle | 0.405 ✗ |

Lepton-only closing exactly validates the framework. The others are biased by a
**forward energy excess in the cached E_had**:

| ν range [GeV] | ⟨ν⟩ | ⟨E_had⟩ | E_had/ν |
|---|---:|---:|---:|
| 0–5 | 4.0 | 5.8 | **1.47** |
| 5–20 | 11.0 | 13.8 | 1.17 |
| 50–150 | 85.3 | 96.6 | 1.05 |
| 500–5000 | 790 | 845 | 1.02 |

E_had − ν has a roughly constant **few-GeV forward offset** (median +3.4 GeV,
and pz ≫ pT in the low-ν bin, so it is beam-directed). Physically E_had should
be ν + M. The offset is negligible at large ν but **dominates at large x, where
ν is small** — exactly the region the whole study is about. On the x ≥ 0.2 charm
sample it reaches **E_had/ν = 2.5**, which is what biases the JB ν by ~2.5 and
hence x_JB by ~0.4.

**Cause:** the hadronic sum keeps beam-remnant particles. The muon flux is
implemented as a beam PDF of a fictitious 7 TeV beam, and removing only the
status-62 remnant photon (as the FASERν study does, which is sufficient for
*its* observables) leaves other remnant activity in the sum.

**This is a defect in my extraction, not in the methods.** No conclusion about
JB / Σ / double-angle can be drawn until it is fixed.

> **Correction to an earlier statement.** I said this study needed no
> re-showering because `e_had` and the hadronic vector were cached. That was
> wrong: the cache stores only the *summed* hadronic four-vector, so the remnant
> cannot be removed retroactively. Re-showering is required.

## The fix, now running

`shower_dis.py` now records **three** hadronic definitions per event:

| key | definition |
|---|---|
| `e_had`, `p*_had` | all non-neutrino non-muon final state (the old, contaminated one) |
| `e_had_nr`, `p*_nr` | additionally drops **all** beam remnants (statusAbs 61–63) |
| `e_had_co`, `p*_co` | additionally drops **very forward** particles (θ < 1 mrad), which in a real detector go down the beam line and are not associated with the vertex |

Re-showering both PDF hypotheses (10 seeds each). Once done, the closure test is
repeated: whichever definition returns 1.000 for JB with a perfect detector is
the physically correct hadronic system, and the smearing study proceeds from
there.

## Next

1. Re-run the perfect-detector closure on all three definitions; pick the one that closes.
2. Truth-x → reco-x migration matrices, four methods, at the real FASERcal resolutions (σ_p from the CDR curve, σ_E_had ≈ 9 %, σ_θ ≈ 0.3 mrad).
3. Surviving fitted/perturbative ratio at x_reco ≥ 0.2 and ≥ 0.4 — the direct answer to the showstopper question.
4. Scan σ_E_had, since it is the parameter the calorimetric methods live or die by, and the ~9 % is quoted for the full 3DCAL+ECAL+AHCAL chain rather than the 3DCAL alone.
