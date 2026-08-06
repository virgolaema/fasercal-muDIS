# FASERcal muon-DIS charm study — status

*Updated 2026-08-03, after the FASERCal Bern CM talk (15 Jul 2026) and CDR v0 (16 Feb 2026).*
*Resolution rule: **the talk wins on conflicts, the CDR fills gaps**.*
*Repo: [github.com/virgolaema/fasercal-muDIS](https://github.com/virgolaema/fasercal-muDIS) (pushed).*

**Where we are:** the full chain from muon-DIS interactions to charge-signed charm
muons is simulated end to end. The Bern talk replaced most of the remaining
assumptions with designed values — and corrected several of them substantially.

> ### 2026-08-06 — the events have been regenerated
>
> The samples are now produced in **20 fixed-energy bins** and reweighted to the
> flux, instead of handing the flux to POWHEG as a beam PDF. Bins 12–19, which
> carry most of the rate, have 80 k events each; the merged sample is
> **3 335 495 events**. See **[BINNED_PRODUCTION.md](BINNED_PRODUCTION.md)** and
> **[CHARM_PAIR.md](CHARM_PAIR.md)**.
>
> **The large-x answer given to A. Rubbia is unchanged.** The regeneration moves
> the *low*-x bins by up to 2.9σ and leaves the signal region alone — which is
> exactly where the defect was, since 56.8 % of the charm below 20 GeV in the old
> sample had `W < 2 m_D`, too little invariant mass to have produced it.
>
> | truth x | flux-as-beam-PDF | binned | tension |
> |---|---|---|---|
> | < 0.05 | +0.011 ± 0.029 | −0.075 ± 0.027 | 2.9σ |
> | 0.05–0.10 | −0.371 ± 0.088 | −0.054 ± 0.079 | 2.7σ |
> | **> 0.2** | **−0.198 ± 0.036** | **−0.264 ± 0.060** | **0.9σ** |
>
> Total DIS yield closes to **3.4 %** between the two strategies. Energy balance:
> 80 % of events negative → **4.5 %**; `E_had(particles)/E_had(conservation)` =
> **0.9997**, and 0.9996 at x > 0.2 where the old ratio reached 2.5.
>
> **N_eff, not N_raw, is the limit.** 37 % of charm events carry a negative POWHEG
> weight against 4 % of non-charm ones, so `N_eff/N_raw = 0.035` for charm against
> 0.48 inclusively — extra events buy ~14× less than the raw count suggests.

---

## 1. Headline numbers

Run 4 at **780 fb⁻¹** (FASERCal's own nominal), **3DCAL only** (AHCAL excluded),
as-designed geometry, **1 mm W baseline**, fiducial z < 1150 mm,
F_flux = 1.10 (measured at the confirmed position), **F_optics = 2.0** (Run 4 beam optics, see below),
identification+linking efficiency ε = 0.9:

| stage | 3DCAL only | + AHCAL | fraction |
|---|---:|---:|---|
| Muon DIS in the target volume | **867 600** | **5 785 000** | 100 % |
| … containing charm | 20 220 | 135 800 | 2.33 % |
| … charm → semileptonic µ | 3 589 | 24 160 | 17.8 % of charm |
| … µ punches out | 3 102 | 20 886 | 86.4 % |
| … reaches spectrometer (43.1 %) | 1 337 | 9 016 | 43.1 % |
| … identified + linked (ε = 0.9) | **1 203** | **8 114** | **0.14 % of DIS** |
| **σ(A_c)** | **0.029** | **0.011** | |

*Binned production, merged sample (3 335 495 events). The charm fraction drops
from 2.86 % to 2.33 % because the old sample contained charm the DIS vertex had
too little invariant mass to produce.*

*Full masses — the CDR's "z < 1150 mm" is **not** a containment cut but the
requirement to be inside the 3DCal, in an older centred coordinate system where
it retains 97.7 %. Removing the erroneous 48 % roughly doubled all yields.*

| target volume | mass | tagged | σ(A_c) |
|---|---:|---:|---:|
| 3DCal only (1 mm W) | 581 kg | 1 203 | 0.029 |
| + ECAL | 1 351 kg | 2 807 | 0.019 |
| + AHCAL | 3 908 kg | 8 114 | 0.011 |
| **+ ECAL + AHCAL** | **4 678 kg** | **9 718** | **0.010** |

Charge-tag purity **99.2 %** (CDR §2.6.2 measured misidentification).

**Two new multiplicative conditions**, both documented in the report's conditions page:

- **F_optics = 2.0** — the muon flux does *not* simply scale with luminosity.
  CDR §1.3.2: the 2024 optics change (crossing angle reversed, Q4 off, looser
  collimators) *"increas[es] the rate by a factor of 2 and increas[es] the energy
  of the muons… a large increase in high-momentum, positively charged muons"*,
  and Run 4 goes to a 250 µrad horizontal crossing angle with mitigation *"not
  clear"*. The FLUKA sample behind both the flux grid and the off-axis map is
  the **Run 3** production, so this is not already in the weights. Two reasons it
  is likely **conservative**: it is quoted for the 2024 configuration, not Run 4;
  and the spectrum also *hardens*, while DIS and charm both rise with E_µ.
  Set `F_OPTICS = 1.0` to recover pure luminosity scaling.

- **AHCAL as an optional second target** — 3.33 t total (~98 % iron), 1.60 t
  fiducial, i.e. **5.8× the 3DCAL fiducial mass**, giving 6.7× the DIS rate.
  Being a sampling calorimeter, interactions in the steel *are* measured, so the
  whole mass counts; the iron composition (w_p = 0.464) is applied. **Caveat:**
  its 4×4 cm² granularity makes vertex finding and muon linking much harder than
  1 cm³ voxels, so the AHCAL curve should be read at a **lower ε** than the
  3DCAL one — the report plots both against ε precisely so they are not compared
  at the same value.

| absorber (3DCal only) | mass | tagged | σ(A_c) |
|---|---:|---:|---:|
| **1 mm W/module** (baseline) | 581 kg | 1 203 | **0.029** |
| 5 mm W/module | 896 kg | 1 855 | 0.023 |
| 10 mm W/module | 1 118 kg | 2 333 | 0.021 |

With the AHCAL included the absorber choice becomes almost irrelevant — the
AHCAL mass dominates.

---

## 2. What is DONE

**Simulation chain**
- [x] Per-event extraction from POWHEG+Pythia8 (`shower_dis.py`): DIS kinematics, hadronic final state, charm hadrons, semileptonic dimuon
- [x] Muon classification by **mother-chain ancestry** (charm-decay vs scattered)
- [x] Charm-sign ≡ decay-muon-charge identity **verified in code**, holds 100 %
- [x] **3 335 495 showered events**, binned production (20 log-E bins × 4 samples; bins 12–19 at 80 k each)
- [x] Both PDF hypotheses: fitted charm (`production_v1`) and perturbative (`production_pc`)

**Normalisation**
- [x] Target composition: scintillator CH (w_p = 0.538) ≠ tungsten (0.402); weights stored target-agnostic
- [x] Reference mass derived from the flux definition (0.724 t)
- [x] Run 4 luminosity 780 fb⁻¹; cross-checked against the TP's Run 4 muon rate

**Detector — now as designed, not assumed**
- [x] 3DCAL geometry from the Bern talk: 10 modules × 20 layers of 1 cm³ cubes, 48×48 cm face
- [x] Absorber options **per module** (1 mm / 5 mm), with material budget, dE/dx and Highland scattering
- [x] Fiducial mass 0.514 t (1 mm) / 0.692 t (5 mm)
- [x] Spectrometer acceptance corroborated against slide 33

**Off-axis flux — measured from FLUKA**
- [x] Built the transverse fluence map (`fluka_offaxis_map.py`) from the EOS FLUKA files
- [x] Evaluated at the **designed** LoS shift (452, 236) mm

**Deliverables** — 10-page PDF report, `docs/TECHNOTE_chain.md`, `docs/REPORT.md`

---

## 2b. DECISIVE: there is no c/c̄ asymmetry in the PDF

See [`docs/CHARM_ASYMMETRY.md`](CHARM_ASYMMETRY.md). Both NNPDF4.0 sets fit the
total c + c̄ and set c − c̄ = 0: the max |c−c̄|/(c+c̄) over the whole grid is
6.4×10⁻³ and the number sum rule gives −1.5×10⁻⁴. **`production_v1` therefore
cannot produce a parton-level A_c — it is zero by construction**, consistent
with the light-cone argument that the minimal |uudcc̄⟩ Fock state gives
A_c(x) = 0 at every x.

What the MC *does* produce is a **fragmentation** asymmetry, via the
leading-baryon effect: Λ_c outnumbers Λ̄_c by **6.8 : 1** (the c picks up two
valence quarks from the target remnant; the c̄ cannot make an anti-baryon), and
Λ_c has the lowest semileptonic BR of any charm hadron (4.8 %). So the c is
preferentially hidden from the muon tag → the tagged sample is c̄-enriched,
A = −0.087 ± 0.024 inclusively, −0.264 ± 0.060 at x > 0.2.

Split by charm side (see [CHARM_PAIR.md](CHARM_PAIR.md)) the mechanism is
explicit: at x > 0.2 the c side is **68.1 % Λ_c** while the c̄ side is **0.8 %**,
and the c̄ side reproduces the world-average fragmentation fractions at every x.
Only one side moves — the c can hadronise against the valence diquark, the c̄
only against the sea.

**Consequences:** (i) A_c is not measurable as an IC observable with current PDF
inputs — it needs an asymmetric-charm set or a model; (ii) the fragmentation
asymmetry is *the same size as our sensitivity*, so it is a competing background,
not a small correction; (iii) **the real IC signal is the large-x enhancement**
(fitted/perturbative = 5.7× at x = 0.4, 10.9× at 0.5, 44× at 0.7), which makes
the x-reconstruction study the decisive question.

---

## 2c. RESULT: the calorimeter rescues large-x — see [`docs/XRECO.md`](XRECO.md)

The showstopper question is **answered**. Large-x signal retained
(efficiency × purity), at the CDR's own resolutions:

| method | x ≥ 0.2 | x ≥ 0.4 |
|---|---:|---:|
| lepton-only (the paper's) | 0.145 | 0.031 |
| Jacquet–Blondel | 0.268 | 0.051 |
| **Σ (fixed-target)** | **0.501 (×3.5)** | **0.145 (×4.7)** |
| double-angle | 0.572 (×3.9) | 0.223 (×7.2) |

**Σ is the method to use** — double-angle scores marginally higher at the nominal
point but degrades steeply with the hadronic angular resolution, while Σ is flat
across 0–100 mrad (it depends on E + p_z, first-order insensitive to θ_h) and is
the only method not needing the incoming muon energy measured.

All four methods close at **1.0000** with a perfect detector, so the framework is
validated before smearing. Report: `reports/fasercal_xreco.pdf` (2 pages:
migration matrices; efficiency × purity vs σ_θh).

Main caveat: **neutrino energy loss from charm semileptonic decays is not
modelled**, and it works against the calorimetric methods. And σ_θh = 20 mrad is
an assumption — though the conclusion is stable because Σ is insensitive to it.

---

## 2d. Why lepton-only fails

The CDR (§2.6.2) gives the measured momentum resolution of the FASERCal muon
spectrometer, from a full Geant4 + GenFit fit. Applied to the muon momenta in
our own sample:

| muon | median p | **σ_p/p** |
|---|---:|---:|
| incoming | 665 GeV | **54.5 %** |
| scattered (DIS) | 469 GeV | **47.0 %** |
| charm decay | 11 GeV | 20.0 % |

(yield-weighted average over the scattered-muon spectrum: **46.8 %**)

Compare with arXiv:2506.13889 §4, the showstopper question:
**σ_p = 10 % already drops the x ≥ 0.2 excess from ×8 to ~×2, and σ_p = 30 %
removes it entirely.** FASERCal's spectrometer is at ~47 %, i.e. well inside the
regime where the paper says the large-x signal is gone.

This is *inherent*, not a flaw: the spectrometer is iron-core (10 planes × 15 cm
Fe ≈ 85 X₀), so it is multiple-scattering dominated — 20 % even at 20 GeV. It is
built for muon identification and **charge**, which it does excellently (0.7 %
misidentification), not for precision momentum on TeV muons.

**But this is exactly the argument for the calorimetric methods.** The Bern talk
quotes the *hadronic* resolution as **~9 %** (p_jet, NC), five times better than
the muon momentum. Jacquet–Blondel reconstructs x from the hadronic system alone
(ν_JB = E_had, Q²_JB from hadronic p_T) and **needs no muon momentum at all**.
So the honest hypothesis is:

> FASERcal may reconstruct large-x *better* than FASERν despite worse angular
> resolution, because it has an energy measurement emulsion cannot provide —
> turning the showstopper into the differentiator.

**Status: not yet tested.** This is the top open item (§5), and it needs no
re-showering: `e_in`, `E_had`, the hadronic vector and Σ(E−p_z) are all cached.

---

## 2e. A_c(x) at large x is NOT measurable — kinematically, not statistically

Meson-cloud models (unlike minimal BHPS) *do* predict c(x) ≠ c̄(x), peaking near
x ≈ 0.5–0.7 with an observable asymmetry of −0.36 to −0.60 and a **sign change at
x = 0.318** — which incidentally *explains* the inclusive null result, since the
sum rule forces the integral to cancel.

We evaluated whether FASERcal could see it. **It cannot**, and the reason is
kinematic:

| true x | charm evts | semilept. µ | punch through | **accepted** |
|---|---:|---:|---:|---:|
| 0.00–0.05 | 140 057 | 25 483 | 21 929 | 6 169 (4.4 %) |
| 0.10–0.20 | 6 942 | 1 140 | 558 | 5 (0.07 %) |
| 0.20–0.32 | 4 583 | 726 | 285 | 1 (0.02 %) |
| **0.32–0.50** | 2 684 | 401 | 165 | **0** |
| **0.50–1.00** | 299 | 51 | 28 | **0** |

**Of 6 213 tagged MC events, not one has x ≥ 0.32.** Large x means small ν, so the
charm is soft and the decay muon is soft and wide-angle — median **1.9 GeV** and
**257 mrad** at x ≥ 0.32, against a **9.9 mrad** acceptance cone (vs 13.1 GeV and
29 mrad at low x). **The dimuon charge-tag and the large-x signal region are
mutually exclusive.**

A dedicated asymmetric-charm-PDF production was therefore **considered and not
run**: it would change the generated c/c̄ content but not the acceptance, so the
tagged yield above x = 0.32 would still be zero.

**What this does not kill:** the large-x charm *rate* measurement needs no charge
information, hence no decay muon in the spectrometer — §2c stands. The natural
route to large-x charm is the **c→hadrons** channel (the FASERcal ML tagger
already does it, 21 % eff / 82 % purity), which measures the excess but gives no
asymmetry.

---

## 3. Key physics findings

1. **The chain is branching-limited, not detector-limited.** 2.83 % × 17.6 % ≈ 200× attrition before any detector effect.
2. **94 % of charm events contain exactly two charm hadrons** — γg→cc̄. Explains the small truth asymmetry and the ~2× per-event semileptonic rate.
3. **Acceptance is angular, not momentum.** Decay µ median 33 mrad vs scattered µ 4.3 mrad; 83 % already pass the momentum cut.
4. **Tungsten helps, but modestly** — the designed 5 mm option is only ×1.34 in tagged yield over 1 mm, because the absorber is per module (see correction #8).
5. **Charge ID is not the bottleneck** — median decay µ is 11.5 GeV, so η ≈ 3 % for any plausible magnet. The MDT/1.5 T design is comfortably sufficient.
6. **Off-axis is not a penalty, but FASERCal sits in the quieter lobe.** The flux is strongly asymmetric about the LoS (up to ×4.9 on the −x side, where it is also harder and µ⁺-enriched). The detector is on the +x side, where F_flux = 0.96 — the conservative branch.

---

## 4. Corrections made along the way

| # | was | now | impact |
|---|---|---|---|
| 1 | charge confusion from a single MC draw | analytic per-muon expectation | removed noise letting \|A_meas\|>\|A_true\| |
| 2 | A_c ≈ −0.16 (243 muons) | **+0.002 ± 0.020** at 10× stats | the −0.16 was a fluctuation |
| 16 | A_c treated as an IC observable | **zero by construction** in NNPDF4.0 | the study's original premise does not hold |
| 17 | error 1/√N_raw on the tagged asymmetry | **1/√N_eff**, N_eff = 778 not 12 019 | −0.061 is 1.7σ, not 6.7σ |
| 3 | "40 % ≡ θ < 17 mrad" (unweighted) | **8.7 mrad** (yield-weighted) | ~2× tighter aperture |
| 4 | reference mass 1.0 t | **0.724 t** (flux bakes in the W target) | yields ×1.38 |
| 5 | "off-axis is an upper bound" | flux **rises** off-axis | reversed |
| 6 | predicted off-axis spectrum softer | measured **harder** | my prediction was wrong |
| 7 | "2× deficit traced to var2" | **var2 is correct** (25×30 cm) | redirects the investigation |
| **8** | **W between every 1 cm layer** | **W per MODULE (every 20 layers)** | **mass 5–10× too high; 5 mm option 7.11 t → 0.692 t** |
| 9 | fiducial mass assumed 1 t | designed **0.514 t** | yields ×0.51 |
| 10 | luminosity 680 fb⁻¹ (TP print) | **780 fb⁻¹** (FASERCal's own) | yields ×1.15 |
| 11 | mass computed from first principles (514 kg) | **CDR Table 5: 581 kg** | ×1.13 — my calc omitted Al enclosures, WLS fibres, glue, Tyvek |
| 12 | ~~z < 1150 mm read as a containment cut (×0.48)~~ | **no fiducial reduction** — the cut just means "inside the 3DCal", in an older centred frame where it keeps 97.7 % | yields **×2.1** back up |
| 13 | toy charge-confusion sigmoid | **CDR §2.6.2 measured**: 0.7 % below 100 GeV | purity 97.3 % → **99.2 %** — the real spectrometer is *better* than my toy |
| 14 | flux scaled by luminosity alone | **× F_optics = 2.0** (CDR §1.3.2) | yields **×2** |
| 15 | 3DCAL only | **+ AHCAL option** (1.60 t fiducial, ~98 % Fe) | DIS rate **×6.7**, σ(A_c) 0.041 → **0.016** |

---

## 5. What is MISSING / open

**Physics not yet simulated**
- [x] ~~x-reconstruction comparison~~ — **DONE**, see §2c and [`docs/XRECO.md`](XRECO.md). Σ retains ×3.5 (x≥0.2) to ×4.7 (x≥0.4) more large-x signal than lepton-only.
- [ ] Fold in **neutrino energy loss** from charm semileptonic decays — the one unmodelled effect that works *against* the calorimetric methods
- [ ] Replace the assumed σ_θh = 20 mrad with a real FASERcal reconstruction number
- [ ] More seeds: only **12 perturbative charm events above x = 0.4**, so fitted-vs-perturbative ratios there are statistics-limited
- [ ] Intrinsic charm **differentially in x** (inclusive ratio is only 1.08 and washes out)
- [ ] Backgrounds: π/K decay-in-flight muons faking the soft second muon — **not simulated**
- [ ] Pile-up / linking with the ~1 Hz/cm² background muon flux
- [ ] Detector reconstruction — all kinematics still truth-level

**Open questions**
- [x] ~~Run 4 muon flux ~2× above luminosity scaling~~ — **applied** as `F_OPTICS = 2.0`. Residual: the *spectral* hardening is still not modelled (needs a Run 4 FLUKA sample), and it would push yields further up.
- [ ] **Sign convention of X** between the FASERCal CAD and the FLUKA ntuple — worth a **factor 2** in every yield
- [ ] Residual **~1.7×** normalisation gap vs the paper's Table 2.1, after the reference-mass fix. Not the flux variant. Unexplained.
- [x] ~~The talk quotes 4.8 λ for the 3DCAL, my model gave 3.0~~ — **resolved**: the quoted figure is the nuclear **collision** length λ_T, mine was the **interaction** length λ_I. With λ_T the 5 mm option gives 5.00 vs the quoted 4.8 (and X₀ 19.1 vs 18.3). A definition mismatch, not missing material.
- [ ] The 4.5° **tilt** is not yet applied to the acceptance/angle calculation

**Toy values still standing in**
- [ ] Punch-through threshold (2 GeV). The identification+linking efficiency is no longer assumed: the main result is now **plotted against it** across 0.05–1.0 (report page 20), so the answer can be read off once the real number is known.
- [ ] Charge-confusion curve — should be replaced by the MDT/1.5 T GenFit result

---

## 6. WHAT I STILL NEED FROM YOU

Much shorter than before — the talk resolved most of it.

| # | need | why it matters |
|---|---|---|
| **3** | Hadronic energy / jet resolution for **3DCAL alone** — the talk gives them for the full 3DCAL+ECAL+AHCAL+MuSpect chain (p_jet: bias 9 %, res 9 % for NC) | needed for the x-reconstruction study; if AHCAL is excluded the resolution will be worse |
| **4** | Confirmation that muon DIS should use 3DCAL only, or 3DCAL+ECAL+AHCAL | the talk's own numbers show AHCAL has 3.3 t vs 3DCAL's 0.5 t — a **6× larger** target |

---

## 7. What I can do next WITHOUT new input

- x-reconstruction study (all four methods) — `e_in`, `E_had`, hadronic vector and Σ(E−p_z) are **already cached**, no re-showering
- Truth-x → reco-x migration matrix and large-x IC yields
- Apply the 4.5° tilt to the acceptance geometry
- More FLUKA statistics for the off-axis map (3 of ~15 files used)
- π/K decay-in-flight background estimate from the existing showered events
- **Extend the chain to ECAL + AHCAL** — worth doing regardless, since AHCAL's 3.3 t dwarfs the 3DCAL's 0.5 t

---

## 8. Repo map

| file | role |
|---|---|
| `python/shower_dis.py` | Pythia8 → per-event cache (DIS, hadronic, charm, dimuon) |
| `python/fasercal_chain.py` | normalisation, composition, response, the chain, absorber options |
| `python/geometry.py` | 3DCAL as designed: modules, material budget, multiple scattering |
| `python/fluka_offaxis_map.py` | FLUKA → transverse fluence map → F_flux |
| `python/make_chain_report.py` | cache → 10-page PDF report |
| `python/shower_charm.py`, `detector.py`, `make_report.py` | earlier per-charm-hadron feasibility study |
| `docs/TECHNOTE_chain.md` | full method + results note |

**Data (EOS):** `fasercal_dis_v1.npz`, `fasercal_dis_pc.npz`, `fluka_muon_map.npz`, `reports/fasercal_chain.pdf`

**To rescale:** `M_FID_T`, `F_FLUX`, `LUMI_RUN4` in `fasercal_chain.py`; absorber options in `geometry.SCENARIOS`.
