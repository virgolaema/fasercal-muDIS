# FASERcal muon-DIS charm study — status

*Updated 2026-08-03, after the FASERCal Bern CM talk (15 Jul 2026) and CDR v0 (16 Feb 2026).*
*Resolution rule: **the talk wins on conflicts, the CDR fills gaps**.*
*Repo: [github.com/virgolaema/fasercal-muDIS](https://github.com/virgolaema/fasercal-muDIS) (pushed).*

**Where we are:** the full chain from muon-DIS interactions to charge-signed charm
muons is simulated end to end. The Bern talk replaced most of the remaining
assumptions with designed values — and corrected several of them substantially.

---

## 1. Headline numbers

Run 4 at **780 fb⁻¹** (FASERCal's own nominal), **3DCAL only** (AHCAL excluded),
as-designed geometry, **1 mm W baseline**, fiducial z < 1150 mm,
F_flux = 0.96 (measured), **F_optics = 2.0** (Run 4 beam optics, see below),
identification+linking efficiency ε = 0.9:

| stage | 3DCAL only | + AHCAL | fraction |
|---|---:|---:|---|
| Muon DIS in fiducial volume | **373 655** | **2 504 165** | 100 % |
| … containing charm | 10 580 | 71 546 | 2.83 % |
| … charm → semileptonic µ | 1 856 | 12 551 | 17.6 % of charm |
| … µ punches out | 1 529 | 10 341 | 82.4 % |
| … reaches spectrometer (43 %) | 657 | 4 447 | 43 % |
| … identified + linked (ε = 0.9) | **590** | **3 992** | **0.158 % of DIS** |
| **σ(A_c)** | **0.041** | **0.016** | |

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

| absorber | 3DCAL fiducial | σ(A_c) 3DCAL | σ(A_c) +AHCAL |
|---|---:|---:|---:|
| **1 mm W** (baseline) | 277 kg | 0.041 | **0.016** |
| 5 mm W | 421 kg | 0.034 | 0.015 |
| 10 mm W | 514 kg | 0.030 | 0.015 |

With the AHCAL included the absorber choice becomes almost irrelevant — the
AHCAL mass dominates.

---

## 2. What is DONE

**Simulation chain**
- [x] Per-event extraction from POWHEG+Pythia8 (`shower_dis.py`): DIS kinematics, hadronic final state, charm hadrons, semileptonic dimuon
- [x] Muon classification by **mother-chain ancestry** (charm-decay vs scattered)
- [x] Charm-sign ≡ decay-muon-charge identity **verified in code**, holds 100 %
- [x] 763 848 showered events per PDF hypothesis (10 seeds × 4 samples × 20 k)
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
A = −0.061 ± 0.036.

**Consequences:** (i) A_c is not measurable as an IC observable with current PDF
inputs — it needs an asymmetric-charm set or a model; (ii) the fragmentation
asymmetry is *the same size as our sensitivity*, so it is a competing background,
not a small correction; (iii) **the real IC signal is the large-x enhancement**
(fitted/perturbative = 5.7× at x = 0.4, 10.9× at 0.5, 44× at 0.7), which makes
the x-reconstruction study the decisive question.

---

## 2c. The σ_p inputs are now known — and they threaten lepton-only x reconstruction

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
| 12 | no fiducial cut | **z < 1150 mm (CDR Table 8)** | **×0.48** — the dominant new effect |
| 13 | toy charge-confusion sigmoid | **CDR §2.6.2 measured**: 0.7 % below 100 GeV | purity 97.3 % → **99.2 %** — the real spectrometer is *better* than my toy |
| 14 | flux scaled by luminosity alone | **× F_optics = 2.0** (CDR §1.3.2) | yields **×2** |
| 15 | 3DCAL only | **+ AHCAL option** (1.60 t fiducial, ~98 % Fe) | DIS rate **×6.7**, σ(A_c) 0.041 → **0.016** |

---

## 5. What is MISSING / open

**Physics not yet simulated**
- [ ] **[TOP PRIORITY]** x-reconstruction comparison (lepton-only vs Jacquet–Blondel vs Σ vs double-angle) and the truth-x → reco-x migration matrix vs arXiv:2506.13889 §4, at the CDR's real σ_p (47 %, not the paper's 10/30 %) and σ_E_had (~9 %). See §2c — this is the answer to the showstopper question.
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
| **1** | One-line confirmation that the FLUKA ntuple's `truth_prod_x` is the standard right-handed FASER/ATLAS x (so that "left looking downstream" = +x) | if reversed, every yield ×2.1 |
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
