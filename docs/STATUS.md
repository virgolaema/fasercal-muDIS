# FASERcal muon-DIS charm study — status

*Updated 2026-08-03, after the FASERCal Bern CM talk (15 July 2026). Repo `muondis-fasercal`, all local, nothing pushed.*

**Where we are:** the full chain from muon-DIS interactions to charge-signed charm
muons is simulated end to end. The Bern talk replaced most of the remaining
assumptions with designed values — and corrected several of them substantially.

---

## 1. Headline numbers

Run 4 at **780 fb⁻¹** (FASERCal's own nominal), **3DCAL only** (AHCAL excluded),
as-designed geometry, on-axis-equivalent flux:

| stage | 1 mm W/module | fraction |
|---|---:|---|
| Muon DIS in 3DCAL | **361 520** | 100 % |
| … containing charm | 10 215 | 2.83 % |
| … charm → semileptonic µ | 1 794 | 17.6 % of charm |
| … µ punches out | 1 477 | 82.4 % |
| … reaches spectrometer (43 %) | 635 | 43 % |
| … identified + linked | **572** | **0.158 % of DIS** |

Charge-tag purity **97.3 %**. Reach, by absorber option and off-axis flux:

| configuration | mass | F=1.0 (+X side) | F=2.0 (−X side) |
|---|---:|---:|---:|
| no W (reference) | 0.470 t | σ(A_c) = 0.044 | 0.031 |
| **1 mm W/module** | **0.514 t** | **0.042** | **0.030** |
| 5 mm W/module | 0.692 t | 0.036 | **0.026** |

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

## 3. Key physics findings

1. **The chain is branching-limited, not detector-limited.** 2.83 % × 17.6 % ≈ 200× attrition before any detector effect.
2. **94 % of charm events contain exactly two charm hadrons** — γg→cc̄. Explains the small truth asymmetry and the ~2× per-event semileptonic rate.
3. **Acceptance is angular, not momentum.** Decay µ median 33 mrad vs scattered µ 4.3 mrad; 83 % already pass the momentum cut.
4. **Tungsten helps, but modestly** — the designed 5 mm option is only ×1.34 in tagged yield over 1 mm, because the absorber is per module (see correction #8).
5. **Charge ID is not the bottleneck** — median decay µ is 11.5 GeV, so η ≈ 3 % for any plausible magnet. The MDT/1.5 T design is comfortably sufficient.
6. **Off-axis is not a penalty.** The flux is asymmetric about the LoS; at the designed position it is ×1.0 or ×2.0 depending on the X sign convention.

---

## 4. Corrections made along the way

| # | was | now | impact |
|---|---|---|---|
| 1 | charge confusion from a single MC draw | analytic per-muon expectation | removed noise letting \|A_meas\|>\|A_true\| |
| 2 | A_c ≈ −0.16 (243 muons) | **+0.002 ± 0.020** at 10× stats | the −0.16 was a fluctuation |
| 3 | "40 % ≡ θ < 17 mrad" (unweighted) | **8.7 mrad** (yield-weighted) | ~2× tighter aperture |
| 4 | reference mass 1.0 t | **0.724 t** (flux bakes in the W target) | yields ×1.38 |
| 5 | "off-axis is an upper bound" | flux **rises** off-axis | reversed |
| 6 | predicted off-axis spectrum softer | measured **harder** | my prediction was wrong |
| 7 | "2× deficit traced to var2" | **var2 is correct** (25×30 cm) | redirects the investigation |
| **8** | **W between every 1 cm layer** | **W per MODULE (every 20 layers)** | **mass 5–10× too high; 5 mm option 7.11 t → 0.692 t** |
| 9 | fiducial mass assumed 1 t | designed **0.514 t** | yields ×0.51 |
| 10 | luminosity 680 fb⁻¹ (TP print) | **780 fb⁻¹** (FASERCal's own) | yields ×1.15 |

---

## 5. What is MISSING / open

**Physics not yet simulated**
- [ ] x-reconstruction comparison (lepton-only vs Jacquet–Blondel vs Σ vs double-angle) and the truth-x → reco-x migration matrix vs arXiv:2506.13889 §4
- [ ] Intrinsic charm **differentially in x** (inclusive ratio is only 1.08 and washes out)
- [ ] Backgrounds: π/K decay-in-flight muons faking the soft second muon — **not simulated**
- [ ] Pile-up / linking with the ~1 Hz/cm² background muon flux
- [ ] Detector reconstruction — all kinematics still truth-level

**Open questions**
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
| **1** | **X sign convention**: does FASERCal's +452 mm shift correspond to +x or −x in the FLUKA ntuple (`truth_prod_x`)? | **factor 2 on every absolute yield** — the single biggest remaining lever |
| **2** | Momentum resolution and charge-confusion vs p from the MDT/GenFit study (slides 35–38 show the machinery, not the resolved curve) | replaces the toy charge-confusion sigmoid |
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
