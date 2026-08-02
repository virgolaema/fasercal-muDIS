# FASERcal muon-DIS charm study — status

*Last updated: 2026-07-29. Repo: `muondis-fasercal` (branch `master`, all local, nothing pushed).*

**Where we are:** the full chain from muon-DIS interactions to charge-signed charm
muons is simulated end to end, with the tungsten-absorber options and the
off-axis muon flux both now measured rather than assumed. What is missing is
almost entirely **detector input from the collaboration**, not simulation work.

---

## 1. Headline numbers

Run 4 (680 fb⁻¹), 1 t scintillator fiducial, on-axis-equivalent flux:

| stage | events | fraction |
|---|---:|---|
| Muon DIS in FASERcal fiducial volume | **613 174** | 100 % |
| … containing charm | 17 326 | 2.83 % |
| … charm → semileptonic µ | 3 042 | 17.6 % of charm |
| … µ punches out of calorimeter | 2 506 | 82.4 % |
| … reaches spectrometer (40 %) | 1 002 | 40 % |
| … identified + linked | **902** | **0.147 % of DIS** |

Charge-tag purity **97.3 %**. Statistical reach:

| configuration | on-axis | best measured off-axis (−100,+50) cm |
|---|---:|---:|
| baseline (no W) | σ(A_c) = 0.033 | **0.015** |
| 1 mm W / layer | 0.021 | 0.009 |
| 5 mm W / layer | 0.013 | **0.006** |

---

## 2. What is DONE

**Simulation chain**
- [x] Per-event extraction from POWHEG+Pythia8 (`shower_dis.py`): DIS kinematics, hadronic final state, charm hadrons, semileptonic dimuon
- [x] Muon classification by **mother-chain ancestry** (charm-decay vs scattered) — no kinematic cuts needed
- [x] Charm-sign ≡ decay-muon-charge identity **verified in code**, holds 100 %
- [x] 763 848 showered events per PDF hypothesis (10 distinct seeds × 4 samples × 20 k)
- [x] Both PDF hypotheses: fitted charm (`production_v1`) and perturbative charm (`production_pc`)

**Normalisation**
- [x] Target composition done properly: scintillator CH (w_p = 0.538) ≠ tungsten (0.402); weights stored target-agnostic so any material can be applied at analysis time
- [x] Reference mass derived from the flux definition (0.724 t, see §4)
- [x] Run 4 luminosity from the TP; cross-checked against the TP's own Run 4 muon rate (§5.7.1) — consistent
- [x] Tungsten cross-check reproduces the FASERν normalisation

**Detector scenarios**
- [x] Tungsten absorber options (1 mm, 5 mm/layer) at fixed 1 m envelope, with material budget, ionisation loss and Highland multiple scattering (`geometry.py`)
- [x] Absorber-thickness optimisation scan
- [x] Angular acceptance mapped to an equivalent aperture cone
- [x] Charge confusion applied analytically (no MC noise)

**Off-axis flux — measured from FLUKA**
- [x] Located the FLUKA source on EOS and built the transverse fluence map (`fluka_offaxis_map.py`)
- [x] F_flux, ⟨p⟩ and µ⁺ fraction tabulated across the transverse plane

**Deliverables**
- [x] 10-page PDF report (`/eos/home-e/evilla/faser/reports/fasercal_chain.pdf`)
- [x] Detailed tech note (`docs/TECHNOTE_chain.md`)
- [x] Earlier feasibility study (`docs/REPORT.md`)

---

## 3. Key physics findings

1. **The chain is branching-limited, not detector-limited.** 2.83 % × 17.6 % ≈ 200× attrition before any detector effect. No instrumentation recovers it.
2. **94 % of charm events contain exactly two charm hadrons** — photon-gluon fusion γg→cc̄. Explains the small truth asymmetry and the ~2× per-event semileptonic rate.
3. **Acceptance is angular, not momentum.** Decay µ median 33 mrad vs scattered µ 4.3 mrad; 83 % already pass the momentum cut. The assumed 40 % ≡ a **θ < 8.7 mrad** cone.
4. **Tungsten pays for itself.** ×2.6 (1 mm) to ×7.0 (5 mm) target mass; the multiple-scattering penalty is only a few % because the decay-muon angular spread is far wider than the aperture, so scattering diffuses rather than depletes. The scattered DIS muon (~470 GeV) is essentially unaffected, so DIS kinematics survive.
5. **Absorber thickness is limited by calorimetry, not statistics** — yield rises monotonically to 3 cm/layer; containment is 1.2/2.1/4.2 λ_int, all below the ~6 λ needed.
6. **Off-axis is a bonus, not a penalty.** Flux rises up to ×4.9, the spectrum gets *harder* (⟨p⟩ 931 → ~2000 GeV) which compounds the gain, and the µ⁺/µ⁻ lobes are charge-separated (µ⁺ fraction 0.27 → 0.75), which helps the dimuon tag.
7. **Charge ID is not the bottleneck** — median decay µ is 11.5 GeV, so η ≈ 3 % for any plausible magnet.

---

## 4. Corrections made along the way

Recorded because several changed the numbers:

| # | was | now | impact |
|---|---|---|---|
| 1 | charge confusion from a single MC draw | analytic per-muon expectation | removed noise that let \|A_meas\| > \|A_true\| |
| 2 | A_c ≈ −0.16 (single seed, 243 muons) | **+0.002 ± 0.020** at 10× stats | the −0.16 was a fluctuation |
| 3 | "40 % ≡ θ < 17 mrad" (unweighted) | **8.7 mrad** (yield-weighted) | ~2× *tighter* aperture requirement |
| 4 | reference mass 1.0 t | **0.724 t** (flux bakes in 25×30 cm² × 50 cm W) | all yields **×1.38** |
| 5 | "off-axis yields are an upper bound" | flux **rises** off-axis, up to ×4.9 | reversed |
| 6 | predicted off-axis spectrum *softer* | measured **harder** (931 → 2000 GeV) | my prediction was wrong |
| 7 | "2× deficit traced to var2 flux variant" | **var2 is correct** (= 25×30 cm, the paper's geometry); var1 is r<9 cm | redirects the investigation |

---

## 5. What is MISSING / open

**Physics not yet simulated**
- [ ] x-reconstruction comparison (lepton-only vs Jacquet–Blondel vs Σ vs double-angle) and the truth-x → reco-x migration matrix, vs arXiv:2506.13889 §4
- [ ] Intrinsic charm **differentially in x** (inclusive ratio is only 1.08 and washes out; the IC signal lives at large x)
- [ ] Backgrounds: π/K decay-in-flight muons faking the soft second muon — **not simulated at all**
- [ ] Pile-up / linking with the ~1 Hz/cm² background muon flux (currently just a 90 % efficiency factor)
- [ ] Detector reconstruction: voxel granularity, hadronic energy resolution, spectrometer momentum resolution — all truth-level now

**Normalisation**
- [ ] Residual **~1.7×** gap vs the paper's Table 2.1, after the reference-mass fix. Not the flux variant (see correction #7). Unexplained.

**Assumptions still standing in for real numbers**
- [ ] Fiducial mass (1 t) — my assumption, taken from the TP's *on-axis tungsten* benchmark, not the 3DCAL
- [ ] Detector envelope (1 m² × 1 m) in `geometry.py` — my construction
- [ ] Punch-through threshold, linking efficiency, charge-confusion curve — toy values

---

## 6. WHAT I NEED FROM YOU

Ordered by how much it changes the answer.

| # | need | why it matters | blocks |
|---|---|---|---|
| **1** | **3DCAL position relative to the LoS** (x, y in cm) | picks the row of the measured F_flux table — worth up to ×4.9 in yield | all absolute numbers |
| **2** | **3DCAL fiducial mass / geometry** (TP §6.2 is empty) | absolute yields scale **linearly**; also fixes the envelope for the tungsten scenarios | all absolute numbers |
| **3** | **Spectrometer**: B·L (or B and lever arm), aperture, σ_p/p, charge-ID reach vs p (TP §6.5/§6.6 empty) | replaces the toy acceptance + charge confusion; lets me check the 40 % against the measured 8.7 mrad equivalent | the A_c reach |
| **4** | **Calorimeter**: hadronic σ_E/E, cell/voxel size, angular resolution (TP §5.1/§5.2 empty) | required for the x-reconstruction study (JB and Σ both need E_had) | the whole IC / large-x programme |
| **5** | **The CDR** — never arrived (two failed attachments). Save it to `/afs/cern.ch/work/e/evilla/private/faser/` and give me the path | may supply 2–4 above | — |
| **6** | **Confirm on-axis vs off-axis for muon DIS.** TP Table 1 checkmarks muon DIS under **Pixel + Spectrometer** (on-axis), but you chose the off-axis 3DCAL | changes geometry, flux and acceptance | scenario definition |
| **7** | Is the **40 % acceptance** inclusive of material effects, or purely geometric? | I calibrate the aperture cone on it | acceptance modelling |
| **8** | Any existing **FASERcal off-axis flux file** used in the collaboration's studies | would supersede my FLUKA map and let me cross-check it | validation |

---

## 7. What I can do next WITHOUT any new input

- x-reconstruction study (all four methods) — `e_in`, `E_had`, hadronic vector and Σ(E−p_z) are **already cached**, no re-showering needed
- Truth-x → reco-x migration matrix and large-x IC yields
- More FLUKA statistics for the off-axis map (currently 3 of ~15 files; 585 raw muons on-axis)
- π/K decay-in-flight background estimate from the existing showered events
- Spectral reweighting to the measured off-axis (harder) spectrum — currently only the rate is applied, not the shape

---

## 8. Repo map

| file | role |
|---|---|
| `python/shower_dis.py` | Pythia8 → per-event cache (DIS, hadronic, charm, dimuon) |
| `python/fasercal_chain.py` | normalisation, composition, response model, the chain, tungsten scenarios |
| `python/geometry.py` | sampling geometry, material budget, multiple scattering |
| `python/fluka_offaxis_map.py` | FLUKA → transverse fluence map → F_flux |
| `python/make_chain_report.py` | cache → 10-page PDF report |
| `python/shower_charm.py`, `detector.py`, `make_report.py` | earlier per-charm-hadron feasibility study |
| `docs/TECHNOTE_chain.md` | full method + results note |
| `docs/REPORT.md` | earlier tagging feasibility report |

**Data (EOS):** `fasercal_dis_v1.npz`, `fasercal_dis_pc.npz`, `fluka_muon_map.npz`, `reports/fasercal_chain.pdf`

**To rescale:** `M_FID_T` (mass), `F_FLUX` (off-axis), `LUMI_RUN4` (780 fb⁻¹ variant) — all in `fasercal_chain.py`.
