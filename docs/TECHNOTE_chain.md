# Technical note — muon-DIS charm chain in FASERcal

**Question:** how many muon-DIS interactions occur in the FASERcal fiducial volume
(3DCAL only, AHCAL excluded), how many contain charm, how many of those decay
semileptonically to a muon, how many of those muons reach the magnetised
spectrometer, and how well can we identify, reconstruct and link them?

**Answer in one line:** at Run 4 nominal luminosity, a 1 t scintillator fiducial
mass and **on-axis-equivalent flux**, **≈6.1×10⁵ muon-DIS interactions → ≈900
charge-signed, vertex-linked charm muons**, i.e. **0.15 %** of DIS — giving a
statistical reach on the charm asymmetry of **σ(A_c) ≈ ±0.033**.

> **The off-axis flux is the largest single unknown and is a REQUIRED INPUT
> (§5.2).** It is *not* a suppression factor: the LHC lattice sweeps µ⁺/µ⁻ in
> opposite directions, and FLUKA studies report the muon rate *rising* by up to
> an order of magnitude in some directions beyond ~1 m from the line of sight.
> Yields scale linearly with it and σ(A_c) as 1/√F_flux.

**With a tungsten absorber** (§8.2) the reach improves substantially and the
multiple-scattering penalty is far smaller than expected: **1 mm W/layer →
σ(A_c) = 0.024 (×1.6), 5 mm W/layer → 0.015 (×2.6)**. The absorber choice is
limited by *calorimetric* performance, not by charm statistics.

---

## 1. Detector context — what is and is not known

The analysis is anchored on the **FASER Run 4 Technical Proposal (v0.01)**. That
document is a skeleton, and the sections carrying the numbers this study would
most want are *empty placeholders*:

| TP section | content | status |
|---|---|---|
| §6.2 3DCAL | FASERcal geometry, mass | **empty** |
| §6.5 Spectrometer Magnet | B, B·L, aperture | **empty** |
| §6.6 Muon Tracking Detector | σ_p/p, charge-ID reach | **empty** ("depends on final muon detector technology") |
| §5.1/5.2 | calorimeter performance, scintillator+SiPM | **empty** |

What the TP *does* supply and this note uses:

- **Run 4 nominal luminosity 680 fb⁻¹** (§2.1, Fig. 1). A margin note flags a
  possible update to 780 fb⁻¹; yields scale linearly (×1.15).
- **Table 1** validates the physics case verbatim: the row *"Muon DIS: proton
  structure, intrinsic charm"* requires `d³σ/dE_µ dE_µ' dθ_µ` with additional
  requirements *"Muon momentum, muon angle, hadronic energy, charm tagging"*,
  contributed by **Pixel** and **Spectrometer**.
- **Lepton charge identification** is listed as "event-by-event **with
  spectrometer**" — precisely the mechanism this chain relies on.
- Run 4 background muon flux ≈ 1 Hz/cm² (§5.7.1).

Consequently every detector-response number in this note is a **transparent toy
parameter**, isolated in `DEFAULTS` in [`python/fasercal_chain.py`](../python/fasercal_chain.py),
**not** a FASERcal specification. The two exceptions are supplied externally:
the **40 % spectrometer acceptance** (user-supplied) and the luminosity (TP).

---

## 2. Physics of the tag

In emulsion (FASERν) charm is identified by imaging the sub-mm decay kink. In a
3D scintillator this is impossible: the charm flight length is ~mm, at or below
cell size, and in the interleaved-tungsten configuration the decay happens
inside an absorber slab. **There is no displaced-vertex charm tag in FASERcal.**

The surviving handle is the **semileptonic decay muon**, whose charge *is* the
charm sign:

```
c  → s W⁺ , W⁺ → µ⁺ ν      (charm)      → µ⁺   (PDG id −13)
c̄  → s̄ W⁻ , W⁻ → µ⁻ ν̄      (anticharm)  → µ⁻   (PDG id +13)
```

This yields the **opposite-sign dimuon** signature: the scattered DIS muon
(beam charge, known) plus a second, softer muon of the opposite sign. FASERcal
itself has no field and cannot sign either muon — the sign must come from the
**downstream magnetised spectrometer** (TP §6.5).

The correctness of this identity is *verified in the code*, not assumed: the
extraction asserts `charm_sign_mu == mu2_q` for every muon matched to a charm
ancestor, and it holds for 100 % of events.

---

## 3. Monte Carlo inputs

Reused from the sibling FASERν study ([`generatoroutputanalysis`](../../generatoroutputanalysis)),
so the generator setup is already validated against the muon-DIS paper.

| production | LHAPDF `lhans` | PDF | role |
|---|---|---|---|
| `production_v1` | 331100 | NNPDF4.0 NNLO, **fitted charm** | IC allowed — nominal |
| `production_pc` | 332100 | NNPDF4.0 NNLO **`pch`** | perturbative charm — null hypothesis |

POWHEG configuration: `ih1 13` (muon beam), `fixed_target 1`, `ebeam1 7000`,
`Qmin 1.65 GeV`. Four per-nucleon samples per production —
`{mum,mup} × {proton,neutron}`.

**Statistics used:** 10 distinct seeds × 4 samples × 20 000 events =
**763 848 showered events** per production. The 25 `pythia8_seedN`
directories hold *distinct* POWHEG samples (verified by md5), so multi-seed
running is a genuine statistics increase, not re-hadronisation.

Showering: Pythia 8 via the LCG_104 python bindings, `POWHEG:nFinal = 2`
matching (identical to the FASERν study's `main31.cmnd`).

---

## 4. Extraction algorithm

Implemented in [`python/shower_dis.py`](../python/shower_dis.py); **one row per DIS
event** (the previous `shower_charm.py` wrote one row per charm hadron, which
cannot express "events with ≥1 semileptonic muon").

### 4.1 DIS kinematics
Taken from the Pythia **hard-process record** (`pythia.process`), i.e. LHE
truth, *not* the showered event — the showered record contains decay and
radiation muons that would corrupt the reconstruction. The incoming muon is the
`statusAbs()==21` muon, the scattered muon the first outgoing one.
Q² = −(p_in − p_out)²; θ_µ from the 3-vector opening angle.

### 4.2 Charm hadrons
A charm hadron is any hadron whose PDG code carries a charm quark
(`(|pid|//100)%10 == 4` or `(|pid|//1000)%10 == 4`). Only **weakly-decaying**
charm hadrons are kept — those whose daughters are *not* themselves charmed —
so a D\*→Dπ cascade is counted once, at the D.

### 4.3 Muon assignment — the key step
Every final-state muon is classified by walking its **mother chain**
(`charm_ancestor()`, following both `mother1`/`mother2` so merges are covered):

- a muon **with** a charmed-hadron ancestor → semileptonic charm-decay muon,
  attached to that charm hadron (leading one kept per hadron);
- a muon **without** one → the scattered DIS muon (its chain leads to the beam).

This ancestry test is what makes the dimuon tag well-defined and cleanly
separates the two muons without any kinematic cut.

### 4.4 Hadronic final state
Summed over final-state particles excluding neutrinos and the scattered muon,
recording `E_had`, the momentum vector, and `Σ_h(E−p_z)` (for the Σ-method
x-reconstruction, prepared for the follow-up study).

> **Critical subtlety carried over from the FASERν study:** the muon-beam
> remnant photon (`status 62`, `id 22`) **must** be dropped. The muon flux is
> implemented as a beam PDF of a fictitious 7 TeV beam, so Pythia parks the
> unused balance in a ~5 TeV photon that is not physical. Including it gives
> E_had ≈ 6.5 TeV instead of ≈ 64 GeV. Note this is deliberately *not* a blanket
> status-61–63 cut: the **target-nucleon** remnant is status 63 and is physical.

### 4.5 Weights stored target-agnostic
The four samples are per-nucleon, so the raw POWHEG weight is stored together
with an `is_neutron` flag and the target composition is applied **at analysis
time**. Baking tungsten in — as the FASERν study did — would be wrong for a
scintillator detector (see §6).

---

## 5. Normalisation

The POWHEG weights carry the muon flux through the LHAPDF "flux as beam PDF"
construction (arXiv:2506.13889 Eq. 2.1), normalised to the **FASERν setup at
Run 3 (250 fb⁻¹)**. Converting that to "FASERcal at Run 4" is a product of
explicit, separately-checkable factors:

```
N = N_ref(250 fb⁻¹, FASERν)
    × (L_Run4 / 250)     luminosity        = 680/250 = 2.72   [TP, solid]
    × (M_fid / M_ref)    fiducial mass     = 1.0 t / 1.0 t    [ASSUMPTION]
    × F_flux             flux at detector  = 1.0              [ASSUMPTION]
    × [CH vs W recombination]              per event          [computed]
```

**Mass is a pure multiplicative factor** — the entire report rescales by editing
`M_FID_T` alone.

### 5.1 Is there a Run 4 muon flux in the documents? — cross-check

**Yes, but only a normalisation, not a spectrum.** TP §5.7.1 states the Run 4
background muon flux as **≈1 Hz/cm², ≈500 Hz over the full detector area**.
Two observations:

- 1 Hz/cm² × 500 cm² = 500 Hz, self-consistent; and 500 cm² is exactly the
  25 × 20 cm² pixel active area of §5.3 — so **the TP figure is quoted for the
  ON-AXIS detector**, reinforcing the off-axis caveat below.
- It is a **rate**, with no energy dependence. The TP gives Run 4 *neutrino*
  spectra (Fig. 1) but no Run 4 *muon* spectrum.

**Why using the Run 3 flux shape is nevertheless correct.** The muon flux
*shape* is essentially luminosity-independent: it is set by the LHC optics, the
absorber/rock geometry along the line of sight, and √s — not by how many
collisions per second. Between Run 3 and Run 4 only √s changes (13.6 → 14 TeV),
a few-percent effect on forward production that is not modelled here. What *does*
scale with luminosity is the **normalisation**, and that is applied explicitly
(×680/250). So the treatment is: **Run 3 flux shape + Run 4 luminosity
normalisation**, which is the physically right decomposition.

**Independent normalisation check.** Converting the TP rate to a fluence requires
an assumed instantaneous luminosity (HL-LHC nominal levelled
5×10³⁴ cm⁻²s⁻¹ — external knowledge, not in the TP):

    1 Hz/cm² / 5×10³⁴ cm⁻²s⁻¹  =  2.0×10⁴ muons/cm²/fb⁻¹
                               →  1.36×10⁷ muons/cm² at 680 fb⁻¹

For our 1 m² × 1 m scintillator (areal density 102 g/cm² ≈ 6.1×10²⁵ nucleons/cm²):

| assumed σ_DIS(Q²>1.65) | implied N_DIS | vs our MC (4.53×10⁵) |
|---|---:|---:|
| 50 nb | 4.2×10⁵ | 0.92 |
| 100 nb | 8.4×10⁵ | 1.85 |
| 200 nb | 1.7×10⁶ | 3.69 |

A naive electroweak estimate 4πα²/Q²_min gives ~1.6×10² nb before structure-function
and kinematic suppression, so the realistic σ sits in the 50–100 nb range —
**bracketing our MC value**. The TP's own Run 4 muon rate is therefore
**consistent with our luminosity-scaled Run 3 normalisation to within the O(2)
band already flagged**, and does not indicate a missing factor. It does not
sharpen the normalisation further, because the conversion needs an assumed
L_inst and an assumed σ_DIS.

### 5.2 Off-axis flux, and a reference-mass correction

Reading arXiv:2506.13889 §2 in detail turned up **two things**, one a correction
to this note and one a reversal of a stated caveat.

**(a) The reference mass was wrong — absolute yields were 1.38× too low.**
Eq. (2.1) of the paper defines the flux as

    f_µ(x_µ) ≡ n_T L_T dN_µ/dx_µ(x_µ)

i.e. **the tungsten target is baked into the flux normalisation**: the nucleon
density n_T times L_T = 50 cm, over the 25 × 30 cm² FASERν face. That reference
mass is

    25 × 30 × 50 cm³ × 19.3 g/cm³ = **0.724 t**,  not the 1.0 t assumed here.

Correcting `M_REF_T` raises every absolute yield by **1.382**. It also partly
closes the long-standing normalisation gap with the paper:

| | our MC (1 t W, 250 fb⁻¹) | paper Table 2.1 | gap |
|---|---:|---:|---:|
| before correction | 1.59×10⁵ | 2.7×10⁵ (0.724 t) | **2.34×** |
| after correction | 2.20×10⁵ | 3.73×10⁵ (scaled to 1 t) | **1.69×** |

So **a factor 1.38 of the known "~2× flux discrepancy" was a target-geometry
bookkeeping error, not a flux problem.** A residual ~1.7× remains unexplained
and is still attributed to the flux variant.

**(b) The "off-axis is an upper bound" caveat was wrong — off-axis can be
*higher*.** The paper places the detector axis at **(x, y) = (1 cm, −3.3 cm)**
relative to the nominal line of sight, i.e. the flux set is effectively
**on-axis**. But going off-axis is *not* a monotonic suppression:

- The muon flux is strongly **asymmetric** about the LoS because the LHC magnetic
  lattice sweeps µ⁺ and µ⁻ in **opposite** directions — the paper's Fig. 2.1
  shows exactly this as a large µ/µ̄ asymmetry, with µ⁻ dominating at large x_µ.
- FLUKA studies for the FPF report the muon rate **rising by up to an order of
  magnitude in some directions beyond ~1 m** from the LoS, and being
  substantially higher at ~2 m in the horizontal (bending) plane.

A plausible working range for the off-axis factor is therefore **F_flux ≈ 0.3
(shielded/vertical) to ≈ 10 (horizontal, beam-pipe side)**, spanning a factor 30.
Yields scale linearly with it and σ(A_c) as 1/√F_flux:

| F_flux | baseline tagged | σ(A_c) | 5 mm W tagged | σ(A_c) |
|---|---:|---:|---:|---:|
| ×0.3 | 271 | 0.061 | 1 852 | 0.023 |
| **×1.0 (on-axis-equiv.)** | **902** | **0.033** | **6 175** | **0.013** |
| ×3.0 | 2 712 | 0.019 | 18 524 | 0.007 |
| ×10 | 9 039 | 0.011 | 61 746 | 0.004 |

**REQUIRED INPUT:** the FLUKA muon fluence at the actual 3DCAL position (or the
off-axis flux file used in the FASERcal studies), together with the position
itself relative to the LoS. This is the **single largest unquantified factor in
the study** — larger than the fiducial mass, the residual normalisation gap, and
every detector assumption combined.

> **A second-order consequence worth checking when that flux arrives.** Magnetic
> sweeping is momentum-dependent (deflection ∝ 1/p), so an off-axis position does
> not just rescale the flux — it *reshapes the spectrum*, preferentially
> selecting lower-momentum muons. Since DIS and charm production both rise
> steeply with E_µ, a softer off-axis spectrum would partly offset a higher
> integrated rate. The per-event cache stores `e_in`, so a spectral reweighting
> can be applied directly once the off-axis flux shape is available — no
> re-showering needed.
>
> Relatedly, off-axis in the bending plane enriches **one beam charge**. That is
> harmless for inclusive DIS (γ\* exchange is charge-blind, as the paper notes)
> but is actually *helpful* for the dimuon tag, where a known beam charge makes
> the opposite-sign/same-sign classification cleaner.

### Two normalisation caveats, stated plainly

1. **Off-axis flux — see §5.2.** Superseded: off-axis is *not* a simple
   suppression, and the naive "upper bound" statement was wrong.
2. **Residual ~1.7× normalisation gap** after the reference-mass correction of
   §5.2. Absolute yields carry that uncertainty; fractions do not.

**The efficiencies and fractions in the chain are ratios and are immune to both.**
They are the robust deliverable; absolute yields should be read with the above.

---

## 6. Target composition — scintillator, not tungsten

DIS is per-nucleon, so the proton fraction of the target matters:

| target | protons | neutrons | w_p | ρ [g/cm³] |
|---|---|---|---|---|
| plastic scintillator (polystyrene C₈H₈) | 56 | 48 | **0.538** | 1.02 |
| tungsten (W-184) | 74 | 110 | 0.402 | 19.3 |

Per **unit mass** the DIS rate differs by only **2.3 %** (scintillator slightly
higher, being marginally proton-rich — protons have the larger CC/EM DIS
cross-section here). Per **unit volume** scintillator is ~19× worse, purely from
density. So the fiducial *mass* is the number that matters, and the material
change is a small correction — but it is applied correctly rather than assumed away.

*Cross-check:* re-running with tungsten at 250 fb⁻¹ and 1 t reproduces
1.59×10⁵ DIS interactions, consistent with the FASERν study's normalisation.

---

## 7. Detector response model

All in `DEFAULTS`, [`python/fasercal_chain.py`](../python/fasercal_chain.py):

| parameter | value | origin |
|---|---|---|
| `p_punch_gev` | 2.0 GeV | muon must punch out of the calorimeter to reach the magnet — toy |
| `acc_spectrometer` | **0.40** | **user-supplied** |
| `eta0`, `p_half_gev`, `p_width_gev` | 0.005, 1000, 300 | charge-confusion sigmoid η(p) — toy, scanned |
| `eff_link` | 0.90 | linking the spectrometer track back to the DIS vertex — toy |

Charge confusion is applied as its **analytic per-muon expectation**, not a
single random draw: purity = ⟨1−η⟩_w and A_meas = Σw(1−2η)s/Σw. (A single MC
realisation adds spurious noise and can even make |A_meas| > |A_true|, which the
dilution forbids.)

---

## 8. Results

Run 4, 680 fb⁻¹, 1 t scintillator fiducial, `production_v1` (fitted charm):

| stage | events | fraction |
|---|---:|---|
| Muon DIS in FASERcal fiducial volume | **613 174** | 100 % |
| … containing charm (≥1 charm hadron) | **17 326** | 2.83 % of DIS |
| … charm decaying semileptonically to µ | **3 042** | 17.6 % of charm |
| … µ punches out of calorimeter (p > 2 GeV) | **2 506** | 82.4 % |
| … reaches spectrometer (40 % acceptance) | **1 002** | 40 % |
| … identified + linked to DIS vertex (90 %) | **902** | **0.147 % of DIS** |

*(Normalisation corrected per §5.2(a): ×1.382 relative to the first version of
this note. Fractions are unchanged — they are ratios.)*

| quantity | value |
|---|---|
| charge-tag purity (c vs c̄) | 97.3 % |
| mean charge confusion η | 2.7 % |
| truth charm asymmetry A_c | −0.022 |
| **statistical reach σ(A_c)** | **±0.033** |
| raw MC events behind the tagged sample | 10 002 (statistically solid) |

Supporting distributions:

| quantity | value |
|---|---|
| incoming muon energy | median 665 GeV |
| scattered muon momentum | median 469 GeV |
| scattered muon angle θ_µ | median 4.3 mrad |
| hadronic energy E_had | median 64 GeV |
| **charm hadrons per charm event** | **mean 2.05 — 94 % have exactly 2** |
| decay-muon momentum | median 11.5 GeV, mean 31.7 |
| decay-muon p > 2 / 5 / 10 GeV | 83.2 % / 67.5 % / 53.2 % |
| **decay-muon angle θ** | **median 33 mrad, 90th pct 285 mrad** |

### 8.1 Fitted vs perturbative charm through the same chain

| stage | fitted (IC allowed) | perturbative | ratio |
|---|---:|---:|---:|
| DIS interactions | 443 785 | 445 971 | 0.995 |
| charm events | 12 540 | 11 596 | **1.081** |
| semileptonic µ | 2 202 | 2 058 | 1.070 |
| tagged | 653 | 664 | 0.983 |
| truth A_c | −0.022 | −0.033 | — |

The DIS ratio being 0.995 is a **closure check**: the two productions differ only
in the charm PDF, so the inclusive DIS rate must agree — it does, to 0.5 %.

The inclusive charm excess is only **8 %**, and after the ~200× attrition of the
chain it is statistically invisible (the 0.983 at "tagged" is consistent with
noise on 10⁴ raw MC events). **This is expected and is not a negative result:**
intrinsic charm is a *large-x* effect, and integrating over all x washes it out.
Extracting it requires the differential-in-x analysis of §11, not the inclusive
yield. What this chain establishes is the *sample size* available to do that.

### 8.2 Tungsten absorber scenarios

Two additional configurations were studied: **1 mm** and **5 mm** of tungsten per
layer, interleaved with the scintillator. The comparison is made at a **fixed 1 m
detector envelope** (the trench space constraint, TP §3.1), so tungsten
*displaces* scintillator rather than lengthening the detector — the physically
constrained comparison. Geometry model in
[`python/geometry.py`](../python/geometry.py).

**Material budget** (1 m² × 1 m envelope, 1 cm scintillator per layer):

| scenario | layers | M_scint | M_W | **M_total** | ⟨X₀⟩ | ⟨dE⟩ | λ_int |
|---|---:|---:|---:|---:|---:|---:|---:|
| baseline (no W) | 100 | 1.02 t | — | **1.02 t** | 1.2 | 98 MeV | 1.2 |
| 1 mm W | 90.9 | 0.93 t | 1.75 t | **2.68 t** | 14.1 | 190 MeV | 2.1 |
| 5 mm W | 66.7 | 0.68 t | 6.43 t | **7.11 t** | 48.4 | 434 MeV | 4.2 |

Tungsten's density means even 1 mm plates *dominate* the nucleon count: at 1 mm
the detector is already 65 % tungsten by mass, at 5 mm it is 90 %.

**Chain results** (Run 4, 680 fb⁻¹, acceptance cone calibrated on the baseline):

| scenario | N_DIS | charm | semilep µ | accepted | **tagged** | acceptance | θ_MS | **σ(A_c)** | gain |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 452 660 | 12 790 | 2 246 | 727 | **654** | 39.5 % | 1.9 mrad | **0.039** | ×1.00 |
| 1 mm W | 1 172 393 | 33 646 | 5 866 | 1 888 | **1 699** | 39.3 % | 6.9 mrad | **0.024** | **×1.61** |
| 5 mm W | 3 091 679 | 89 260 | 15 520 | 4 965 | **4 469** | 39.1 % | 13.4 mrad | **0.015** | **×2.61** |

**The mass gain wins outright.** Multiple scattering grows steeply — θ_MS at the
median decay-muon momentum goes 1.9 → 6.9 → 13.4 mrad, exceeding the 8.7 mrad
aperture at 5 mm — yet acceptance falls only 39.5 % → 39.1 %.

The reason is worth stating because it is counter-intuitive: **the decay-muon
angular distribution is intrinsically far broader than the aperture** (median
33 mrad vs an 8.7 mrad cone). Acceptance is therefore "does the muon happen to
point into a small cone", which samples the angular *density near zero*.
Scattering diffuses that density but does not deplete it — it moves nearly as
many muons into the cone as out. Explicit migration at 5 mm W: **5.3 % scattered
out, 2.4 % scattered in**, net −3 % absolute (−7.5 % relative).

**Cost that this toy does *not* penalise.** Scattering destroys the *measurement*
of the decay-muon angle even where it preserves the *count*. That is harmless for
charge tagging and counting (what this note computes) but would matter for any
analysis using the decay-muon direction. By contrast the **scattered DIS muon is
barely affected** — being ~470 GeV, its θ_MS is 0.045/0.169/0.327 mrad, i.e. only
1.1 / 3.9 / 7.6 % of its 4.3 mrad scattering angle. **So tungsten does not
compromise the DIS kinematics**, which is the measurement Table 1 actually asks for.

**Hadronic containment is the real limit.** The interaction-length budget is
1.2 / 2.1 / 4.2 λ_int, all below the ~6 λ typically needed to contain a hadronic
shower. The absorber-thickness scan shows the tagged yield rising
**monotonically** out to 3 cm/layer with no statistical optimum — meaning the
choice is *not* set by charm statistics but by calorimetric performance
(containment, sampling fluctuations, E_had resolution), which this toy does not
model. That is the right place to make the decision.

---

## 9. Key findings

**1. The chain is branching-limited, not detector-limited.** The two largest
losses are physics, not instrumentation: only 2.83 % of DIS makes charm, and
only 17.6 % of charm gives a muon. Together that is a factor ~200 before any
detector effect. No detector improvement can recover it.

**2. 94 % of charm events contain exactly two charm hadrons.** This is the
signature of photon-gluon fusion γg→cc̄, and it explains why the semileptonic
fraction per *event* (17.6 %) is roughly double the ~9 % per *hadron*: two
chances per event. It also explains why the truth asymmetry is small (−0.022) —
pair production is charge-symmetric, and intrinsic charm rides on top as a
small large-x effect.

**3. Acceptance is an *angular* problem, not a momentum one.** The decay muon is
an order of magnitude wider than the scattered muon (median 33 mrad vs 4.3), with
a tail to 285 mrad. 83 % already pass the punch-through momentum cut, so
essentially all of the loss is geometric.

**The assumed 40 % acceptance is equivalent to a θ < 8.7 mrad cone**, a directly
checkable statement against the real §6.5 aperture.

> **Correction.** An earlier version of this note quoted 17 mrad. That was the
> *unweighted* 40th percentile of the MC events. The event weights correlate
> with angle (harder, more forward interactions carry more weight), so the
> yield-weighted value — the correct one for an event count — is **8.7 mrad**,
> a factor ~2 *tighter*, i.e. a more demanding aperture requirement. All
> acceptance calibration in the code now uses the weighted quantile
> (`calibrate_cone()`).

**4. Charge ID is not the bottleneck.** With a median decay-muon momentum of
11.5 GeV, charge confusion is ~3 % for any plausible spectrometer, and the
asymmetry dilution (1−2η) ≈ 0.95. The magnet needs modest performance at
few-tens of GeV, not TeV-scale reach.

**5. Tungsten pays for itself, and the scattering penalty is much smaller than
naively expected.** Adding absorber multiplies the target mass (×2.6 at 1 mm,
×7.0 at 5 mm) and the tagged yield tracks it almost exactly, while acceptance
falls by only a few percent relative. The naive worry — "MCS exceeds the
aperture, so acceptance collapses" — is wrong, because the decay-muon angular
distribution is intrinsically much wider than the aperture and scattering
diffuses rather than depletes it. Crucially, the *scattered DIS muon* (~470 GeV)
is essentially unaffected (θ_MS ≤ 7.6 % of θ_µ even at 5 mm), so the DIS
kinematics that Table 1 actually asks for survive.

**6. σ(A_c) ≈ ±0.04 at nominal assumptions.** Whether that is sufficient depends
on the predicted intrinsic-charm asymmetry, which this note does not yet compute
differentially in x — that is the follow-up (§11).

---

## 10. Caveats

- **Fiducial mass is assumed (1 t).** TP §6.2 is empty. Everything scales
  linearly; see the mass-scan figure.
- **Off-axis flux suppression not applied** — yields are an upper bound.
- **O(2) absolute flux normalisation uncertainty** inherited from the sibling study.
- **Toy response.** Punch-through, linking efficiency and charge confusion are
  illustrative; only the 40 % acceptance and the luminosity are external inputs.
- **No reconstruction simulation.** Muon momenta/angles are truth-level; no
  voxel granularity, no multiple scattering in the calorimeter, no pile-up from
  the ~1 Hz/cm² background muon flux (which is the real threat to *linking*, and
  is currently just a 90 % efficiency factor).
- **Truth-level charm.** No fake/mis-ID background from π/K decay-in-flight
  muons, which are a genuine background to a soft second muon and are not
  simulated here.

## 11. Follow-up

The per-event cache already stores `e_had`, the hadronic momentum vector and
`Σ_h(E−p_z)`, so the x-reconstruction study (lepton-only vs Jacquet–Blondel vs
Σ vs double-angle, and the truth-x → reco-x migration matrix against
arXiv:2506.13889 §4) can be built without re-showering. Note that the
**double-angle method does not close in muon DIS** — it needs a known beam
energy, which the broad LHC muon flux does not provide; the **Σ method** is the
beam-energy-independent alternative and is the one expected to beat lepton-only.

---

## 12. Reproduction

```bash
# heavy: shower both PDF hypotheses (~2 h), writes per-event caches to EOS
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
python3 python/shower_dis.py --prod production_v1 --seeds 1-10 --npz <eos>/fasercal_dis_v1.npz
python3 python/shower_dis.py --prod production_pc --seeds 1-10 --npz <eos>/fasercal_dis_pc.npz

# fast: rebuild the report from cache, iterate on assumptions freely
python3 python/make_chain_report.py --png-dir docs/figures
```

| file | role |
|---|---|
| [`python/shower_dis.py`](../python/shower_dis.py) | Pythia8 → per-event cache (DIS, hadronic, charm, dimuon) |
| [`python/fasercal_chain.py`](../python/fasercal_chain.py) | normalisation, target composition, response model, the chain, tungsten scenarios |
| [`python/geometry.py`](../python/geometry.py) | sampling geometry, material budget, multiple scattering |
| [`python/make_chain_report.py`](../python/make_chain_report.py) | cache → PDF report |
| [`python/config.py`](../python/config.py) | path resolution (`json/ev.json`) |

To rescale to the real 3DCAL: edit `M_FID_T` in `fasercal_chain.py`; for the
tungsten scenarios edit `AREA_CM2`/`LENGTH_CM`/`D_SCINT_CM` in `geometry.py`.
Absorber thicknesses studied live in `geometry.SCENARIOS`.
To use 780 fb⁻¹: set `LUMI_RUN4 = LUMI_RUN4_ALT` (×1.15 on all yields).
