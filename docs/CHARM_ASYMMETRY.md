# Is there a c/c̄ asymmetry in the intrinsic-charm PDF? — No.

*Checked 2026-08-03. Decisive for whether A_c is measurable at all.*

## 1. The check

Both PDF sets used in this study, evaluated at Q = 10 GeV:

| x | fitted charm (IC) | | | perturbative | |
|---|---|---|---|---|---|
| | x·c | x·c̄ | **(c−c̄)/(c+c̄)** | x·c | **(c−c̄)/(c+c̄)** |
| 0.01 | 2.043e−1 | 2.047e−1 | **−1.0e−3** | 2.057e−1 | −1.2e−3 |
| 0.10 | 2.711e−2 | 2.741e−2 | **−5.6e−3** | 2.682e−2 | −6.6e−3 |
| 0.30 | 8.635e−3 | 8.669e−3 | **−2.0e−3** | 2.729e−3 | −7.2e−3 |
| 0.50 | 3.616e−3 | 3.619e−3 | **−3.6e−4** | 3.303e−4 | −4.7e−3 |
| 0.70 | 5.190e−4 | 5.191e−4 | **−7.3e−5** | 1.172e−5 | −3.9e−3 |

Direct grid inspection of `NNPDF40_nnlo_as_01180_0000.dat`: over all 2352 grid
points the maximum of |c−c̄|/(|c|+|c̄|) is **6.4×10⁻³**. The number sum rule gives

    ∫(c − c̄) dx / ∫(c + c̄) dx = −1.5×10⁻⁴

**Conclusion: NNPDF4.0 has no charm asymmetry.** It fits the total c + c̄ and sets
c − c̄ = 0; the residual sub-percent wobble is fit/interpolation noise. Note the
"asymmetry" is *larger* in the perturbative set, which is the signature of
numerical noise rather than physics (that set has less charm at large x, so the
same absolute noise is a bigger fraction).

**Therefore `production_v1` cannot generate a parton-level A_c. It is zero by
construction.** This is consistent with the light-cone argument: the minimal
|uudcc̄⟩ Fock state contains one c and one c̄ with identical distributions, so
BHPS intrinsic charm predicts A_c(x) = 0 at every x. A non-zero A_c requires
something beyond the minimal Fock state (e.g. a meson–baryon decomposition), is
model-dependent, constrained by ∫(c−c̄)dx = 0, and expected to be small.

## 2. But the MC *does* produce an asymmetry — from fragmentation

This is the part worth bringing to the discussion. Even with an exactly symmetric
charm PDF, the **observable** (muon-tagged) sample is not symmetric, because
hadronisation in the presence of a target remnant is not charge-symmetric:

| species | N(c) | N(c̄) | ratio | asymmetry | BR(→µ) |
|---|---:|---:|---:|---:|---:|
| **Λ_c** | 13 430 | 1 963 | **6.84** | **+0.745** | **4.8 %** |
| D⁺ | 16 683 | 20 355 | 0.82 | −0.099 | 16.7 % |
| D⁰ | 31 730 | 38 820 | 0.82 | −0.100 | 6.6 % |
| D_s | 6 091 | 5 944 | 1.02 | +0.012 | 8.8 % |

**The leading-baryon effect.** The c quark can pick up two valence quarks from
the target remnant to form Λ_c; the c̄ cannot easily form an anti-baryon. Hence
Λ_c outnumbers Λ̄_c by **6.8 : 1**. Λ_c also has the *lowest* semileptonic
branching ratio of any charm hadron (4.8 %, vs 16.7 % for D⁺). So the c quark is
preferentially sequestered into a state that is **hidden from the muon tag**,
while the c̄ goes into D̄ mesons that are not. The muon-tagged sample is therefore
c̄-enriched.

Measured in the muon-tagged sample:

    A = −0.061 ± 0.036  (1.7σ — not significant, but see below)

> **Correction to an earlier number in this repo.** A first pass quoted
> −0.061 ± 0.009 (6.7σ) using 1/√N_raw with N_raw = 12 019. That error is wrong:
> the event weights are highly non-uniform, and the effective sample size is
> **N_eff = (Σw)²/Σw² = 778**, not 12 019. The correct error is ±0.036. The
> underlying *species* asymmetry (Λ_c/Λ̄_c = 6.8) is nevertheless highly
> significant — it is the propagation to the muon-tagged A that is statistics-limited.

## 3. What this means

1. **A_c is not measurable as an intrinsic-charm observable with the PDF inputs we
   have.** The signal is zero by construction in both sets. Studying it requires a
   PDF set with an explicit charm asymmetry, or a model.
2. **Fragmentation generates a competing asymmetry of the same order as our
   sensitivity.** Our reach is σ(A_c) ≈ 0.04 (3DCAL) to 0.016 (+AHCAL); the
   fragmentation-induced central value is ≈ −0.06. So even with an asymmetric PDF,
   a measured A_c would have to be disentangled from the leading-baryon effect —
   which depends on the hadronisation model, not on the proton's charm content.
3. **The real intrinsic-charm signal is the large-x enhancement, not the asymmetry.**
   Same two sets, Q = 10 GeV:

   | x | fitted / perturbative |
   |---|---:|
   | 0.2 | 1.8× |
   | 0.3 | 3.2× |
   | 0.4 | 5.7× |
   | 0.5 | **10.9×** |
   | 0.7 | **44.3×** |

   That is a large, unambiguous signal — and it is exactly the quantity that
   arXiv:2506.13889 §4 shows is destroyed by muon-momentum smearing. Which makes
   the **x-reconstruction study** (lepton-only vs Jacquet–Blondel vs Σ vs
   double-angle) the decisive question, not A_c.

## 4. Reproduce

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current:$LHAPDF_DATA_PATH
python3 -c "
import lhapdf; p = lhapdf.mkPDF('NNPDF40_nnlo_as_01180', 0)
for x in [0.01,0.1,0.3,0.5]:
    c, cb = p.xfxQ(4,x,10.), p.xfxQ(-4,x,10.)
    print(x, c, cb, (c-cb)/(c+cb))"
```
