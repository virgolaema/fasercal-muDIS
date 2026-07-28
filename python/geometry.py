"""
FASERcal geometry / material model for the tungsten-absorber scenarios.

Scenario: FASERcal is a sampling structure of N layers, each
    [ d_scint cm plastic scintillator | t_W cm tungsten absorber ]
inside a FIXED detector envelope of length L (the trench space constraint,
Run 4 TP Sec. 3.1).  Adding tungsten therefore DISPLACES scintillator rather
than extending the detector -- this is the physically constrained comparison.
The alternative (fixed scintillator, detector grows) is reported as a variant.

Three competing effects as t_W increases:

  (+) TARGET MASS.  Tungsten is ~19x denser, so even thin plates dominate the
      nucleon count -> many more muon-DIS interactions.
  (-) MULTIPLE SCATTERING.  X0(W) = 3.5 mm, so the charm-decay muon is angularly
      smeared on its way out.  Since acceptance into the spectrometer is an
      ANGULAR problem (see TECHNOTE_chain.md Sec. 9), this directly eats the
      signal.
  (-) ENERGY LOSS.  ~22 MeV/cm in W raises the punch-through threshold.

There is therefore an optimum absorber thickness, which this module lets the
analysis find rather than assume.

Material constants (PDG):
  tungsten      X0 = 3.504 g/cm^2 -> 0.350 cm ; dE/dx ~ 22.1 MeV/cm ; lambda_I 9.95 cm
  polystyrene   X0 = 43.79 g/cm^2 -> 41.3 cm  ; dE/dx ~ 1.97 MeV/cm ; lambda_I 80 cm
"""
import numpy as np

RHO_W,     X0_W,     DEDX_W,     LAMBDA_W = 19.30, 0.3504, 22.1, 9.95
RHO_SC,    X0_SC,    DEDX_SC,    LAMBDA_SC = 1.02, 41.3,   1.97, 80.0

# Detector envelope.  TP Sec. 6.2 (3DCAL) is an empty placeholder, so this is an
# explicit assumption chosen so the pure-scintillator case gives ~1 t, matching
# the TP benchmark target mass and the baseline of TECHNOTE_chain.md.
AREA_CM2   = 1.0e4    # 1 m^2 transverse
LENGTH_CM  = 100.0    # 1 m longitudinal envelope
D_SCINT_CM = 1.0      # scintillator thickness per layer


def configure(t_w_cm, area_cm2=AREA_CM2, length_cm=LENGTH_CM,
              d_scint_cm=D_SCINT_CM, fixed_envelope=True):
    """
    Return the material budget for an absorber thickness t_w_cm per layer.

    fixed_envelope=True : total length is fixed, tungsten displaces scintillator
    fixed_envelope=False: scintillator is fixed, the detector grows
    """
    pitch = d_scint_cm + t_w_cm
    if fixed_envelope:
        n_layers = length_cm / pitch
        len_sc = n_layers * d_scint_cm
        len_w = n_layers * t_w_cm
        total_len = length_cm
    else:
        n_layers = length_cm / d_scint_cm
        len_sc = length_cm
        len_w = n_layers * t_w_cm
        total_len = len_sc + len_w

    m_sc = area_cm2 * len_sc * RHO_SC / 1.0e6      # tonnes
    m_w = area_cm2 * len_w * RHO_W / 1.0e6

    # A muon produced uniformly along the detector traverses, on average, half
    # of the remaining material.
    x0_mean = 0.5 * (len_w / X0_W + len_sc / X0_SC)
    dedx_mean_mev = 0.5 * (len_w * DEDX_W + len_sc * DEDX_SC)
    lam_total = len_w / LAMBDA_W + len_sc / LAMBDA_SC

    return dict(
        t_w_cm=t_w_cm, n_layers=n_layers, total_len_cm=total_len,
        len_sc=len_sc, len_w=len_w, m_sc=m_sc, m_w=m_w, m_tot=m_sc + m_w,
        x0_mean=x0_mean, dedx_mean_mev=dedx_mean_mev, lambda_int=lam_total,
    )


def theta_ms(p_gev, x0):
    """
    Highland multiple-scattering RMS plane angle [rad] for momentum p after x0
    radiation lengths.  Returned as the space angle (sqrt(2) x plane angle).
    """
    p_mev = np.asarray(p_gev, dtype=float) * 1.0e3
    x0 = max(float(x0), 1e-9)
    plane = (13.6 / np.maximum(p_mev, 1e-9)) * np.sqrt(x0) * (1 + 0.038 * np.log(x0))
    return np.sqrt(2.0) * np.maximum(plane, 0.0)


SCENARIOS = {
    "baseline_no_W": 0.0,
    "W_1mm":         0.1,
    "W_5mm":         0.5,
}
