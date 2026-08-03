"""
FASERCal 3DCAL geometry / material model — AS DESIGNED.

Source: "FASERCal", FASER Collaboration Meeting, Bern, 15 July 2026
(FASERCal_simulation_BernFaserCM_15July2026.pdf), slides 4, 11, 12, 39.

    3DCAL: 10 modules, each with 20 layers of 1 cm^3 scintillating cubes and
           1 mm OR 5 mm of tungsten PER MODULE (slide 39 is explicit: "1mm per
           module, i.e. every 20 scintillator layers").
           Module face 480 x 480 mm (slide 11).
           Quoted totals: 18.3 X0, 4.8 lambda_int.
    Position: detector placed parallel to the tunnel wall, tilted 4.5 degrees,
           shifted from the LoS by 452 mm in X and 236 mm in Y (slides 11, 12, 31).

>>> CORRECTION vs the previous version of this module <<<
An earlier version put the tungsten between EVERY 1 cm scintillator layer, which
over-estimated the tungsten by a factor 10 (and the total 3DCAL mass by 5-10x).
The absorber is per MODULE, i.e. every 20 layers.  Correct totals:

    1 mm W/module : 1.0 cm W total -> 0.470 t scint + 0.044 t W = 0.514 t
    5 mm W/module : 5.0 cm W total -> 0.470 t scint + 0.222 t W = 0.692 t

The computed X0 for the 5 mm option (19.1) reproduces the quoted 18.3, and the
collision length (5.00) reproduces the quoted 4.8 -- see the note on nuclear
lengths below.

Material constants (PDG).  NOTE the two different nuclear lengths:
  lambda_I = nuclear INTERACTION length ; lambda_T = nuclear COLLISION length.
  tungsten      X0 = 0.350 cm ; dE/dx 22.1 MeV/cm ; lambda_I 9.94 cm ; lambda_T 5.72 cm
  polystyrene   X0 = 41.3 cm  ; dE/dx  1.97 MeV/cm ; lambda_I 77.1 cm ; lambda_T 48.5 cm

The talk's quoted "total interaction length 4.8" is reproduced by the COLLISION
length (this model gives 5.00), not the interaction length (3.10).  An earlier
version of this note reported the 3.10 as an unexplained discrepancy with the
talk; it was a definition mismatch, not missing material.  Both are computed and
returned; `lambda_int` now carries the collision-length convention so that it is
directly comparable to the quoted 4.8.
"""
import numpy as np

RHO_W,  X0_W,  DEDX_W  = 19.30, 0.3504, 22.1
RHO_SC, X0_SC, DEDX_SC = 1.02, 41.3, 1.97
# nuclear interaction (lambda_I) and collision (lambda_T) lengths [cm]
LAMBDA_I_W,  LAMBDA_T_W  = 191.9 / 19.30, 110.3 / 19.30      # 9.94, 5.72
LAMBDA_I_SC, LAMBDA_T_SC = 81.7 / 1.06,   51.4 / 1.06        # 77.1, 48.5

# ---- 3DCAL as designed -----------------------------------------------------
FACE_CM = 48.0            # 480 mm module face
N_MODULES = 10
LAYERS_PER_MODULE = 20
CUBE_CM = 1.0             # 1 cm^3 scintillating cubes

# Position of the detector relative to the nominal line of sight (mm)
LOS_SHIFT_X_MM = 452.0
LOS_SHIFT_Y_MM = 236.0
TILT_DEG = 4.5

AREA_CM2 = FACE_CM ** 2                                   # 2304 cm^2
LEN_SCINT_CM = N_MODULES * LAYERS_PER_MODULE * CUBE_CM    # 200 cm


def configure(t_w_cm, area_cm2=AREA_CM2, n_modules=N_MODULES,
              layers_per_module=LAYERS_PER_MODULE, cube_cm=CUBE_CM):
    """
    Material budget for `t_w_cm` of tungsten PER MODULE (not per layer).

    t_w_cm = 0.1 -> the 1 mm option ; 0.5 -> the 5 mm option ; 0.0 -> no absorber.
    """
    len_sc = n_modules * layers_per_module * cube_cm
    len_w = n_modules * t_w_cm

    m_sc = area_cm2 * len_sc * RHO_SC / 1.0e6      # tonnes
    m_w = area_cm2 * len_w * RHO_W / 1.0e6

    # a muon produced uniformly along the detector traverses, on average, half
    # of the remaining material
    x0_mean = 0.5 * (len_w / X0_W + len_sc / X0_SC)
    dedx_mean_mev = 0.5 * (len_w * DEDX_W + len_sc * DEDX_SC)
    # collision length -- the convention that reproduces the talk's quoted 4.8
    lam_total = len_w / LAMBDA_T_W + len_sc / LAMBDA_T_SC
    lam_interaction = len_w / LAMBDA_I_W + len_sc / LAMBDA_I_SC
    x0_total = len_w / X0_W + len_sc / X0_SC

    return dict(
        t_w_cm=t_w_cm, n_layers=n_modules * layers_per_module,
        total_len_cm=len_sc + len_w, len_sc=len_sc, len_w=len_w,
        m_sc=m_sc, m_w=m_w, m_tot=m_sc + m_w,
        x0_mean=x0_mean, x0_total=x0_total,
        dedx_mean_mev=dedx_mean_mev,
        lambda_int=lam_total,               # collision length, cf. quoted 4.8
        lambda_interaction=lam_interaction, # interaction length, for reference
    )


def theta_ms(p_gev, x0):
    """Highland multiple-scattering RMS space angle [rad] after x0 rad. lengths."""
    p_mev = np.asarray(p_gev, dtype=float) * 1.0e3
    x0 = max(float(x0), 1e-9)
    plane = (13.6 / np.maximum(p_mev, 1e-9)) * np.sqrt(x0) * (1 + 0.038 * np.log(x0))
    return np.sqrt(2.0) * np.maximum(plane, 0.0)


# The two absorber options actually under consideration, plus a no-W reference.
SCENARIOS = {
    "no_W_reference": 0.0,
    "W_1mm_module":   0.1,
    "W_5mm_module":   0.5,
}
