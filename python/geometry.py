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
LOS_SHIFT_X_MM = 572.0   # CDR/talk 452 mm + a further 12 cm (user, 2026-08)
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
# ---------------------------------------------------------------------------
# CDR Table 5 (v0, 16 Feb 2026) gives the AUTHORITATIVE material budget, which
# supersedes the values computed from first principles above (those omit the
# aluminium enclosures, WLS fibres, glue and Tyvek, and came out 13-29% light).
# Resolution rule for this study: the Bern CM talk (15 Jul 2026) wins on
# conflicts, the CDR fills gaps.  Mass was a gap.
#
#   t_W/module   length[mm]  weight[kg]   X0    lambda   nu-int per ab^-1 [k]
#      1 mm         2410        581       6.9     4.4          15.1
#      5 mm         2450        896      18.3     4.8          19.6
#     10 mm         2500       1118      32.5     5.3          25.2
CDR_TABLE5 = {
    0.1: dict(length_mm=2410, mass_kg=581,  x0=6.9,  lam=4.4, nu_per_ab=15.1),
    0.5: dict(length_mm=2450, mass_kg=896,  x0=18.3, lam=4.8, nu_per_ab=19.6),
    1.0: dict(length_mm=2500, mass_kg=1118, x0=32.5, lam=5.3, nu_per_ab=25.2),
}

# FIDUCIAL REGION.  CDR Table 8 quotes rates "restricted to the fiducial region
# z < 1150 mm", i.e. only the upstream ~48% of the 2.4 m detector, so that the
# hadronic shower is contained downstream of the vertex.  The same containment
# argument applies to muon DIS, so the fiducial fraction is applied here too.
FIDUCIAL_Z_MM = 1150.0


def cdr_config(t_w_cm, fiducial=True):
    """
    Material budget from CDR Table 5, optionally restricted to the fiducial
    z < 1150 mm region.  Returns the same keys as configure() so the two are
    interchangeable; use this one for absolute yields.
    """
    if t_w_cm not in CDR_TABLE5:
        raise KeyError(f"CDR Table 5 has no entry for t_W = {t_w_cm} cm")
    c = CDR_TABLE5[t_w_cm]
    f_fid = min(FIDUCIAL_Z_MM / c["length_mm"], 1.0) if fiducial else 1.0
    g = configure(t_w_cm)                       # for the derived quantities
    m_tot = c["mass_kg"] / 1.0e3 * f_fid
    # split the fiducial mass between scintillator and W in the computed ratio
    frac_w = g["m_w"] / g["m_tot"]
    return dict(g, m_tot=m_tot, m_w=m_tot * frac_w, m_sc=m_tot * (1 - frac_w),
                x0_total=c["x0"], lambda_int=c["lam"], f_fiducial=f_fid,
                mass_full_kg=c["mass_kg"], length_mm=c["length_mm"])


# ---------------------------------------------------------------------------
# AHCAL, as an OPTIONAL additional target volume.
#
# Bern CM talk: ~3.3 t.  CDR Table 3: "~3 tons; 40 layers of 2 cm Fe + 0.3 cm
# scintillator; granularity 4x4 cm^2".  CDR Sec. 2.5.2 quotes 5.5 t for the
# prototype including structure, and a sensitive area of 720 x 720 mm^2.
# Per the resolution rule the talk's 3.3 t is used; it is also what the layer
# structure gives from first principles:
#     Fe    : 40 x 2.0 cm x (72 x 72 cm^2) x 7.87 = 3.26 t
#     scint : 40 x 0.3 cm x (72 x 72 cm^2) x 1.02 = 0.06 t
# i.e. the AHCAL target is ~98% IRON by mass -- hence the "iron" composition.
#
# Being a sampling calorimeter, interactions in the absorber ARE measured, so
# the whole mass is target, not just the scintillator tiles.  The catch is
# granularity: 4x4 cm^2 tiles over 40 layers, versus 1 cm^3 voxels in the
# 3DCAL, so vertex finding and muon-to-vertex linking are much harder.  That
# enters the chain as a LOWER identification+linking efficiency, which is why
# the report shows the result against efficiency rather than assuming one.
AHCAL = dict(
    n_layers=40, fe_cm=2.0, sc_cm=0.3, face_cm=72.0,
    rho_fe=7.87, rho_sc=1.02,
)


# ---------------------------------------------------------------------------
# ECAL, as a further optional target volume.
#
# Bern talk: 770 kg, 40 layers of 0.3 cm Pb + 0.3 cm scintillator.  (CDR Table 3
# instead says 66 layers of 2 mm Pb + 4 mm scintillator with 12x12 cm^2
# granularity; per the resolution rule the talk wins on the mass.)  That layer
# structure reproduces a 72 x 72 cm face -- the same as the AHCAL -- and ~21 X0,
# consistent with the CDR's "about 18 radiation lengths", so the two documents
# are describing the same object.
#
# The ECAL is ~92% LEAD by mass (Pb-207: 82 p, 125 n -> w_p = 0.396).
#
# TWO CAVEATS, both worse than for the AHCAL:
#   * granularity is 12x12 cm^2 -- 9x coarser than the AHCAL and 144x coarser
#     than the 3DCAL -- so vertex finding and muon-to-vertex linking are harder
#     still, i.e. a yet lower identification efficiency;
#   * it is only ~0.7 interaction lengths deep, so it cannot contain a hadronic
#     shower by itself.  It does however sit UPSTREAM of the AHCAL (the layout
#     is 3DCAL -> ECAL -> AHCAL -> spectrometer), so an interaction in the ECAL
#     has the AHCAL's ~5 lambda behind it for containment.  For that reason the
#     ECAL is given full fiducial acceptance while the AHCAL, which is the last
#     calorimeter, keeps the 48% of CDR Table 8.
ECAL = dict(n_layers=40, pb_cm=0.3, sc_cm=0.3, face_cm=72.0,
            rho_pb=11.35, rho_sc=1.02, mass_kg=770.0)


def ecal_config(fiducial_frac=1.0):
    """ECAL target mass, split by material.  See the note above on fiducial_frac."""
    e = ECAL
    f_pb = (e["pb_cm"] * e["rho_pb"]) / (e["pb_cm"] * e["rho_pb"]
                                         + e["sc_cm"] * e["rho_sc"])
    m_tot = e["mass_kg"] / 1.0e3 * fiducial_frac
    return dict(m_pb=m_tot * f_pb, m_sc=m_tot * (1 - f_pb), m_tot=m_tot,
                m_full_t=e["mass_kg"] / 1.0e3, f_pb=f_pb,
                f_fiducial=fiducial_frac)


def ahcal_config(fiducial_frac=0.48):
    """
    AHCAL target mass, split by material.

    `fiducial_frac` defaults to the same 0.48 as the 3DCAL (CDR Table 8) for
    consistency.  NOTE this is an assumption: the AHCAL is only ~5 lambda deep
    in total, so a containment-driven cut could be considerably harsher.
    """
    a = AHCAL
    area = a["face_cm"] ** 2
    len_fe = a["n_layers"] * a["fe_cm"]
    len_sc = a["n_layers"] * a["sc_cm"]
    m_fe = area * len_fe * a["rho_fe"] / 1.0e6 * fiducial_frac
    m_sc = area * len_sc * a["rho_sc"] / 1.0e6 * fiducial_frac
    return dict(m_fe=m_fe, m_sc=m_sc, m_tot=m_fe + m_sc,
                m_full_t=(area * len_fe * a["rho_fe"]
                          + area * len_sc * a["rho_sc"]) / 1.0e6,
                len_fe=len_fe, len_sc=len_sc, f_fiducial=fiducial_frac,
                x0_total=len_fe / 1.757 + len_sc / X0_SC,
                lambda_int=len_fe / (132.1 / 7.87) + len_sc / LAMBDA_T_SC)


SCENARIOS = {
    "W_1mm_module":  0.1,      # CDR/talk baseline
    "W_5mm_module":  0.5,
    "W_10mm_module": 1.0,      # CDR Table 5 third option
}
