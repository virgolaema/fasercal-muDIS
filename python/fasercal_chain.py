"""
The FASERcal muon-DIS charm chain: yields, efficiencies and the c/c~ tag.

Chain requested:
    N(muon DIS in the FASERcal fiducial volume, no AHCAL)
      -> N(events containing charm)
        -> N(charm events decaying semileptonically to a muon)
          -> N(those muons reaching the magnetised spectrometer)   [~40%]
            -> N(identified, momentum-reconstructed, charge-signed and
                 linked back to the DIS vertex)

Normalisation, stated explicitly because it is the least certain part
--------------------------------------------------------------------
The POWHEG weights carry the muon flux through the LHAPDF "flux as a beam PDF"
trick (arXiv:2506.13889 Eq. 2.1), normalised to the FASERnu setup at Run 3
(250 fb^-1).  Converting that to "FASERcal at Run 4" is a product of explicit
factors, each of which is a separate, checkable assumption:

    N = N_ref(250/fb, FASERnu)
        x (L_RUN4 / 250)          luminosity        -- solid (TP)
        x (M_FID / M_REF)         fiducial mass     -- LINEAR, trivially rescaled
        x F_FLUX                  flux at detector  -- see caveat below
        x [CH vs W recombination] target composition -- done per event

MASS is a pure multiplicative factor, so the whole report can be rescaled to the
real 3DCAL fiducial mass by changing M_FID_T alone.

F_FLUX caveat: the flux set is evaluated with the detector axis at
(x,y) = (1 cm, -3.3 cm) w.r.t. the line of sight, i.e. effectively ON-AXIS.
An off-axis detector does NOT simply see less: the LHC magnetic lattice sweeps
mu+ and mu- in opposite directions, so the flux is strongly asymmetric about the
LoS and FLUKA studies report it RISING by up to an order of magnitude in some
directions beyond ~1 m.  Going off-axis can therefore increase the muon-DIS
yield.  See the F_FLUX block below; the real off-axis fluence is a REQUIRED
INPUT that is not yet available.

Independently, the sibling FASERnu study found the absolute normalisation ~2x
below the published rate.  Correcting M_REF_T (the flux definition bakes in a
25x30 cm^2 x 50 cm tungsten target = 0.724 t, not the 1 t previously assumed)
accounts for a factor 1.38 of that; a residual ~1.7x remains unexplained.
The EFFICIENCIES and FRACTIONS in the chain are ratios and are unaffected by all
of this.
"""
import numpy as np

# ---------------------------------------------------------------- luminosity
LUMI_REF   = 250.0     # fb^-1, the normalisation baked into the weights (Run 3)
# 780 fb^-1 is what the FASERCal team itself quotes for Run 4 nominal
# (Bern CM talk, 15 July 2026, slides 5/9/14/23/24), and matches the TP margin
# note "[FK: update to 780]".  Switched from the TP's printed 680.
LUMI_RUN4  = 780.0     # fb^-1, Run 4 nominal as used by FASERCal
LUMI_RUN4_TP = 680.0   # the value printed in the Run 4 TP v0.01

# ------------------------------------------------------------ target material
# Per-nucleon samples are recombined with the proton fraction of the target.
# Plastic scintillator (polystyrene, C8H8): per monomer 8 C (6p+6n) + 8 H (1p),
#   protons  = 8*6 + 8*1 = 56 ; neutrons = 8*6 = 48 ; nucleons = 104
# Tungsten (W-184): 74 p, 110 n.
COMPOSITION = {
    "scintillator_CH": dict(w_p=56.0 / 104.0, w_n=48.0 / 104.0, rho=1.02),
    "tungsten":        dict(w_p=74.0 / 184.0, w_n=110.0 / 184.0, rho=19.3),
    # Fe-56: 26 protons, 30 neutrons.  Needed for the AHCAL, whose mass is
    # dominated by its steel absorber.
    "iron":            dict(w_p=26.0 / 56.0,  w_n=30.0 / 56.0,   rho=7.87),
    # Pb-207: 82 protons, 125 neutrons.  The ECAL is ~92% lead by mass.
    "lead":            dict(w_p=82.0 / 207.0, w_n=125.0 / 207.0, rho=11.35),
}

# ------------------------------------------------------------------- geometry
# 3DCAL fiducial mass, now the DESIGNED value rather than an assumption
# (FASERCal Bern CM talk, 15 July 2026, slides 4/11/39):
#   10 modules x 20 layers of 1 cm^3 cubes, 48 x 48 cm face
#   -> 0.470 t scintillator, plus 0.044 t (1 mm W/module) or 0.222 t (5 mm)
# The 1 mm option is the baseline used in the FASERCal simulations and costing,
# so it is taken as nominal here.  geometry.configure() computes both.
# CDR Table 5: 581 kg total for the 1 mm baseline; CDR Table 8 restricts rates
# to the fiducial z < 1150 mm (48% of the 2410 mm length), giving 277 kg.
# (The first-principles calculation gave 514 kg, 13% light: it omits the
#  aluminium enclosures, WLS fibres, glue and Tyvek.)
M_FID_T = 0.277    # tonne, 3DCAL 1 mm W, fiducial z < 1150 mm (AHCAL excluded)

# REFERENCE MASS implicit in the flux normalisation.  arXiv:2506.13889 Eq. (2.1)
# defines the flux as f_mu(x_mu) = n_T L_T dN_mu/dx_mu, i.e. the TUNGSTEN TARGET
# IS BAKED IN: a 25 x 30 cm^2 face over L_T = 50 cm of tungsten.
#     25*30*50 cm^3 x 19.3 g/cm^3 = 0.724 t
# (An earlier version of this code assumed 1.0 t, under-predicting all absolute
#  yields by a factor 1.38.  See TECHNOTE_chain.md Sec. 5.2.)
M_REF_T = 25.0 * 30.0 * 50.0 * 19.30 / 1.0e6     # = 0.724 t

# ------------------------------------------------------------------ off-axis
# The flux set is evaluated with the detector axis at (x,y) = (1 cm, -3.3 cm)
# w.r.t. the nominal line of sight (paper Sec. 2), i.e. effectively ON-AXIS.
#
# F_FLUX is the ratio of the muon fluence per cm^2 at the actual detector
# position to that on-axis value.  It is NOT a simple suppression:
#   * the muon flux is strongly ASYMMETRIC about the LoS because the LHC
#     magnetic lattice sweeps mu+ and mu- in opposite directions (paper Fig 2.1
#     shows the resulting large mu/mubar asymmetry);
#   * FLUKA studies for the FPF report the rate RISING by up to an order of
#     magnitude in some directions beyond ~1 m from the LoS, and being
#     substantially higher at ~2 m in the horizontal (bending) plane.
# So going off-axis can INCREASE the muon-DIS yield.  A plausible working range
# is ~0.3 (shielded/vertical) to ~10 (horizontal, beam-pipe side).
#
# RESOLVED, and now CONFIRMED against the FLUKA convention.
# The detector sits on the LEFT of the LoS looking downstream, which in the
# standard right-handed frame (+z downstream, +y up) means +x -- and the FLUKA
# coordinate system has been confirmed to place FASERCal at POSITIVE x (user,
# 2026-08).  The earlier caveat that this might be reversed, worth a factor 2.1,
# is therefore CLOSED.
#
# Position: the CDR/talk shift is (452, 236) mm, plus a further 12 cm reported by
# the collaboration, giving (572, 236) mm.  The FLUKA map over the 480x480 mm
# 3DCAL face there gives
#
#     F_flux = 1.10 ,  <p_mu> = 1188 GeV ,  mu+ fraction = 0.18
#
# The 12 cm shift matters little: across the four possible directions F_flux
# spans only 0.90-1.17, so the conclusion is insensitive to which axis it is on.
#
#     position (mm)      F_flux
#     (452, 236) CDR       0.96
#     (572, 236) +12 x     1.10   <- used
#     (332, 236) -12 x     0.97
#     (452, 356) +12 y     1.17
#     (452, 116) -12 y     0.90
#
# Full transverse map, F_flux for a FASERnu-sized window at (x0,y0) cm:
#        y0=0    y0=+50  y0=+100        <p> [GeV]      mu+ fraction
#   x0=+100  1.49   0.98    1.88          1973            0.09
#   x0=   0  1.00   1.42    1.43           931            0.27
#   x0= -50  1.88   2.99    1.76          1969            0.78
#   x0=-100  3.89   4.94    2.43          2079            0.74
# The flux RISES off-axis on the -x side (magnetic sweeping) and the lobes are
# charge-separated; the +x side, where FASERCal sits, is close to on-axis.
F_FLUX = 1.10                    # measured at (+572,+236) mm

# ------------------------------------------------------- Run 4 beam optics
# The muon flux does NOT simply scale with luminosity from Run 3 to Run 4: the
# LHC optics change too, and they make the muon background substantially worse.
# FASERCal CDR v0 Sec. 1.3.2, verbatim:
#   "For 2022 and 2023 Run 3 physics, the LHC used a 160 urad downward crossing
#    angle. For 2024 running the crossing angle has moved upwards; these reverse
#    polarity optics reduce the radiation incident on magnets close to IP1. In
#    order to match these changes, the Q4 magnet ... has been switched off. The
#    collimators for the outgoing beam also have looser settings as a result and
#    are less effective at removing potential background.
#    An early test of the 2024 configuration suggests an increased rate of
#    background muons at FASER, increasing the rate by a factor of 2 and
#    increasing the energy of the muons. These tests particularly suggest a
#    large increase in high-momentum, positively charged muons. ...
#    The plan for Run 4 is also to have a horizontal crossing angle, increased
#    to 250 urad. Various optics configurations are being tested; however, it is
#    not clear whether the dramatic increase in muon flux can be mitigated."
#
# The FLUKA sample used to build both the LHAPDF flux grid and the off-axis map
# is the Run 3 (MC22) production, so this enhancement is NOT already contained
# in the weights and must be applied explicitly.
#
# F_OPTICS = 2.0 is the CDR's own "factor of 2".  Two caveats, both in the
# direction of this being CONSERVATIVE:
#   * it is quoted for the 2024 configuration, while Run 4 goes further (250 urad
#     horizontal) with mitigation explicitly "not clear";
#   * the spectrum also HARDENS, and both the DIS cross-section and the charm
#     fraction rise with E_mu, so a pure rate rescaling under-counts the gain.
# Set F_OPTICS = 1.0 to recover pure luminosity scaling (the previous behaviour).
F_OPTICS = 2.0
F_FLUX_RANGE = (0.90, 1.17)   # spread over the possible 12 cm shift directions

# --------------------------------------------------- detector response (toy)
# The Run 4 TP Sec. 6.5/6.6 (spectrometer magnet, muon tracker) are empty, so
# these are transparent toy parameters to be scanned, NOT FASERcal specs.
DEFAULTS = dict(
    # a muon must punch through the remaining calorimeter to reach the magnet
    p_punch_gev=2.0,
    # Spectrometer acceptance.  The FASERCal Bern talk (slide 33) reports, for
    # 3DCAL interactions, leading muons crossing >=2 MuSpect stations 43%,
    # >=3 stations 36%, 4 stations 31%.  The user-supplied 40% sits right in
    # that range, so it is CONFIRMED rather than assumed; 0.43 is used as the
    # >=2-station working point.
    acc_spectrometer=0.43,
    # charge-confusion model eta(p) = eta0 + (0.5-eta0)*sigmoid((p-p_half)/w)
    eta0=0.005, p_half_gev=1000.0, p_width_gev=300.0,
    # Probability of correctly linking the spectrometer track back to the DIS
    # vertex through the calorimeter.  REALITY CHECK: the FASERCal ML charm
    # tagger (Bern talk slide 24) achieves, for the c->mu channel in neutrino
    # events, tagging efficiency 0.15 at purity 0.71 -- far below this toy 0.90.
    # Their number folds in the full ML chain on a harder (neutrino) topology,
    # so it is not directly transferable, but it shows 0.90 is optimistic.
    eff_link=0.90,
)


def event_weight(d, material="scintillator_CH", lumi=LUMI_RUN4,
                 m_fid_t=M_FID_T, f_flux=F_FLUX, f_optics=None):
    """
    Per-event expected count in the FASERcal fiducial volume at `lumi`.

        N = w_raw x composition
              x (lumi / 250)      luminosity scaling
              x (m_fid / M_REF)   target mass
              x f_flux            off-axis position (FLUKA map)
              x f_optics          Run 4 beam-optics enhancement (CDR 1.3.2)
    """
    c = COMPOSITION[material]
    f_optics = F_OPTICS if f_optics is None else f_optics
    w_nucl = np.where(d["is_neutron"] > 0.5, c["w_n"], c["w_p"])
    scale = (lumi / LUMI_REF) * (m_fid_t / M_REF_T) * f_flux * f_optics
    return d["w_raw"] * w_nucl * scale


# Charge misidentification, MEASURED.  FASERCal CDR v0 Sec. 2.6.2 reports, from
# a full Geant4 + GenFit Kalman fit of the baseline spectrometer (10 iron planes
# at 1.5 T, 11 SciFi stations at 100 um, tracks required to cross >=8 stations):
#     20 - 100  GeV : misidentification below 0.7%
#    200 - 600  GeV : rises gradually to the 1-2% level
#    800 - 1500 GeV : 2-5%
# Interpolated log-linearly below.  This SUPERSEDES the previous toy sigmoid,
# whose parameters (eta0, p_half_gev, p_width_gev) are now ignored but kept in
# DEFAULTS so old calls do not break.
_CC_P    = np.array([20.,  100.,  200.,  600.,  800., 1500., 5000.])
_CC_ETA  = np.array([0.007, 0.007, 0.010, 0.020, 0.020, 0.050, 0.120])


def charge_confusion(p, *_a, **_kw):
    """Probability of assigning the wrong charge sign, vs muon momentum [GeV]."""
    p = np.asarray(p, dtype=float)
    return np.interp(np.log10(np.clip(p, 1e-3, None)), np.log10(_CC_P), _CC_ETA)


# ------------------------------------------------------- tungsten scenarios
def yield_per_tonne(d, material, lumi=LUMI_RUN4, f_flux=F_FLUX, f_optics=None):
    """Per-event weights for 1 tonne of `material` at `lumi`."""
    return event_weight(d, material=material, lumi=lumi, m_fid_t=1.0,
                        f_flux=f_flux, f_optics=f_optics)


def calibrate_cone(d, geo_base, params=None, rng=None):
    """
    Find the spectrometer angular half-aperture that reproduces the assumed
    acceptance (DEFAULTS['acc_spectrometer'], user-supplied) in the BASELINE
    pure-scintillator geometry.  The same physical cone is then applied to the
    tungsten scenarios, so the comparison isolates the multiple-scattering
    penalty rather than re-assuming the acceptance.
    """
    p = dict(DEFAULTS); p.update(params or {})
    rng = rng or np.random.default_rng(12345)
    th = smeared_theta(d, geo_base, rng)
    w = (yield_per_tonne(d, "scintillator_CH") * geo_base["m_sc"]
         + yield_per_tonne(d, "tungsten") * geo_base["m_w"])
    ok = ((d["n_semilep_mu"] >= 1) & np.isfinite(th)
          & (d["mu2_p"] > p["p_punch_gev"] + geo_base["dedx_mean_mev"] / 1e3))
    # WEIGHTED quantile: the chain sums weights, so the cone must be defined on
    # the same weighted distribution or the baseline will not reproduce the
    # assumed acceptance.
    order = np.argsort(th[ok])
    ts, ws = th[ok][order], w[ok][order]
    cdf = np.cumsum(ws) / ws.sum()
    return float(np.interp(p["acc_spectrometer"], cdf, ts))


def smeared_theta(d, geo, rng):
    """
    Charm-decay muon polar angle after multiple scattering in the detector.
    Small-angle 2D treatment: the true direction (theta,0) picks up an
    independent Gaussian deflection in each transverse plane.
    """
    import geometry as G
    th = d["mu2_theta"].astype(float)
    pmu = d["mu2_p"].astype(float)
    sig_plane = G.theta_ms(pmu, geo["x0_mean"]) / np.sqrt(2.0)
    dx = rng.normal(0.0, 1.0, size=th.shape) * sig_plane
    dy = rng.normal(0.0, 1.0, size=th.shape) * sig_plane
    return np.sqrt((th + dx) ** 2 + dy ** 2)


def chain_scenario(d, t_w_cm, params=None, lumi=LUMI_RUN4, cone=None,
                   rng=None, f_flux=F_FLUX, include_ahcal=False,
                   include_ecal=False, f_optics=None):
    """
    Full chain for a sampling geometry with `t_w_cm` tungsten per layer.

    Yields are built per material: the scintillator mass and the tungsten mass
    each contribute muon-DIS interactions with their OWN proton fraction, then
    are summed.  Detector effects (energy loss, multiple scattering, acceptance,
    charge ID, linking) are common and applied afterwards.
    """
    import geometry as G
    p = dict(DEFAULTS); p.update(params or {})
    rng = rng or np.random.default_rng(2024)
    geo = G.cdr_config(t_w_cm)

    # ---- yields per material, correct composition each ----
    w = (yield_per_tonne(d, "scintillator_CH", lumi, f_flux, f_optics) * geo["m_sc"]
         + yield_per_tonne(d, "tungsten", lumi, f_flux, f_optics) * geo["m_w"])
    # Optionally add the AHCAL as a second target volume.  It is a sampling
    # calorimeter, so interactions in the steel absorber are measured too and
    # the whole mass counts -- ~98% of it iron.  Its coarse 4x4 cm^2 granularity
    # makes vertex finding and muon linking harder, which enters as a lower
    # identification+linking efficiency, not modelled here (see the eff scan).
    if include_ecal:
        e = G.ecal_config()
        w = w + (yield_per_tonne(d, "lead", lumi, f_flux, f_optics) * e["m_pb"]
                 + yield_per_tonne(d, "scintillator_CH", lumi, f_flux, f_optics) * e["m_sc"])
    if include_ahcal:
        a = G.ahcal_config()
        w = w + (yield_per_tonne(d, "iron", lumi, f_flux, f_optics) * a["m_fe"]
                 + yield_per_tonne(d, "scintillator_CH", lumi, f_flux, f_optics) * a["m_sc"])

    n_dis = w.sum()
    is_charm = d["n_charm"] >= 1
    n_charm = w[is_charm].sum()
    has_mu = d["n_semilep_mu"] >= 1
    n_semi = w[has_mu].sum()

    # ---- punch-through: must survive ionisation loss on the way out ----
    dedx_gev = geo["dedx_mean_mev"] / 1e3
    mu_p = np.where(np.isfinite(d["mu2_p"]), d["mu2_p"], -1.0)
    punch = has_mu & (mu_p > p["p_punch_gev"] + dedx_gev)
    n_punch = w[punch].sum()

    # ---- acceptance: MCS-smeared angle against the calibrated cone ----
    th = smeared_theta(d, geo, rng)
    inside = punch & np.isfinite(th) & (th < cone) if cone is not None else punch
    n_acc = w[inside].sum()

    wa = w[inside]
    eta = charge_confusion(mu_p[inside], **p)
    purity_q = (wa * (1 - eta)).sum() / n_acc if n_acc > 0 else 0.0
    n_tag = n_acc * p["eff_link"]

    sgn = d["charm_sign_mu"][inside]
    A_true = (wa * sgn).sum() / n_acc if n_acc > 0 else 0.0
    A_meas = (wa * (1 - 2 * eta) * sgn).sum() / n_acc if n_acc > 0 else 0.0

    # representative MCS width at the median decay-muon momentum
    med_p = float(np.median(mu_p[has_mu & (mu_p > 0)]))
    return dict(
        geo=geo, t_w_cm=t_w_cm, f_flux=f_flux, include_ahcal=include_ahcal,
        include_ecal=include_ecal,
        n_dis=n_dis, n_charm=n_charm, n_semi=n_semi,
        n_punch=n_punch, n_acc=n_acc, n_tag=n_tag,
        f_charm=n_charm / n_dis if n_dis else 0.0,
        f_semi=n_semi / n_charm if n_charm else 0.0,
        f_punch=n_punch / n_semi if n_semi else 0.0,
        f_acc=n_acc / n_punch if n_punch else 0.0,
        f_overall=n_tag / n_dis if n_dis else 0.0,
        purity_q=purity_q, A_true=A_true, A_meas=A_meas,
        n_tag_raw=int(inside.sum()),
        sigma_A=1.0 / np.sqrt(n_tag) if n_tag > 0 else np.inf,
        theta_ms_med_mrad=1e3 * float(G.theta_ms(med_p, geo["x0_mean"])),
        cone_mrad=1e3 * cone if cone else np.nan,
    )


def chain(d, params=None, **kw):
    """
    Run the full chain and return a dict of absolute yields and fractions.
    Charge confusion is taken as its analytic per-muon expectation (no MC draw).
    """
    p = dict(DEFAULTS); p.update(params or {})
    w = event_weight(d, **kw)

    n_dis = w.sum()
    is_charm = d["n_charm"] >= 1
    n_charm = w[is_charm].sum()
    has_mu = d["n_semilep_mu"] >= 1
    n_semi = w[has_mu].sum()

    # muon must punch out of the calorimeter, then be inside the spectrometer
    mu_p = np.where(np.isfinite(d["mu2_p"]), d["mu2_p"], -1.0)
    punch = has_mu & (mu_p > p["p_punch_gev"])
    n_punch = w[punch].sum()
    n_acc = n_punch * p["acc_spectrometer"]

    # charge identification on the accepted muons (weighted mean correct-tag)
    wa = w[punch] * p["acc_spectrometer"]
    eta = charge_confusion(mu_p[punch], **p)
    purity_q = (wa * (1 - eta)).sum() / n_acc if n_acc > 0 else 0.0

    # link the spectrometer track back to the DIS vertex in the calorimeter
    n_tag = n_acc * p["eff_link"]

    # charm asymmetry, per-muon diluted by charge confusion
    sgn = d["charm_sign_mu"][punch]
    A_true = (wa * sgn).sum() / n_acc if n_acc > 0 else 0.0
    A_meas = (wa * (1 - 2 * eta) * sgn).sum() / n_acc if n_acc > 0 else 0.0
    n_tag_raw = int(punch.sum())          # raw MC events, for statistics

    return dict(
        n_dis=n_dis, n_charm=n_charm, n_semi=n_semi,
        n_punch=n_punch, n_acc=n_acc, n_tag=n_tag,
        f_charm=n_charm / n_dis if n_dis else 0.0,
        f_semi=n_semi / n_charm if n_charm else 0.0,
        f_punch=n_punch / n_semi if n_semi else 0.0,
        f_overall=n_tag / n_dis if n_dis else 0.0,
        purity_q=purity_q, A_true=A_true, A_meas=A_meas,
        n_tag_raw=n_tag_raw,
        sigma_A=1.0 / np.sqrt(n_tag) if n_tag > 0 else np.inf,
    )
