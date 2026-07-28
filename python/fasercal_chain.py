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

F_FLUX caveat: the flux set describes the on-axis line-of-sight muon flux.  An
OFF-AXIS detector sees a reduced flux, and that suppression factor is NOT in the
Run 4 TP (Sec. 6 is an empty placeholder in v0.01).  F_FLUX is therefore left at
1.0 = "on-axis-equivalent flux per unit area", and every absolute yield in this
report must be read as an UPPER BOUND for an off-axis placement until the real
off-axis flux is supplied.

Independently, the sibling FASERnu study found the absolute flux normalisation
to sit ~2x below the published rate (unresolved, traced to the flux variant).
Absolute yields therefore carry an O(2) normalisation uncertainty.  The
EFFICIENCIES and FRACTIONS in the chain are ratios and are NOT affected.
"""
import numpy as np

# ---------------------------------------------------------------- luminosity
LUMI_REF   = 250.0     # fb^-1, the normalisation baked into the weights (Run 3)
LUMI_RUN4  = 680.0     # fb^-1, nominal Run 4 (FASER Run 4 TP, Fig. 1 / Sec. 2.1)
LUMI_RUN4_ALT = 780.0  # TP margin note "[FK: update to 780]" -- reported as a variant

# ------------------------------------------------------------ target material
# Per-nucleon samples are recombined with the proton fraction of the target.
# Plastic scintillator (polystyrene, C8H8): per monomer 8 C (6p+6n) + 8 H (1p),
#   protons  = 8*6 + 8*1 = 56 ; neutrons = 8*6 = 48 ; nucleons = 104
# Tungsten (W-184): 74 p, 110 n.
COMPOSITION = {
    "scintillator_CH": dict(w_p=56.0 / 104.0, w_n=48.0 / 104.0, rho=1.02),
    "tungsten":        dict(w_p=74.0 / 184.0, w_n=110.0 / 184.0, rho=19.3),
}

# ------------------------------------------------------------------- geometry
# NOMINAL fiducial mass.  Sec. 6.2 (3DCAL) of the Run 4 TP is an empty
# placeholder, so this is an explicit assumption, chosen equal to the TP
# benchmark target mass (1 t) so the numbers are directly comparable.
# EVERY absolute yield scales linearly with this number.
M_FID_T   = 1.0    # tonne, FASERcal 3DCAL fiducial volume (AHCAL excluded)
M_REF_T   = 1.0    # tonne, reference mass implicit in the flux normalisation
F_FLUX    = 1.0    # on-axis-equivalent flux; see module docstring

# --------------------------------------------------- detector response (toy)
# The Run 4 TP Sec. 6.5/6.6 (spectrometer magnet, muon tracker) are empty, so
# these are transparent toy parameters to be scanned, NOT FASERcal specs.
DEFAULTS = dict(
    # a muon must punch through the remaining calorimeter to reach the magnet
    p_punch_gev=2.0,
    # geometric+tracking acceptance of the spectrometer: USER-SUPPLIED 40%
    acc_spectrometer=0.40,
    # charge-confusion model eta(p) = eta0 + (0.5-eta0)*sigmoid((p-p_half)/w)
    eta0=0.005, p_half_gev=1000.0, p_width_gev=300.0,
    # probability of correctly linking the spectrometer track back to the DIS
    # vertex through the calorimeter (track-following through 3D voxels)
    eff_link=0.90,
)


def event_weight(d, material="scintillator_CH", lumi=LUMI_RUN4,
                 m_fid_t=M_FID_T, f_flux=F_FLUX):
    """Per-event expected count in the FASERcal fiducial volume at `lumi`."""
    c = COMPOSITION[material]
    w_nucl = np.where(d["is_neutron"] > 0.5, c["w_n"], c["w_p"])
    scale = (lumi / LUMI_REF) * (m_fid_t / M_REF_T) * f_flux
    return d["w_raw"] * w_nucl * scale


def charge_confusion(p, eta0, p_half_gev, p_width_gev, **_):
    return eta0 + (0.5 - eta0) / (1.0 + np.exp(-(p - p_half_gev) / p_width_gev))


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
