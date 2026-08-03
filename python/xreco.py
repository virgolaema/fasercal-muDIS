"""
Bjorken-x reconstruction in FASERcal: four methods under realistic smearing.

THE QUESTION (arXiv:2506.13889 Sec. 4, and the "showstopper" it raises):
the intrinsic-charm signal is a large-x excess of charm -- fitted/perturbative
= 1.8x at x=0.2 rising to 44x at x=0.7.  The paper reconstructs x from the muon
alone, and finds that a muon-momentum smearing of 10% already drops the x>=0.2
excess from x8 to ~x2, and that 30% removes it entirely.

FASERcal's own muon spectrometer is iron-core and multiple-scattering dominated,
giving sigma_p/p ~ 47% at the median scattered-muon momentum (CDR Sec. 2.6.2).
That is far inside the regime where the paper says the signal is gone.

BUT FASERcal is a calorimeter, which FASERnu is not.  The hadronic energy
resolution is ~9%, five times better than the muon momentum.  Methods that
reconstruct x from the HADRONIC system should therefore survive where the
lepton-only method dies.  That is the hypothesis this module tests.

THE FOUR METHODS (fixed-target kinematics, masses neglected, M = nucleon mass;
nu = E_in - E_mu' is the energy transfer, x = Q^2 / (2 M nu)):

  1. LEPTON-ONLY  -- what the paper uses.
        Q^2 = 2 E_in E' (1 - cos th_mu)   ;   nu = E_in - E'
     Needs BOTH muon energies.  nu is a small difference of two large numbers,
     so a fractional error on E blows up as (E/nu) x sigma_E -- this is the
     mechanism that kills it.

  2. JACQUET-BLONDEL -- hadronic only.
        nu_JB = E_had   ;   Q^2_JB = pT_had^2 / (1 - y)   ;   y = nu_JB / E_in
     nu comes DIRECTLY from the calorimeter, not from a difference, so a 9%
     hadronic resolution gives a 9% error on nu rather than an amplified one.
     The muon momentum is not used at all (E_in enters only the mild (1-y)).

  3. SIGMA (fixed-target form) -- reconstructs the incoming energy from the
     final state, so it does not need the beam energy measured.
     With the beam along +z, the initial state has E + p_z = 2 E_in + M, hence
        2 E_in = Sum_had(E + p_z) + E'(1 + cos th_mu) - M
     then proceed as lepton-only with that E_in.  (Note: the collider Sigma
     method uses E - p_z; for a fixed target with the beam along +z that
     combination is identically M and carries no information, so the E + p_z
     form is the correct analogue.)

  4. DOUBLE-ANGLE -- angles only, plus the beam energy.
        E'_DA = E_in sin(th_h) / sin(th_mu + th_h)
     Immune to energy scale errors but NOT to the beam energy, which in muon
     DIS is not known event-by-event (unlike a collider) and must be measured.

SMEARING (all switchable; see RESOLUTIONS):
  * muon momentum : CDR Sec. 2.6.2, interpolated in log p  (20% at 20 GeV
    rising to 63% at 1 TeV) -- applied to BOTH muon legs.
  * hadronic energy : ~9% (Bern CM talk, p_jet NC).
  * angles : the 3DCAL measures a track with 1 cm voxels over ~200 layers, so
    the straight-line fit error is sigma_hit sqrt(12/N)/L ~ 0.3 mrad, small
    against the ~4.3 mrad median scattering angle.  The naive "cm voxels are
    much worse than um emulsion" worry is largely defused by the number of
    samples; multiple scattering adds ~0.06 mrad at 470 GeV.
"""
import numpy as np

M_N = 0.9383           # nucleon mass [GeV]

# ---- CDR Sec. 2.6.2 muon-momentum resolution -------------------------------
_P_NODES  = np.array([20., 100., 200., 1000.])
_SIG_NODES = np.array([0.20, 0.228, 0.29, 0.631])

RESOLUTIONS = dict(
    sigma_had=0.09,        # hadronic energy, Bern talk (p_jet, NC)
    sigma_theta=3.0e-4,    # 0.3 mrad, 3DCAL track fit (the MUON)
    # Direction of the HADRONIC system.  This is a shower axis in a calorimeter,
    # not a track, so it is far worse than the muon: 3DCAL cubes are 1 cm but
    # the AHCAL is 4x4 cm, and the axis is diluted by shower fluctuations.
    # 20 mrad is a placeholder against a median theta_h of ~81 mrad; it is
    # scanned in the report because both the double-angle and (through pT) the
    # Jacquet-Blondel methods depend on it.  NOT taken from any document.
    sigma_theta_had=20.0e-3,
    muon_scale=1.0,        # multiply the CDR muon resolution (1.0 = as measured)
)


def sigma_p_over_p(p):
    """CDR muon-momentum resolution, interpolated in log p."""
    p = np.asarray(p, dtype=float)
    return np.interp(np.log10(np.clip(p, 1.0, None)),
                     np.log10(_P_NODES), _SIG_NODES)


def truth(d):
    """Truth-level nu, x, y from the cached DIS kinematics."""
    nu = d["e_in"] - d["p_out"]
    x = d["q2"] / (2.0 * M_N * np.maximum(nu, 1e-9))
    return dict(nu=nu, x=x, y=nu / np.maximum(d["e_in"], 1e-9), q2=d["q2"])


def hadronic_truth(d):
    """
    The TRUE hadronic four-vector, from momentum conservation rather than from
    summing Pythia particles.

    WHY.  The hadronic system recoiling against the scattered muon is exactly
        X = q + p_target ,   q = p_in - p_out ,   p_target = (M, 0, 0, 0)
    so it is fully determined by the muon kinematics:
        E_had  = nu + M
        pT_had = p_out sin(theta_mu)
        pz_had = |p_in| - p_out cos(theta_mu)

    Summing the Pythia final state instead gives a contaminated answer: the muon
    flux is implemented as a beam PDF of a fictitious 7 TeV beam, and that beam's
    remnant deposits a few GeV of hadronic activity that is not part of the DIS
    vertex.  It is negligible at large nu but dominates at large x (where nu is
    small), and it cannot be removed by particle-level cuts -- beam-remnant
    removal, forward-cone cuts and mother-chain ancestry were all tried and all
    failed, the last because Pythia's colour reconnection leaves the ancestry
    graph too connected to separate the two beam sides.  See docs/XRECO.md.

    Using conservation is exact and makes the perfect-detector closure exact by
    construction for every method, which is what a resolution study needs.

    NOT MODELLED, and flagged as such:
      * energy carried off by neutrinos from charm semileptonic decays, which
        biases the measurable E_had low for exactly the signal events;
      * the charged/neutral composition of the shower and its detector response,
        which is folded into the single sigma_had parameter instead.
    """
    nu = d["e_in"] - d["p_out"]
    return dict(
        e=nu + M_N,
        pt=d["p_out"] * np.sin(d["theta_mu"]),
        pz=d["e_in"] - d["p_out"] * np.cos(d["theta_mu"]),
    )


def smear(d, rng, res=None):
    """Apply the detector response to every measured quantity."""
    r = dict(RESOLUTIONS); r.update(res or {})
    n = len(d["e_in"])

    def gauss(mu, rel):
        return mu * (1.0 + rng.normal(0.0, 1.0, size=n) * rel)

    ein = np.maximum(gauss(d["e_in"], r["muon_scale"] * sigma_p_over_p(d["e_in"])), 1.0)
    pout = np.maximum(gauss(d["p_out"], r["muon_scale"] * sigma_p_over_p(d["p_out"])), 0.5)
    thmu = np.maximum(d["theta_mu"] + rng.normal(0.0, r["sigma_theta"], size=n), 1e-6)

    h = hadronic_truth(d)
    ehad = np.maximum(gauss(h["e"], r["sigma_had"]), 1e-3)
    # smear the hadronic DIRECTION independently of its magnitude
    th_h = np.arctan2(h["pt"], np.maximum(h["pz"], 1e-9))
    th_h = np.clip(th_h + rng.normal(0.0, r["sigma_theta_had"], size=n), 1e-6, np.pi - 1e-6)
    pmag = np.hypot(h["pt"], h["pz"]) * ehad / np.maximum(h["e"], 1e-9)
    pt = pmag * np.sin(th_h)
    pz = pmag * np.cos(th_h)
    return dict(ein=ein, pout=pout, thmu=thmu, ehad=ehad, px=pt, py=0.0 * pt, pz=pz)


def _x_from(q2, nu):
    nu = np.where(nu > 1e-6, nu, np.nan)
    return q2 / (2.0 * M_N * nu)


def reconstruct(d, rng, res=None):
    """Return {method: x_reco} for the four methods."""
    s = smear(d, rng, res)
    ein, pout, th = s["ein"], s["pout"], s["thmu"]
    one_m_cos = 1.0 - np.cos(th)

    out = {}

    # 1. lepton-only
    q2_l = 2.0 * ein * pout * one_m_cos
    out["lepton-only"] = _x_from(q2_l, ein - pout)

    # 2. Jacquet-Blondel: nu straight from the calorimeter.
    # The measured hadronic energy is nu + M (it includes the struck nucleon's
    # rest energy), so the target mass must be subtracted to get the energy
    # transfer.  Omitting it biases x_JB low by M/nu -- 9% on the x>=0.2 sample,
    # where nu is only ~9 GeV.
    nu_jb = np.maximum(s["ehad"] - M_N, 1e-6)
    pt_h = np.hypot(s["px"], s["py"])
    y_jb = np.clip(nu_jb / ein, 1e-6, 0.999)
    q2_jb = pt_h**2 / (1.0 - y_jb)
    out["Jacquet-Blondel"] = _x_from(q2_jb, nu_jb)

    # 3. Sigma (fixed-target): rebuild E_in from the final state via E + p_z
    ein_sig = 0.5 * (s["ehad"] + s["pz"] + pout * (1.0 + np.cos(th)) - M_N)
    ein_sig = np.maximum(ein_sig, pout + 1e-3)
    q2_s = 2.0 * ein_sig * pout * one_m_cos
    out["Sigma"] = _x_from(q2_s, ein_sig - pout)

    # 4. double-angle
    th_h = np.arctan2(np.hypot(s["px"], s["py"]), np.maximum(s["pz"], 1e-9))
    denom = np.sin(th + th_h)
    e_da = np.where(np.abs(denom) > 1e-6, ein * np.sin(th_h) / denom, np.nan)
    e_da = np.clip(e_da, 0.5, None)
    q2_da = 2.0 * ein * e_da * one_m_cos
    out["double-angle"] = _x_from(q2_da, ein - e_da)

    return out


def large_x_survival(d_fit, d_pert, rng, x_cuts=(0.2, 0.4), res=None, n_toys=1):
    """
    The headline number: the fitted/perturbative charm ratio at large x, at
    TRUTH level and after reconstruction by each method.

    A method preserves the intrinsic-charm signal if its reconstructed ratio
    stays close to the truth-level one.
    """
    import fasercal_chain as F
    out = {}
    for cut in x_cuts:
        row = {}
        tf, tp = truth(d_fit), truth(d_pert)
        wf = F.event_weight(d_fit) * (d_fit["n_charm"] >= 1)
        wp = F.event_weight(d_pert) * (d_pert["n_charm"] >= 1)
        nf = wf[tf["x"] >= cut].sum(); npr = wp[tp["x"] >= cut].sum()
        row["truth"] = nf / npr if npr > 0 else np.nan
        acc = {}
        for _ in range(n_toys):
            rf = reconstruct(d_fit, rng, res)
            rp = reconstruct(d_pert, rng, res)
            for meth in rf:
                a = wf[np.nan_to_num(rf[meth], nan=-1) >= cut].sum()
                b = wp[np.nan_to_num(rp[meth], nan=-1) >= cut].sum()
                acc.setdefault(meth, []).append(a / b if b > 0 else np.nan)
        for meth, v in acc.items():
            row[meth] = float(np.nanmean(v))
        out[cut] = row
    return out
