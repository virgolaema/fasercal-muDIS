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
    # Subtract the energy carried off by neutrinos from charm semileptonic
    # decays.  This is an IRREDUCIBLE, ONE-SIDED loss and it is specific to the
    # signal events -- exactly the ones we select.  From the 25-seed sample:
    # nonzero in 35% of charm events, median 13.2 GeV = 8.4% of nu when nonzero,
    # and >10% of nu for 21% of the x>=0.2 charm sample.
    # Caveat: the 9% hadronic resolution quoted by the Bern talk is derived from
    # full simulation of neutrino events, which also contain charm decays, so
    # this may partially double-count.  Switchable for exactly that reason.
    subtract_nu=True,
)


def sigma_p_over_p(p):
    """
    CDR muon-momentum resolution at momentum p [GeV].

    The CDR quotes only four points (20, 100, 200, 1000 GeV); our median
    scattered muon sits at 469 GeV, so the value there is INTERPOLATED and the
    interpolation rule matters.

    For a magnetic spectrometer the resolution has the form
        sigma(1/p)/(1/p) = sqrt( a^2 + (b p)^2 )
    a constant multiple-scattering floor in quadrature with a measurement term
    linear in p.  A fit to the four quoted points gives a = 0.225, b = 5.9e-4,
    confirming the form.  So sigma^2 is linear in p^2, and interpolating
    sigma^2 against p^2 both passes exactly through every quoted point AND has
    the physically correct shape in between.

    An earlier version interpolated linearly in log p, which has no physical
    basis and gave 47% at 469 GeV against 38% here -- i.e. it was ~25% too
    PESSIMISTIC about the spectrometer.
    """
    p = np.asarray(p, dtype=float)
    p2 = np.clip(p, 1.0, None) ** 2
    return np.sqrt(np.interp(p2, _P_NODES ** 2, _SIG_NODES ** 2))


def truth(d):
    """Truth-level nu, x, y from the cached DIS kinematics."""
    nu = d["e_in"] - d["p_out"]
    x = d["q2"] / (2.0 * M_N * np.maximum(nu, 1e-9))
    return dict(nu=nu, x=x, y=nu / np.maximum(d["e_in"], 1e-9), q2=d["q2"])


def hadronic_truth(d, source="auto"):
    """
    The TRUE hadronic four-vector.

    Two definitions are available and, for the BINNED production, they agree:

      "particles"   -- the direct sum over the Pythia final state, excluding the
                       scattered muon and all neutrinos.  This is the honest
                       one: it carries the real fragmentation fluctuations and
                       the real charged/neutral composition.
      "conservation" -- from the muon kinematics alone, using
                            X = q + p_target ,  p_target = (M, 0, 0, 0)
                        so  E = nu + M,  pT = p_out sin(theta),
                            pz = |p_in| - p_out cos(theta).

    HISTORY.  The original samples took the flux as a beam PDF of a fictitious
    7 TeV muon beam.  That beam has a remnant which hadronises into the event,
    so the particle sum was contaminated: energy conservation failed in 80% of
    events and E_had/nu reached 2.5 at large x, where nu is small.  No
    particle-level cut removed it (beam-remnant status codes, forward cones and
    mother-chain ancestry were all tried and all failed -- see docs/XRECO.md),
    which forced "conservation" as a workaround.

    The workaround had a real cost: it made E_had a deterministic function of
    the muon kinematics, so the calorimetric methods (JB, Sigma) were being fed
    an input partly built from the very quantity they are supposed to replace.
    The binned production removes the remnant at source (fixed_lepton_beam 1 =>
    x_lepton = 1 => no lepton-side remnant), the particle sum closes, and the
    circularity is gone.

    "auto" uses the particle sum when the cache provides it (binned production,
    tagged by the `ibin` column) and falls back to conservation otherwise.

    NOT MODELLED, and flagged as such:
      * energy carried off by neutrinos from charm semileptonic decays, which
        biases the measurable E_had low for exactly the signal events (the
        neutrinos are excluded from the sum here, so this is the truth the
        detector could at best see, not the full hadronic system);
      * the detector response, folded into the single sigma_had parameter.
    """
    if source == "auto":
        source = "particles" if ("ibin" in d and "px_had" in d) else "conservation"
    if source == "particles":
        return dict(
            e=d["e_had"],
            pt=np.hypot(d["px_had"], d["py_had"]),
            pz=d["pz_had"],
        )
    nu = d["e_in"] - d["p_out"]
    return dict(
        e=nu + M_N,
        pt=d["p_out"] * np.sin(d["theta_mu"]),
        pz=d["e_in"] - d["p_out"] * np.cos(d["theta_mu"]),
    )


def smear(d, rng, res=None, had_source="auto"):
    """Apply the detector response to every measured quantity.

    `had_source` selects the definition of the true hadronic four-vector; see
    hadronic_truth().  It is exposed so that the two can be compared on the same
    sample, which is the check that retires the momentum-conservation workaround.
    """
    r = dict(RESOLUTIONS); r.update(res or {})
    n = len(d["e_in"])

    def gauss(mu, rel):
        return mu * (1.0 + rng.normal(0.0, 1.0, size=n) * rel)

    def smear_muon(p, rel):
        """
        Smear a muon momentum the way the spectrometer actually measures it.

        The CDR quotes sigma(Delta(1/p)/(1/p)), i.e. the resolution is Gaussian
        in the CURVATURE 1/p, not in p.  At these large resolutions the two are
        not interchangeable: smearing p multiplicatively produces a symmetric
        distribution, whereas smearing 1/p produces the correct long high-p tail
        and can even flip the sign of the curvature (an unmeasurably straight
        track), which is physically what happens.
        """
        inv = (1.0 / np.maximum(p, 1e-6)) * (1.0 + rng.normal(0.0, 1.0, size=n) * rel)
        inv = np.where(np.abs(inv) < 1e-9, 1e-9, inv)     # guard the pole
        out = 1.0 / np.abs(inv)
        return np.clip(out, 0.5, 1.0e5)                   # cap runaway tails

    ein = smear_muon(d["e_in"], r["muon_scale"] * sigma_p_over_p(d["e_in"]))
    pout = smear_muon(d["p_out"], r["muon_scale"] * sigma_p_over_p(d["p_out"]))
    thmu = np.maximum(d["theta_mu"] + rng.normal(0.0, r["sigma_theta"], size=n), 1e-6)

    h = hadronic_truth(d, had_source)
    e_vis = h["e"]
    if r["subtract_nu"] and "e_nu_charm" in d:
        e_vis = np.maximum(e_vis - d["e_nu_charm"], 1e-3)
    ehad = np.maximum(gauss(e_vis, r["sigma_had"]), 1e-3)
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


def reconstruct(d, rng, res=None, ein_known=True, had_source="auto"):
    """
    Return {method: x_reco} for the four methods.

    ein_known=True   the incoming muon energy is measured (e.g. by an upstream
                     spectrometer).  All four methods are available.
    ein_known=False  it is NOT measured -- which is the situation if only the
                     DOWNSTREAM spectrometer exists, since that sees only the
                     outgoing muon.  Then:
                       * lepton-only  : impossible (nu = E_in - E' needs it)
                       * double-angle : impossible (x_DA is proportional to E_in)
                       * Jacquet-Blondel : usable, E_in enters only the mild
                         (1-y) factor; the Sigma-reconstructed E_in is used
                       * Sigma        : unaffected -- it reconstructs E_in from
                         the final state, which is exactly what it is for
                     The two impossible methods return NaN rather than a number,
                     so they cannot silently contribute to a comparison.
    """
    s = smear(d, rng, res, had_source)
    ein, pout, th = s["ein"], s["pout"], s["thmu"]
    one_m_cos = 1.0 - np.cos(th)

    # E_in as reconstructed from the final state (the Sigma estimator); this is
    # what the calorimetric methods fall back on when E_in is not measured.
    ein_sig = 0.5 * (s["ehad"] + s["pz"] + pout * (1.0 + np.cos(th)) - M_N)
    ein_sig = np.maximum(ein_sig, pout + 1e-3)
    ein_for_y = ein if ein_known else ein_sig

    out = {}
    nan = np.full(len(ein), np.nan)

    # 1. lepton-only
    if ein_known:
        q2_l = 2.0 * ein * pout * one_m_cos
        out["lepton-only"] = _x_from(q2_l, ein - pout)
    else:
        out["lepton-only"] = nan

    # 2. Jacquet-Blondel: nu straight from the calorimeter.
    # The measured hadronic energy is nu + M (it includes the struck nucleon's
    # rest energy), so the target mass must be subtracted to get the energy
    # transfer.  Omitting it biases x_JB low by M/nu -- 9% on the x>=0.2 sample,
    # where nu is only ~9 GeV.
    nu_jb = np.maximum(s["ehad"] - M_N, 1e-6)
    pt_h = np.hypot(s["px"], s["py"])
    y_jb = np.clip(nu_jb / ein_for_y, 1e-6, 0.999)
    q2_jb = pt_h**2 / (1.0 - y_jb)
    out["Jacquet-Blondel"] = _x_from(q2_jb, nu_jb)

    # 3. Sigma (fixed-target): E_in rebuilt above from E + p_z
    q2_s = 2.0 * ein_sig * pout * one_m_cos
    out["Sigma"] = _x_from(q2_s, ein_sig - pout)

    # 4b. DOUBLE-ANGLE fed with the SIGMA-RECONSTRUCTED E_in.
    # x_DA is proportional to E_in, so the method is only as good as that input.
    # The Sigma estimator rebuilds E_in from the final state with the ~9%
    # calorimeter rather than the ~38% spectrometer, so this variant is BETTER
    # than double-angle with a measured E_in -- and it needs no incoming-muon
    # measurement at all, so it is available in the configuration where only a
    # downstream spectrometer exists.
    th_h_s = np.arctan2(np.hypot(s["px"], s["py"]), np.maximum(s["pz"], 1e-9))
    den_s = np.sin(th + th_h_s)
    e_da_s = np.where(np.abs(den_s) > 1e-6, ein_sig * np.sin(th_h_s) / den_s, np.nan)
    e_da_s = np.clip(e_da_s, 0.5, None)
    out["double-angle (Sigma E_in)"] = _x_from(2.0 * ein_sig * e_da_s * one_m_cos,
                                               ein_sig - e_da_s)

    # 4. double-angle
    if ein_known:
        th_h = np.arctan2(np.hypot(s["px"], s["py"]), np.maximum(s["pz"], 1e-9))
        denom = np.sin(th + th_h)
        e_da = np.where(np.abs(denom) > 1e-6, ein * np.sin(th_h) / denom, np.nan)
        e_da = np.clip(e_da, 0.5, None)
        q2_da = 2.0 * ein * e_da * one_m_cos
        out["double-angle"] = _x_from(q2_da, ein - e_da)
    else:
        out["double-angle"] = nan

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
