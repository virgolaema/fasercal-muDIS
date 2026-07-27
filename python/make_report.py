#!/usr/bin/env python3
"""
FASERcal charm/anticharm-tagging toy: turn the showered charm-muon cache
(shower_charm.py) into a PDF report.

Answers, with the real Pythia distributions rather than back-of-envelope:
  * how energetic (hence how far-travelling) the semileptonic decay muon is,
  * how far the charm hadron flies before decaying (the FASERnu vertex handle,
    here compared to scintillator cell size),
  * what fraction of charm yields a sign-able muon, and
  * the resulting c-vs-c~ tagging efficiency, purity, and asymmetry dilution
    as a function of the downstream-spectrometer assumptions.

Reads only the .npz cache + detector.py; no re-showering.  Fast to iterate.
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import config
import detector as D

# D-meson PDG ids for the flight-length breakdown
DPM, D0, DS, LC = 411, 421, 431, 4122


def load(npz):
    f = np.load(npz)
    return {k: f[k] for k in f.files}


def summarise_metrics(d, params, rng):
    """Weighted tagging metrics under a given detector-parameter dict."""
    w = d["w"]
    has_mu = d["has_mu"] > 0.5
    n_charm = w.sum()
    n_mu = w[has_mu].sum()

    acc = has_mu & D.accepts(d["mu_p"], d["mu_theta"], **params)
    n_acc = w[acc].sum()

    true_q = d["sign"][acc].astype(int)
    reco_q = D.tag_charge(true_q, d["mu_p"][acc], rng, **params)
    wa = w[acc]
    correct = wa[reco_q == true_q].sum()
    purity = correct / n_acc if n_acc > 0 else 0.0
    eta_eff = 1.0 - purity

    # true and measured charm asymmetry A = (Nc - Nc~)/(Nc + Nc~)
    nc_t = wa[true_q > 0].sum();  nca_t = wa[true_q < 0].sum()
    A_true = (nc_t - nca_t) / (nc_t + nca_t) if (nc_t + nca_t) > 0 else 0.0
    nc_r = wa[reco_q > 0].sum();  nca_r = wa[reco_q < 0].sum()
    A_meas = (nc_r - nca_r) / (nc_r + nca_r) if (nc_r + nca_r) > 0 else 0.0

    return dict(n_charm=n_charm, n_mu=n_mu, n_acc=n_acc,
                br_mu=n_mu / n_charm if n_charm else 0.0,
                acc_frac=n_acc / n_mu if n_mu else 0.0,
                tag_eff=n_acc / n_charm if n_charm else 0.0,
                purity=purity, eta_eff=eta_eff,
                A_true=A_true, A_meas=A_meas)


def page_muon_spectrum(pdf, d):
    """Decay-muon lab momentum + the range it implies in scintillator."""
    has_mu = d["has_mu"] > 0.5
    p = d["mu_p"][has_mu]
    w = d["w"][has_mu]
    p = p[p > 0]
    fig, ax = plt.subplots(figsize=(11, 6))
    bins = np.logspace(np.log10(0.3), np.log10(max(p.max(), 10)), 60)
    ax.hist(p, bins=bins, weights=w[:len(p)], color="#377eb8", alpha=0.85,
            edgecolor="white")
    ax.set_xscale("log")
    ax.set_xlabel("Semileptonic decay-muon lab momentum  $p_\\mu$ [GeV]", fontsize=12)
    ax.set_ylabel("Charm hadrons / bin (tungsten-weighted)", fontsize=12)
    med = np.median(p); q1, q3 = np.percentile(p, [25, 75])
    for v, ls, lab in [(med, "-", f"median {med:.0f} GeV"),
                       (q1, "--", f"25% {q1:.0f}"), (q3, "--", f"75% {q3:.0f}")]:
        ax.axvline(v, color="black", ls=ls, lw=1.3, alpha=0.7)
    ax.set_title("Charm decay muon is 'soft' only within the jet: lab spectrum is tens of GeV\n"
                 f"median $p_\\mu$ = {med:.0f} GeV  (IQR {q1:.0f}-{q3:.0f} GeV)  "
                 f"-> ionisation range {D.muon_range_m(med):.0f} m in plastic scintillator",
                 fontsize=11)
    # secondary axis: range in metres (range = p/dEdx)
    secax = ax.secondary_xaxis("top", functions=(D.muon_range_m,
                                                  lambda r: r * 100 * D.DEDX_MIP_MEV_PER_CM / 1e3))
    secax.set_xlabel("muon ionisation range in plastic scintillator [m]  "
                     "(FASERcal is a few m -> the muon exits)", fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_flight_length(pdf, d):
    """Charm-hadron flight length vs scintillator cell size (the vertex handle)."""
    L = d["Ldecay"] * 1e4  # mm -> um  (report in microns for the short D0, cm for D+)
    apid = np.abs(d["pid"]).astype(int)
    fig, ax = plt.subplots(figsize=(11, 6))
    bins = np.logspace(np.log10(1), np.log10(max(L.max(), 1e5)), 60)  # 1 um .. 10 cm
    for pid, col, lab in [(DPM, "#e41a1c", "$D^\\pm$ (c$\\tau$=312 $\\mu$m)"),
                          (D0,  "#377eb8", "$D^0$ (c$\\tau$=123 $\\mu$m)"),
                          (DS,  "#4daf4a", "$D_s^\\pm$"),
                          (LC,  "#984ea3", "$\\Lambda_c^+$")]:
        m = apid == pid
        if m.sum() == 0:
            continue
        ax.hist(L[m], bins=bins, weights=d["w"][m], histtype="step", lw=2,
                color=col, label=f"{lab}, N={m.sum()}")
    ax.set_xscale("log")
    ax.axvspan(5e3, 3e4, color="gray", alpha=0.18,
               label="scintillator cell ~0.5-3 cm")
    ax.set_xlabel("Charm-hadron flight length before decay [$\\mu$m]", fontsize=12)
    ax.set_ylabel("Charm hadrons / bin (tungsten-weighted)", fontsize=12)
    med_all = np.median(L[L > 0])
    ax.set_title("Decay vertex vs scintillator granularity: charm flies ~mm, below cell size\n"
                 f"median flight length = {med_all:.0f} $\\mu$m "
                 f"({med_all/1e4:.2f} cm) -> emulsion resolves it, FASERcal essentially cannot",
                 fontsize=11)
    ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_cutflow_table(pdf, d, params, rng):
    m = summarise_metrics(d, params, rng)
    fig, ax = plt.subplots(figsize=(11, 6.5)); ax.axis("off")
    rows = [
        ("Charm hadrons (weighted)",            f"{m['n_charm']:.3g}", "100%"),
        ("-> with a semileptonic muon",          f"{m['n_mu']:.3g}",   f"{100*m['br_mu']:.1f}%"),
        (f"-> muon forward & $p_\\mu$>{params['p_min_gev']:.0f} GeV "
         f"(<{params['theta_max_mrad']:.0f} mrad)", f"{m['n_acc']:.3g}",
                                                  f"{100*m['acc_frac']:.1f}% of muons"),
        ("Overall sign-able fraction / charm",   "",                  f"{100*m['tag_eff']:.2f}%"),
        ("", "", ""),
        ("Charge-tag purity (accepted muons)",   "",                  f"{100*m['purity']:.1f}%"),
        ("Mean charge confusion $\\eta$",        "",                  f"{100*m['eta_eff']:.1f}%"),
        ("True charm asymmetry $A$",             "",                  f"{m['A_true']:+.3f}"),
        ("Measured asymmetry $A_{meas}=(1-2\\eta)A$", "",             f"{m['A_meas']:+.3f}"),
    ]
    tab = ax.table(cellText=[[r[0], r[1], r[2]] for r in rows],
                   colLabels=["FASERcal c/c~ tagging cutflow", "yield", "fraction"],
                   colWidths=[0.60, 0.20, 0.20], loc="center", cellLoc="left")
    tab.auto_set_font_size(False); tab.set_fontsize(10.5); tab.scale(1, 1.7)
    for j in range(3):
        tab[0, j].set_facecolor("#333333")
        tab[0, j].set_text_props(color="white", fontweight="bold")
    ax.set_title("FASERcal charm/anticharm tagging — toy detector cutflow\n"
                 f"(spectrometer: $\\theta$<{params['theta_max_mrad']:.0f} mrad, "
                 f"$p_\\mu$>{params['p_min_gev']:.0f} GeV, "
                 f"charge conf. $p_{{1/2}}$={params['p_half_gev']:.0f} GeV)",
                 fontsize=12, pad=18)
    fig.text(0.5, 0.06,
             "Key number: the sign-able fraction per charm hadron sets the tagging efficiency; "
             "the asymmetry dilution (1-2$\\eta$) sets how much of the intrinsic-charm c/c~ "
             "signal survives.", ha="center", fontsize=9, style="italic", wrap=True)
    pdf.savefig(fig); plt.close(fig)


def page_phalf_scan(pdf, d, base_params, rng):
    """How tagging purity & asymmetry survival depend on the spectrometer reach."""
    p_halfs = np.array([200, 400, 700, 1000, 1500, 2500, 5000.0])
    purity, dilution, tageff = [], [], []
    for ph in p_halfs:
        pr = dict(base_params); pr["p_half_gev"] = ph
        m = summarise_metrics(d, pr, np.random.default_rng(7))
        purity.append(100 * m["purity"])
        dilution.append(100 * (1 - 2 * m["eta_eff"]))   # asymmetry survival
        tageff.append(100 * m["tag_eff"])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ax1.plot(p_halfs, purity, "o-", color="#377eb8", label="charge-tag purity")
    ax1.plot(p_halfs, dilution, "s--", color="#e41a1c",
             label="asymmetry survival (1-2$\\eta$)")
    ax1.set_xscale("log")
    ax1.set_xlabel("spectrometer charge-confusion $p_{1/2}$ [GeV]", fontsize=11)
    ax1.set_ylabel("[%]", fontsize=11); ax1.set_ylim(0, 105)
    ax1.legend(fontsize=9); ax1.grid(alpha=0.3, which="both")
    ax1.set_title("Sign quality vs spectrometer reach", fontsize=11)

    ax2.axhline(np.mean(tageff), color="#4daf4a", lw=2)
    ax2.plot(p_halfs, tageff, "d-", color="#4daf4a")
    ax2.set_xscale("log")
    ax2.set_xlabel("spectrometer charge-confusion $p_{1/2}$ [GeV]", fontsize=11)
    ax2.set_ylabel("sign-able fraction / charm [%]", fontsize=11)
    ax2.grid(alpha=0.3, which="both")
    ax2.set_title("Tagging efficiency (set by BR x acceptance,\n~independent of charge reach)",
                  fontsize=11)
    fig.suptitle("What the spectrometer must deliver: purity climbs steeply once "
                 "$p_{1/2}$ exceeds the muon momenta", fontsize=12)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", default=config.get("charm_npz"))
    ap.add_argument("--output", default=config.get("report_pdf"))
    args = ap.parse_args()

    d = load(args.npz)
    params = dict(D.DEFAULTS)
    rng = np.random.default_rng(42)

    import os
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with PdfPages(args.output) as pdf:
        page_muon_spectrum(pdf, d)
        page_flight_length(pdf, d)
        page_cutflow_table(pdf, d, params, rng)
        page_phalf_scan(pdf, d, params, rng)

    m = summarise_metrics(d, params, np.random.default_rng(42))
    print(f"Done -> {args.output}")
    print(f"  charm hadrons: {m['n_charm']:.3g} (weighted)")
    print(f"  semileptonic-mu BR: {100*m['br_mu']:.1f}%")
    print(f"  sign-able fraction/charm: {100*m['tag_eff']:.2f}%")
    print(f"  charge-tag purity: {100*m['purity']:.1f}%  (eta={100*m['eta_eff']:.1f}%)")
    print(f"  asymmetry: A_true={m['A_true']:+.3f}  A_meas={m['A_meas']:+.3f}")


if __name__ == "__main__":
    main()
