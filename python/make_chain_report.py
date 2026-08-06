#!/usr/bin/env python3
"""
PDF report for the FASERcal muon-DIS charm chain (see fasercal_chain.py).

Reads only the per-event .npz caches written by shower_dis.py, so it is fast and
can be re-run freely while detector assumptions are varied.

  python3 make_chain_report.py --npz .../fasercal_dis_v1.npz \
      --npz-pc .../fasercal_dis_pc.npz --output report.pdf --png-dir docs/figures
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import config
import fasercal_chain as F

PNG_DIR = None


def _emit(fig, pdf, name, tight=True):
    if tight:
        fig.tight_layout()
    pdf.savefig(fig)
    if PNG_DIR:
        os.makedirs(PNG_DIR, exist_ok=True)
        fig.savefig(os.path.join(PNG_DIR, f"{name}.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


def load(npz):
    f = np.load(npz)
    return {k: f[k] for k in f.files}


def page_conditions(pdf, d):
    """
    Front page: every condition and input the report depends on, plus a
    to-scale schematic of the FASERCal layout.  Drawn from the numbers in
    geometry.py / fasercal_chain.py so the picture cannot drift from the code.
    """
    import geometry as G
    fig = plt.figure(figsize=(11.5, 9.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.72, 1.28], hspace=0.05,
                          top=0.95, bottom=0.03)

    # ---------------- schematic (side view, z along beam) ----------------
    ax = fig.add_subplot(gs[0])
    # subdetector: (label, length cm, half-height cm, colour)
    subs = [("3DCAL\n10 modules", 255.8, 24.0, "#c9c9a0"),
            ("ECAL", 24.0, 36.0, "#b0b0b0"),
            ("AHCAL", 92.0, 43.0, "#8f8fc0"),
            ("MuSpect\n3 x Fe 1.5 T + MDT", 175.1, 50.0, "#7a7ab5")]
    z = 0.0
    for name, dz, hy, col in subs:
        ax.add_patch(plt.Rectangle((z, -hy), dz, 2 * hy, facecolor=col,
                                   edgecolor="black", lw=1.0, alpha=0.85))
        ax.text(z + dz / 2, 0, name, ha="center", va="center", fontsize=8.5,
                fontweight="bold")
        ax.text(z + dz / 2, -hy - 5, f"{dz:.0f} cm", ha="center", va="top", fontsize=7.5)
        z += dz + 6
    # the 10 3DCAL modules
    for i in range(G.N_MODULES):
        ax.plot([i * 25.58, i * 25.58], [-24, 24], color="black", lw=0.4, alpha=0.5)
    # line of sight
    ax.axhline(0, color="#c44e52", ls="--", lw=1.4)
    ax.text(z + 4, 0, "LoS", color="#c44e52", fontsize=9, va="center", fontweight="bold")
    ax.annotate("", xy=(-14, 0), xytext=(-14, 23.6),
                arrowprops=dict(arrowstyle="<->", color="#1f3b57", lw=1.4))
    ax.text(-18, 12, "LoS shift\n23.6 cm in Y\n45.2 cm in X\n(tilt 4.5$\\degree$)",
            fontsize=8, ha="right", va="center", color="#1f3b57")
    ax.set_xlim(-70, z + 30); ax.set_ylim(-62, 62)
    ax.set_xlabel("z along the beam [cm]", fontsize=10)
    ax.set_aspect("equal"); ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_title("FASERCal layout, to scale (side view) — this study uses the 3DCAL only",
                 fontsize=11.5)

    # ---------------- conditions table ----------------
    ax2 = fig.add_subplot(gs[1]); ax2.axis("off")
    g1, g5 = G.cdr_config(0.1), G.cdr_config(0.5)
    rows = [
        ("Luminosity (Run 4 nominal)", f"{F.LUMI_RUN4:.0f} fb$^{{-1}}$", "FASERCal Bern CM talk"),
        ("Fiducial volume", "3DCAL only, AHCAL excluded", "as requested"),
        ("3DCAL geometry", f"{G.N_MODULES} modules x {G.LAYERS_PER_MODULE} layers of "
                           f"{G.CUBE_CM:.0f} cm$^3$ cubes", "talk slide 4"),
        ("3DCAL face", f"{G.FACE_CM:.0f} x {G.FACE_CM:.0f} cm", "talk slide 11"),
        ("Absorber", "1 mm or 5 mm W per MODULE (every 20 layers)", "talk slide 39"),
        ("3DCAL total mass", f"{g1['mass_full_kg']:.0f} kg (1 mm) / {g5['mass_full_kg']:.0f} kg (5 mm)",
         "CDR Table 5"),
        # The CDR's "z < 1150 mm" is not a containment cut but the requirement to
        # be inside the 3DCal, in an older centred coordinate system.  Using the
        # 3DCal alone already enforces it, so no fiducial cut is applied.
        ("Fiducial cut", "none -- 3DCal volume only", "CDR Table 8"),
        ("Target mass", f"{1e3*g1['m_tot']:.0f} kg (1 mm) / {1e3*g5['m_tot']:.0f} kg (5 mm)", "derived"),
        ("Charge misidentification", "0.7% ($<$100 GeV) $\\to$ 2-5% (TeV)", "CDR Sec. 2.6.2"),
        ("Spectrometer", "10 Fe planes 1.5 T, 11 SciFi stations 100 $\\mu$m", "CDR Table 7"),
        ("Radiation / collision length (5 mm)", f"{g5['x0_total']:.1f} $X_0$ / "
                                                f"{g5['lambda_int']:.1f} $\\lambda$", "CDR Table 5"),
        ("Position", f"LoS shift ({G.LOS_SHIFT_X_MM:.0f}, {G.LOS_SHIFT_Y_MM:.0f}) mm, "
                     f"tilt {G.TILT_DEG:.1f}$\\degree$", "talk slides 11-12"),
        ("Target composition", "polystyrene CH, $w_p$=0.538 (+ W where present)", "computed per event"),
        ("Muon flux", "FASERnu Run 3 LHAPDF grid 'var2' (25x30 cm)", "arXiv:2506.13889 Eq. 2.1"),
        ("Off-axis flux factor", f"$F_{{\\rm flux}}$ = {F.F_FLUX:.2f} at (+452,+236) mm",
         "FLUKA map, this work"),
        ("Run 4 optics enhancement", f"$F_{{\\rm optics}}$ = {F.F_OPTICS:.1f} (rate $\\times$2, harder)",
         "CDR Sec. 1.3.2"),
        ("AHCAL as extra target (optional)", "3.33 t total, 1.60 t fiducial, ~98% Fe",
         "CDR Table 3 / talk"),
        ("Detector side", "left of LoS looking downstream = $+x$", "user; right-handed frame"),
        ("Spectrometer acceptance", f"{100*F.DEFAULTS['acc_spectrometer']:.0f}% "
                                    "($\\geq$2 MDT stations)", "talk slide 33"),
        ("Generator", "POWHEG + Pythia8, NNPDF4.0 fitted vs perturbative charm", "10 seeds, 764k events"),
        ("Charge tag", "semileptonic $\\mu$ charge = charm sign", "verified in code, 100%"),
    ]
    tab = ax2.table(cellText=rows, colLabels=["condition", "value", "source"],
                    colWidths=[0.32, 0.44, 0.24], loc="center", cellLoc="left")
    tab.auto_set_font_size(False); tab.set_fontsize(8.6); tab.scale(1, 1.42)
    for j in range(3):
        tab[0, j].set_facecolor("#1f3b57")
        tab[0, j].set_text_props(color="white", fontweight="bold")
    ax2.set_title("All conditions and inputs used in this report", fontsize=12, pad=26)
    _emit(fig, pdf, "09_conditions", tight=False)


def page_cutflow(pdf, d, m, m_ah=None):
    """The requested chain, with and without the AHCAL as an extra target."""
    import geometry as G
    fig, ax = plt.subplots(figsize=(12.5, 7)); ax.axis("off")
    a = m_ah or {k: float("nan") for k in m}
    def two(k, fmt="{:,.0f}"):
        return [fmt.format(m[k]), fmt.format(a[k])]
    rows = [
        ["Muon DIS interactions in fiducial volume"] + two("n_dis") + ["100%"],
        ["  ... containing charm ($\\geq$1 charm hadron)"] + two("n_charm")
            + [f"{100*m['f_charm']:.2f}%"],
        ["  ... charm decaying semileptonically to $\\mu$"] + two("n_semi")
            + [f"{100*m['f_semi']:.1f}% of charm"],
        [f"  ... $\\mu$ punches out ($p>{F.DEFAULTS['p_punch_gev']:.0f}$ GeV)"] + two("n_punch")
            + [f"{100*m['f_punch']:.1f}%"],
        [f"  ... reaches spectrometer ({100*F.DEFAULTS['acc_spectrometer']:.0f}%)"] + two("n_acc")
            + [f"{100*F.DEFAULTS['acc_spectrometer']:.0f}%"],
        [f"  ... identified + linked ($\\varepsilon$={F.DEFAULTS['eff_link']:.2f})"] + two("n_tag")
            + [f"{100*m['f_overall']:.3f}% of DIS"],
        ["", "", "", ""],
        ["Charge-tag purity (c vs $\\bar{c}$)", f"{100*m['purity_q']:.1f}%",
         f"{100*a['purity_q']:.1f}%", "CDR Sec. 2.6.2"],
        ["Statistical reach on $A_c$", f"$\\pm${m['sigma_A']:.3f}",
         f"$\\pm${a['sigma_A']:.3f}", ""],
    ]
    tab = ax.table(cellText=rows,
                   colLabels=[f"FASERcal chain — Run 4, {F.LUMI_RUN4:.0f} fb$^{{-1}}$",
                              f"3DCAL only\n({1e3*G.cdr_config(0.1)['m_tot']:.0f} kg)",
                              f"3DCAL + AHCAL\n({1e3*(G.cdr_config(0.1)['m_tot']+G.ahcal_config()['m_tot']):.0f} kg)",
                              "fraction"],
                   colWidths=[0.44, 0.19, 0.19, 0.18], loc="center", cellLoc="left")
    tab.auto_set_font_size(False); tab.set_fontsize(10); tab.scale(1, 1.85)
    for j in range(4):
        tab[0, j].set_text_props(fontweight="bold")
    for j in range(4):
        tab[0, j].set_facecolor("#1f3b57")
        tab[0, j].set_text_props(color="white", fontweight="bold")
    fig.text(0.5, 0.09,
             "3DCAL geometry and mass are now the DESIGNED values (Bern CM talk, 15 Jul 2026): "
             "10 modules x 20 layers of 1 cm cubes, 48x48 cm face, W per module.\n"
             "The detector sits LEFT of the LoS looking downstream, i.e. at (+452, +236) mm, where the "
             "FLUKA map gives F_flux = 0.96 — the quieter lobe.\n"
             "(The mirror position would have given 2.02, with a harder, mu+-enriched beam.) "
             "Yields scale linearly with F_flux and with the efficiency; the chain fractions are ratios.",
             ha="center", fontsize=9, style="italic")
    ax.set_title("Muon-DIS charm chain in FASERcal (AHCAL excluded)", fontsize=13, pad=16)
    _emit(fig, pdf, "10_chain_cutflow", tight=False)


def page_funnel(pdf, m):
    """Log-scale waterfall of the chain."""
    labels = ["DIS", "charm", r"semilep. $\mu$", "punch-out", "spectrometer", "tagged"]
    vals = [m["n_dis"], m["n_charm"], m["n_semi"], m["n_punch"], m["n_acc"], m["n_tag"]]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    cols = plt.cm.viridis(np.linspace(0.15, 0.85, len(vals)))
    ax.bar(labels, vals, color=cols, edgecolor="black", lw=0.6)
    for i, v in enumerate(vals):
        ax.text(i, v * 1.25, f"{v:,.0f}", ha="center", fontsize=10, fontweight="bold")
        if i:
            ax.text(i, v * 0.35, f"x{v/vals[i-1]:.3f}", ha="center", fontsize=8.5,
                    color="white", fontweight="bold")
    ax.set_yscale("log"); ax.set_ylabel("Expected events (Run 4)", fontsize=12)
    ax.set_ylim(max(min(vals) * 0.2, 0.1), max(vals) * 6)
    ax.grid(axis="y", alpha=0.3, which="both")
    ax.set_title(f"Chain attrition — {F.LUMI_RUN4:.0f} fb$^{{-1}}$, {F.M_FID_T:.1f} t fiducial. "
                 f"Overall {100*m['f_overall']:.3f}% of DIS events give a signed, linked charm muon",
                 fontsize=11.5)
    _emit(fig, pdf, "11_funnel")


def page_muon_spectrum(pdf, d):
    """Decay-muon momentum vs the thresholds that define the chain."""
    has = d["n_semilep_mu"] >= 1
    p = d["mu2_p"][has]; p = p[np.isfinite(p) & (p > 0)]
    w = F.event_weight(d)[has][: len(p)]
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.hist(p, bins=np.logspace(np.log10(0.2), np.log10(max(p.max(), 10)), 55),
            weights=w, color="#4c72b0", alpha=0.85, edgecolor="white")
    ax.axvline(F.DEFAULTS["p_punch_gev"], color="#c44e52", ls="--", lw=2,
               label=f"punch-through $p>{F.DEFAULTS['p_punch_gev']:.0f}$ GeV")
    ax.set_xscale("log")
    med = np.median(p)
    ax.axvline(med, color="black", ls=":", lw=1.5, label=f"median {med:.1f} GeV")
    ax.set_xlabel(r"semileptonic charm-decay muon $p_\mu$ [GeV]", fontsize=12)
    ax.set_ylabel("Expected events / bin (Run 4)", fontsize=12)
    ax.legend(fontsize=10); ax.grid(alpha=0.3, which="both")
    ax.set_title("The muon that carries the charm charge: this spectrum sets both the\n"
                 "punch-through efficiency and the charge-confusion dilution", fontsize=11.5)
    _emit(fig, pdf, "12_decay_muon_p")


def page_angular_acceptance(pdf, d):
    """
    What the assumed 40% spectrometer acceptance implies geometrically.

    The charm-decay muon is emitted at a much larger angle than the scattered
    DIS muon, so acceptance is dominated by the angular aperture, not by
    momentum.  Mapping the assumed 40% onto an equivalent cone gives the
    collaboration a directly checkable number.
    """
    has = d["n_semilep_mu"] >= 1
    th = d["mu2_theta"][has]; p = d["mu2_p"][has]
    ok = np.isfinite(th) & np.isfinite(p) & (p > F.DEFAULTS["p_punch_gev"])
    th = th[ok] * 1e3                                   # mrad
    th_mu = d["theta_mu"][has][ok] * 1e3

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2))
    bins = np.logspace(np.log10(0.5), np.log10(1000), 55)
    ax1.hist(th_mu, bins=bins, histtype="step", lw=2, color="#4c72b0",
             label=rf"scattered DIS $\mu$ (median {np.median(th_mu):.1f} mrad)")
    ax1.hist(th, bins=bins, histtype="step", lw=2, color="#c44e52",
             label=rf"charm-decay $\mu$ (median {np.median(th):.0f} mrad)")
    ax1.set_xscale("log"); ax1.set_xlabel(r"$\theta$ w.r.t. beam [mrad]", fontsize=11)
    ax1.set_ylabel("events / bin", fontsize=11)
    ax1.legend(fontsize=9); ax1.grid(alpha=0.3, which="both")
    ax1.set_title("The decay muon is an order of magnitude wider\nthan the scattered muon",
                  fontsize=11)

    cones = np.logspace(np.log10(2), np.log10(500), 60)
    frac = [100 * (th < c).mean() for c in cones]
    ax2.plot(cones, frac, "-", color="#1f3b57", lw=2.2)
    a = 100 * F.DEFAULTS["acc_spectrometer"]
    th_eq = np.percentile(th, a)
    ax2.axhline(a, color="#c44e52", ls="--", lw=1.6, label=f"assumed acceptance {a:.0f}%")
    ax2.axvline(th_eq, color="#c44e52", ls=":", lw=1.6,
                label=rf"$\Rightarrow$ equivalent cone $\theta<{th_eq:.0f}$ mrad")
    ax2.set_xscale("log"); ax2.set_xlabel(r"angular half-aperture [mrad]", fontsize=11)
    ax2.set_ylabel("decay muons accepted [%]", fontsize=11)
    ax2.set_ylim(0, 102); ax2.legend(fontsize=9); ax2.grid(alpha=0.3, which="both")
    ax2.set_title("Assumed 40% acceptance $\\equiv$ a "
                  f"$\\theta<{th_eq:.0f}$ mrad cone\n(checkable against the real aperture)",
                  fontsize=11)
    fig.suptitle("Acceptance is an ANGULAR problem, not a momentum one", fontsize=12.5)
    _emit(fig, pdf, "15_angular_acceptance")


def page_mass_and_params(pdf, d):
    """The two dominant assumptions: fiducial mass, and the toy response."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2))

    masses = np.array([0.25, 0.5, 1.0, 2.0, 5.0, 10.0])
    tags = [F.chain(d, m_fid_t=mm)["n_tag"] for mm in masses]
    ax1.plot(masses, tags, "o-", color="#1f3b57", lw=2)
    ax1.axvline(F.M_FID_T, color="#c44e52", ls="--", label=f"nominal {F.M_FID_T:.0f} t")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("FASERcal fiducial mass [t]", fontsize=11)
    ax1.set_ylabel("tagged charm muons (Run 4)", fontsize=11)
    ax1.legend(fontsize=9); ax1.grid(alpha=0.3, which="both")
    ax1.set_title("Yield is LINEAR in fiducial mass\n(rescale when 3DCAL geometry is fixed)",
                  fontsize=11)

    phalf = np.array([50, 100, 200, 400, 800, 1500, 3000.0])
    pur = [100 * F.chain(d, params={"p_half_gev": ph})["purity_q"] for ph in phalf]
    dil = [100 * (2 * F.chain(d, params={"p_half_gev": ph})["purity_q"] - 1) for ph in phalf]
    ax2.plot(phalf, pur, "o-", color="#4c72b0", label="charge-tag purity")
    ax2.plot(phalf, dil, "s--", color="#c44e52", label=r"asymmetry survival $(1-2\eta)$")
    ax2.set_xscale("log"); ax2.set_ylim(0, 105)
    ax2.set_xlabel(r"spectrometer charge-confusion $p_{1/2}$ [GeV]", fontsize=11)
    ax2.set_ylabel("[%]", fontsize=11)
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3, which="both")
    ax2.set_title("Charge ID is NOT the bottleneck:\nthe decay muons are soft in absolute terms",
                  fontsize=11)
    fig.suptitle("Dependence on the two dominant assumptions", fontsize=12.5)
    _emit(fig, pdf, "13_assumptions")


def page_ic(pdf, d, dpc):
    """Fitted (IC-allowed) vs perturbative charm, through the same chain."""
    mf, mp = F.chain(d), F.chain(dpc)
    fig, ax = plt.subplots(figsize=(11, 6))
    labels = ["charm events", r"semilep. $\mu$", "tagged"]
    vf = [mf["n_charm"], mf["n_semi"], mf["n_tag"]]
    vp = [mp["n_charm"], mp["n_semi"], mp["n_tag"]]
    x = np.arange(len(labels)); wd = 0.38
    ax.bar(x - wd/2, vf, wd, label="fitted charm (IC allowed)", color="#4c72b0",
           edgecolor="black", lw=0.6)
    ax.bar(x + wd/2, vp, wd, label="perturbative charm", color="#dd8452",
           edgecolor="black", lw=0.6)
    for i, (a, b) in enumerate(zip(vf, vp)):
        ax.text(i, max(a, b) * 1.15, f"ratio {a/b:.2f}" if b else "--",
                ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_yscale("log"); ax.set_ylabel("Expected events (Run 4)", fontsize=12)
    ax.legend(fontsize=10); ax.grid(axis="y", alpha=0.3, which="both")
    ax.set_title("Intrinsic-charm sensitivity: the two PDF hypotheses through the same chain\n"
                 "(inclusive ratios; the IC signal lives at large x, not in the total)",
                 fontsize=11.5)
    _emit(fig, pdf, "14_ic_comparison")


def page_tungsten_scenarios(pdf, d, cone):
    """Baseline vs 1 mm vs 5 mm tungsten absorber, through the whole chain."""
    import geometry as G
    ms = {nm: F.chain_scenario(d, t, cone=cone) for nm, t in G.SCENARIOS.items()}
    names = list(G.SCENARIOS)
    base = ms[names[0]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.6))
    stages = ["n_dis", "n_charm", "n_semi", "n_punch", "n_acc", "n_tag"]
    slab = ["DIS", "charm", r"semilep $\mu$", "punch", "accept", "tagged"]
    x = np.arange(len(stages)); wd = 0.26
    cols = ["#4c72b0", "#dd8452", "#c44e52"]
    for i, nm in enumerate(names):
        ax1.bar(x + (i - 1) * wd, [ms[nm][s] for s in stages], wd,
                label=nm.replace("_", " "), color=cols[i], edgecolor="black", lw=0.5)
    ax1.set_yscale("log"); ax1.set_xticks(x); ax1.set_xticklabels(slab, fontsize=10)
    ax1.set_ylabel("Expected events (Run 4)", fontsize=11)
    ax1.legend(fontsize=9); ax1.grid(axis="y", alpha=0.3, which="both")
    ax1.set_title("Tungsten multiplies the whole chain", fontsize=11.5)

    mass = [ms[nm]["geo"]["m_tot"] for nm in names]
    tag = [ms[nm]["n_tag"] for nm in names]
    sig = [ms[nm]["sigma_A"] for nm in names]
    ax2.plot(mass, tag, "o-", color="#1f3b57", lw=2, ms=9)
    for nm, mm, tt in zip(names, mass, tag):
        ax2.annotate(nm.replace("_", " "), (mm, tt), textcoords="offset points",
                     xytext=(8, -12), fontsize=9)
    ax2.set_xlabel("total fiducial mass [t]", fontsize=11)
    ax2.set_ylabel("tagged charm muons", fontsize=11)
    ax2.grid(alpha=0.3)
    axb = ax2.twinx()
    axb.plot(mass, sig, "s--", color="#c44e52", lw=1.8, ms=7)
    axb.set_ylabel(r"$\sigma(A_c)$", color="#c44e52", fontsize=11)
    axb.tick_params(axis="y", colors="#c44e52")
    ax2.set_title(r"Yield tracks mass; $\sigma(A_c)$ improves as $1/\sqrt{N}$",
                  fontsize=11.5)
    fig.suptitle("Tungsten absorber options as designed — 10 modules x 20 layers, "
                 "W per module", fontsize=12.5)
    _emit(fig, pdf, "16_tungsten_scenarios")

    # ---- companion table page ----
    fig, ax = plt.subplots(figsize=(12, 6.2)); ax.axis("off")
    rows = []
    for nm in names:
        m = ms[nm]; g = m["geo"]
        rows.append([nm.replace("_", " "),
                     f"{g['m_sc']:.2f} / {g['m_w']:.2f}",
                     f"{g['m_tot']:.2f}",
                     f"{g['x0_mean']:.1f}",
                     f"{m['theta_ms_med_mrad']:.1f}",
                     f"{100*m['f_acc']:.1f}%",
                     f"{m['n_tag']:,.0f}",
                     f"{m['sigma_A']:.3f}",
                     f"x{base['sigma_A']/m['sigma_A']:.2f}"])
    tab = ax.table(cellText=rows,
                   colLabels=["scenario", "M scint/W [t]", "M tot [t]",
                              r"$\langle X_0\rangle$", r"$\theta_{MS}$ [mrad]",
                              "accept.", "tagged", r"$\sigma(A_c)$", "gain"],
                   loc="center", cellLoc="center")
    tab.auto_set_font_size(False); tab.set_fontsize(10); tab.scale(1, 2.0)
    for j in range(9):
        tab[0, j].set_facecolor("#1f3b57")
        tab[0, j].set_text_props(color="white", fontweight="bold")
    ax.set_title("Tungsten scenarios — material budget and physics reach", fontsize=13, pad=18)
    fig.text(0.5, 0.11,
             "Multiple scattering grows steeply ($\\theta_{MS}$ 1.9 → 13.4 mrad) but costs only a few % of "
             "acceptance:\nthe decay-muon angular spread is intrinsically much wider than the aperture, so "
             "scattering moves nearly as many muons in as out.\nThe mass gain therefore wins outright.",
             ha="center", fontsize=9, style="italic")
    _emit(fig, pdf, "17_tungsten_table", tight=False)


def page_tungsten_scan(pdf, d, cone):
    """Is there an optimum absorber thickness?"""
    import geometry as G
    tw = np.array(sorted(G.CDR_TABLE5))
    res = [F.chain_scenario(d, float(t), cone=cone) for t in tw]
    tag = np.array([r["n_tag"] for r in res])
    facc = np.array([100 * r["f_acc"] for r in res])
    mtot = np.array([r["geo"]["m_tot"] for r in res])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.4))
    ax1.plot(tw * 10, tag, "-", color="#1f3b57", lw=2.4)
    for nm, t in G.SCENARIOS.items():
        m = F.chain_scenario(d, t, cone=cone)
        ax1.plot(t * 10, m["n_tag"], "o", ms=10, label=nm.replace("_", " "))
    ax1.set_xscale("log")
    ax1.set_xlabel("tungsten per module [mm]", fontsize=11)
    ax1.set_ylabel("tagged charm muons (Run 4)", fontsize=11)
    ax1.legend(fontsize=9); ax1.grid(alpha=0.3)
    ax1.set_title("Yield rises monotonically with absorber thickness", fontsize=11.5)

    ax2.plot(tw * 10, facc, "-", color="#c44e52", lw=2.2, label="acceptance [%]")
    ax2.set_xlabel("tungsten per module [mm]", fontsize=11)
    ax2.set_ylabel("spectrometer acceptance [%]", color="#c44e52", fontsize=11)
    ax2.tick_params(axis="y", colors="#c44e52")
    ax2.set_ylim(0, 50); ax2.grid(alpha=0.3)
    axb = ax2.twinx()
    axb.plot(tw * 10, mtot, "-", color="#4c72b0", lw=2.2)
    axb.set_ylabel("total mass [t]", color="#4c72b0", fontsize=11)
    axb.tick_params(axis="y", colors="#4c72b0")
    ax2.set_title("The trade-off: mass gain (blue) vs\nacceptance loss (red)", fontsize=11.5)
    fig.suptitle("Absorber thickness per MODULE (every 20 layers, as designed) — "
                 "yield rises monotonically; the limit is calorimetric", fontsize=12.5)
    _emit(fig, pdf, "18_tungsten_scan")


def page_efficiency(pdf, d, cone):
    """
    The main result as a function of the identification + linking efficiency.

    Everything downstream of the spectrometer acceptance is folded into one
    number, eff_link: the probability that a decay muon which reached the
    spectrometer is actually identified, momentum-reconstructed, charge-signed
    and linked back to the DIS vertex.  It is the least constrained ingredient
    of the chain, so rather than assume a value the result is shown across the
    full range and can be read off once the real efficiency is known.

    Scaling is simple and exact: N_tag is linear in eff, so sigma(A_c) goes as
    1/sqrt(eff).  The point of the plot is the absolute numbers, not the shape.
    """
    import geometry as G
    eff = np.linspace(0.05, 1.0, 40)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.6))
    palette = ["#4c72b0", "#dd8452", "#c44e52"]

    for i, (nm, t) in enumerate(G.SCENARIOS.items()):
        c = palette[i % len(palette)]
        lab = nm.replace("_module", "").replace("_", " ")
        for ahcal, ls, alpha in [(False, "-", 1.0), (True, "--", 0.85)]:
            base = F.chain_scenario(d, t, cone=cone, f_flux=F.F_FLUX,
                                    include_ahcal=ahcal,
                                    params={"eff_link": 1.0})
            n = base["n_tag"] * eff
            tag = lab + (" + AHCAL" if ahcal else " (3DCAL)")
            ax1.plot(eff, n, ls, lw=2.2, color=c, alpha=alpha, label=tag)
            ax2.plot(eff, 1.0 / np.sqrt(n), ls, lw=2.2, color=c, alpha=alpha,
                     label=tag)

    for ax in (ax1, ax2):
        ax.set_xlabel("identification + linking efficiency $\\varepsilon$", fontsize=11.5)
        ax.set_xlim(0.05, 1.0)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8.2, ncol=2,
                  title="solid: 3DCAL only   dashed: + AHCAL", title_fontsize=8.2)
    ax1.set_ylabel("charge-signed, linked charm muons", fontsize=11.5)
    ax1.set_title("Tagged yield is LINEAR in the efficiency", fontsize=11.5)
    ax2.set_ylabel(r"$\sigma(A_c)$", fontsize=12)
    ax2.set_yscale("log")
    ax2.set_title(r"Statistical reach: $\sigma(A_c)\propto 1/\sqrt{\varepsilon}$",
                  fontsize=11.5)

    fig.suptitle("Main result vs the identification + linking efficiency  "
                 f"(Run 4, {F.LUMI_RUN4:.0f} fb$^{{-1}}$, $F_{{\\rm flux}}$={F.F_FLUX:.2f}, "
                 f"$F_{{\\rm optics}}$={F.F_OPTICS:.1f})\n"
                 "The AHCAL adds ~5.8x the 3DCAL target mass, but its coarse 4x4 cm$^2$ "
                 "granularity should be read as a LOWER $\\varepsilon$ — compare curves at "
                 "different $\\varepsilon$, not at the same one", fontsize=11)
    _emit(fig, pdf, "20_efficiency_scan")


def page_offaxis(pdf, d, cone):
    """
    Off-axis flux: the single largest unquantified factor.

    The flux set is evaluated essentially ON-AXIS (detector axis at
    (1 cm, -3.3 cm) from the LoS, arXiv:2506.13889 Sec. 2).  Going off-axis is
    NOT a simple suppression -- the LHC lattice sweeps mu+/mu- in opposite
    directions, and FLUKA studies report the rate rising by up to an order of
    magnitude in some directions beyond ~1 m.  This page spans the plausible
    range so the answer can be read off once the real fluence is supplied.
    """
    import geometry as G
    ff = np.logspace(np.log10(0.2), np.log10(12.0), 40)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.4))
    palette = ["#4c72b0", "#dd8452", "#c44e52"]
    cols = {nm: palette[i % len(palette)] for i, nm in enumerate(G.SCENARIOS)}
    for nm, t in G.SCENARIOS.items():
        tag = [F.chain_scenario(d, t, cone=cone, f_flux=float(f))["n_tag"] for f in ff]
        sig = [F.chain_scenario(d, t, cone=cone, f_flux=float(f))["sigma_A"] for f in ff]
        ax1.plot(ff, tag, "-", lw=2.2, color=cols[nm], label=nm.replace("_", " "))
        ax2.plot(ff, sig, "-", lw=2.2, color=cols[nm], label=nm.replace("_", " "))
    for ax in (ax1, ax2):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.axvline(1.0, color="black", ls="--", lw=1.4)
        ax.axvspan(0.2, 1.0, color="gray", alpha=0.10)
        ax.axvspan(1.0, 12.0, color="green", alpha=0.07)
        ax.set_xlabel(r"off-axis flux factor $F_{\rm flux}$ (1 = on-axis)", fontsize=11)
        ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")
    ax1.set_ylabel("tagged charm muons (Run 4)", fontsize=11)
    ax1.set_title("Yield scales linearly with the off-axis flux", fontsize=11.5)
    ax2.set_ylabel(r"$\sigma(A_c)$", fontsize=11)
    ax2.set_title(r"$\sigma(A_c)\propto 1/\sqrt{F_{\rm flux}}$", fontsize=11.5)
    fig.suptitle("Off-axis flux — the largest unquantified factor. FLUKA reports the rate "
                 "RISING off-axis in some\ndirections (green), so this is not necessarily a "
                 "penalty. REQUIRED INPUT: fluence at the 3DCAL position.",
                 fontsize=11.5)
    _emit(fig, pdf, "19_offaxis_flux")


def main():
    global PNG_DIR
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", default="/eos/home-e/evilla/faser/fasercal_dis_v1.npz")
    ap.add_argument("--npz-pc", default="/eos/home-e/evilla/faser/fasercal_dis_pc.npz")
    ap.add_argument("--output", default=config.get("report_pdf",
                    "/eos/home-e/evilla/faser/reports/fasercal_chain.pdf"))
    ap.add_argument("--png-dir", default=None)
    args = ap.parse_args()
    PNG_DIR = args.png_dir

    d = load(args.npz)
    m = F.chain(d)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with PdfPages(args.output) as pdf:
        page_conditions(pdf, d)
        m_ah = F.chain_scenario(d, 0.1, cone=F.calibrate_cone(d, __import__("geometry").cdr_config(0.1)),
                                f_flux=F.F_FLUX, include_ahcal=True)
        page_cutflow(pdf, d, m, m_ah)
        page_funnel(pdf, m)
        page_muon_spectrum(pdf, d)
        page_angular_acceptance(pdf, d)
        page_mass_and_params(pdf, d)
        import geometry as G
        cone = F.calibrate_cone(d, G.configure(0.0))
        page_tungsten_scenarios(pdf, d, cone)
        page_tungsten_scan(pdf, d, cone)
        page_efficiency(pdf, d, cone)
        page_offaxis(pdf, d, cone)
        if os.path.exists(args.npz_pc):
            page_ic(pdf, d, load(args.npz_pc))

    print(f"Done -> {args.output}")
    for k in ["n_dis", "n_charm", "n_semi", "n_punch", "n_acc", "n_tag"]:
        print(f"  {k:10s} {m[k]:>14,.1f}")
    print(f"  f_charm {100*m['f_charm']:.3f}%  f_semi {100*m['f_semi']:.1f}%  "
          f"f_punch {100*m['f_punch']:.1f}%  overall {100*m['f_overall']:.4f}%")
    print(f"  purity_q {100*m['purity_q']:.1f}%  sigma(A_c) +-{m['sigma_A']:.3f}  "
          f"raw MC {m['n_tag_raw']}")


if __name__ == "__main__":
    main()
