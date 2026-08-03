#!/usr/bin/env python3
"""Report pages for the Bjorken-x reconstruction study (see docs/XRECO.md)."""
import argparse, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
import xreco as X, fasercal_chain as F

METHODS = ["lepton-only", "Jacquet-Blondel", "Sigma", "double-angle"]
PNG_DIR = None


def _emit(fig, pdf, name):
    pdf.savefig(fig)
    if PNG_DIR:
        os.makedirs(PNG_DIR, exist_ok=True)
        fig.savefig(os.path.join(PNG_DIR, f"{name}.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


def page_migration(pdf, d, t, w, r):
    """Truth-x -> reco-x migration matrix for each method."""
    bins = np.logspace(-3, 0, 40)
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.8))
    for ax, m in zip(axes, METHODS):
        xr = np.nan_to_num(r[m], nan=1e-9)
        ok = (t["x"] > 0) & (xr > 0) & (d["n_charm"] >= 1)
        ax.hist2d(t["x"][ok], xr[ok], bins=[bins, bins], weights=w[ok],
                  norm=LogNorm(), cmap="viridis")
        ax.plot([1e-3, 1], [1e-3, 1], "w--", lw=1.4)
        for c in (0.2, 0.4):
            ax.axhline(c, color="#ff6666", lw=0.9, ls=":")
            ax.axvline(c, color="#ff6666", lw=0.9, ls=":")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$x_{\rm true}$", fontsize=11)
        if m == METHODS[0]:
            ax.set_ylabel(r"$x_{\rm reco}$", fontsize=11)
        ax.set_title(m, fontsize=11.5)
    fig.suptitle("Truth-$x$ $\\to$ reco-$x$ migration, charm events. "
                 "Dashed = perfect; dotted = the $x\\geq0.2$ and $0.4$ signal cuts.\n"
                 "Lepton-only shows almost NO correlation — $x$ information is destroyed. "
                 "$\\Sigma$ retains a clear diagonal.",
                 fontsize=12)
    fig.tight_layout(); _emit(fig, pdf, "30_migration")


CUT = 0.2       # only x >= 0.2: above 0.4 the perturbative sample has 12 raw
                # events, so nothing there is statistically meaningful.


def _effpur(d, t, w, m, cut=CUT, res=None, ein_known=True, seed=3):
    r = X.reconstruct(d, np.random.default_rng(seed), res=res, ein_known=ein_known)
    if not np.isfinite(r[m]).any():
        return np.nan
    T = t["x"] >= cut
    R = np.nan_to_num(r[m], nan=-1) >= cut
    eff = w[T & R].sum() / w[T].sum()
    pur = w[T & R].sum() / max(w[R].sum(), 1e-30)
    return eff * pur


def page_effpur(pdf, d, t, w):
    """Signal retention vs the hadronic angular resolution, and vs E_in knowledge."""
    cols = dict(zip(METHODS, ["#c44e52", "#dd8452", "#4c72b0", "#55a868"]))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.8, 5.4))

    sths = np.array([0., 5., 10., 20., 35., 50., 75., 100.]) * 1e-3
    for m in METHODS:
        v = [_effpur(d, t, w, m, res=dict(sigma_theta_had=s)) for s in sths]
        ax1.plot(sths * 1e3, v, "o-", color=cols[m], lw=2.2, label=m)
    ax1.axvline(20, color="black", ls="--", lw=1.2)
    ax1.set_xlabel(r"hadronic angular resolution $\sigma_{\theta_h}$ [mrad]", fontsize=11)
    ax1.set_ylabel(r"efficiency $\times$ purity", fontsize=11)
    ax1.set_title(r"vs hadronic angular resolution", fontsize=11.5)
    ax1.grid(alpha=0.3); ax1.legend(fontsize=9)

    xs = np.arange(len(METHODS)); wd = 0.36
    a = [_effpur(d, t, w, m, ein_known=True) for m in METHODS]
    b = [_effpur(d, t, w, m, ein_known=False) for m in METHODS]
    ax2.bar(xs - wd/2, a, wd, label=r"$E_{\rm in}$ measured", color="#4c72b0",
            edgecolor="black", lw=0.6)
    ax2.bar(xs + wd/2, np.nan_to_num(b), wd, label=r"$E_{\rm in}$ NOT measured",
            color="#dd8452", edgecolor="black", lw=0.6)
    for i, v in enumerate(b):
        if not np.isfinite(v):
            ax2.text(i + wd/2, 0.02, "N/A", ha="center", fontsize=9,
                     rotation=90, fontweight="bold")
    ax2.set_xticks(xs); ax2.set_xticklabels([m.replace("-", "-\n") for m in METHODS],
                                            fontsize=9)
    ax2.set_ylabel(r"efficiency $\times$ purity", fontsize=11)
    ax2.set_title("vs whether the incoming muon energy is measured", fontsize=11.5)
    ax2.grid(axis="y", alpha=0.3); ax2.legend(fontsize=9)

    fig.suptitle(rf"Large-$x$ signal retention, $x \geq {CUT}$.  "
                 r"$\Sigma$ is flat in $\sigma_{\theta_h}$ (it uses $E+p_z$) AND "
                 r"unaffected by $E_{\rm in}$ (it reconstructs it)." "\n"
                 "If the incoming muon is not measured — only a downstream spectrometer — "
                 "lepton-only and double-angle are not available at all.", fontsize=11)
    fig.tight_layout(); _emit(fig, pdf, "31_effpur")


def main():
    global PNG_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="/eos/home-e/evilla/faser/fasercal_dis_v1_had.npz")
    ap.add_argument("--output", default="/eos/home-e/evilla/faser/reports/fasercal_xreco.pdf")
    ap.add_argument("--png-dir", default=None)
    a = ap.parse_args(); PNG_DIR = a.png_dir
    d = {k: v for k, v in np.load(a.npz).items()}
    t = X.truth(d); w = F.event_weight(d) * (d["n_charm"] >= 1)
    r = X.reconstruct(d, np.random.default_rng(3))
    os.makedirs(os.path.dirname(a.output), exist_ok=True)
    with PdfPages(a.output) as pdf:
        page_migration(pdf, d, t, w, r)
        page_effpur(pdf, d, t, w)
    print("Done ->", a.output)


if __name__ == "__main__":
    main()
