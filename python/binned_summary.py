#!/usr/bin/env python3
"""
Every number the binned production changed, in one place.

Run after shower_binned.py.  Prints the closure against the old flux-as-beam-PDF
production, the validation tests that justify the change, the chain yields, and
the asymmetry with its EFFECTIVE (not raw) statistics.

  python3 binned_summary.py [--new <npz>] [--old <npz>]
"""
import argparse
import sys
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fasercal_chain as F      # noqa: E402
import xreco as X               # noqa: E402

M_N = 0.9383
NEW = "/eos/home-e/evilla/faser/fasercal_binned_v1.npz"
OLD = "/eos/home-e/evilla/faser/fasercal_dis_v1_had.npz"


def load(p):
    return {k: v for k, v in np.load(p).items()}


def n_eff(w):
    """Effective entries.  For weighted MC with negative weights this is far
    below the raw count, and it is what a statistical error must be built on."""
    s2 = (w ** 2).sum()
    return w.sum() ** 2 / s2 if s2 > 0 else 0.0


def wquant(v, w, q):
    i = np.argsort(v)
    v, w = v[i], w[i]
    c = np.cumsum(w) - 0.5 * w
    return np.interp(q, c / w.sum(), v)


def hadronic_W(d):
    """Invariant mass of the hadronic system: W^2 = M^2 + 2 M nu - Q^2."""
    nu = d["e_in"] - d["p_out"]
    return np.sqrt(np.maximum(M_N ** 2 + 2 * M_N * nu - d["q2"], 0.0))


def sec(title):
    print(f"\n{'=' * 74}\n{title}\n{'=' * 74}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--new", default=NEW)
    ap.add_argument("--old", default=OLD)
    a = ap.parse_args()
    new = load(a.new)
    old = load(a.old) if os.path.exists(a.old) else None

    sec("1. GENERATION CLOSURE  (the two strategies share only flux and PDF set)")
    for tag, d in (("flux-as-beam-PDF", old), ("binned", new)):
        if d is None:
            continue
        w = F.event_weight(d)
        c = d["n_charm"] >= 1
        print(f"  {tag:18s} N={len(w):7d}  yield={w.sum():.4g}  "
              f"<E_mu>={(w * d['e_in']).sum() / w.sum():7.1f} GeV  "
              f"charm frac={w[c].sum() / w.sum():.4f}")
    if old is not None:
        r = F.event_weight(new).sum() / F.event_weight(old).sum()
        print(f"  --> total DIS yield ratio binned/old = {r:.3f}")

    sec("2. THE REMNANT IS GONE  (energy conservation, particle-level)")
    bal = new["ebal"]
    ratio = new["e_had"] / np.maximum(new["e_had_cons"], 1e-9)
    t = X.truth(new)
    hi = t["x"] > 0.2
    print(f"  energy balance      median {np.median(bal):+.3f} GeV, "
          f"{100 * (bal < -1).mean():.1f}% below -1 GeV   (was 80% negative)")
    print(f"  E_had(part)/E_had(cons)   median {np.median(ratio):.4f}")
    print(f"    restricted to x > 0.2   median {np.median(ratio[hi]):.4f}   "
          f"(the old sample reached 2.5 here)")

    sec("3. CHARM RESPECTS ITS KINEMATIC THRESHOLD  (needs W > 2 m_D = 3.73 GeV)")
    for tag, d in (("flux-as-beam-PDF", old), ("binned", new)):
        if d is None:
            continue
        W = hadronic_W(d)
        c = np.isfinite(d["charm_sign_mu"]) | (d["n_charm"] >= 1)
        for lab, sel in ((f"{tag}, all charm", c),
                         (f"{tag}, E_mu < 20 GeV", c & (d["e_in"] < 20))):
            if sel.sum() == 0:
                print(f"  {lab:34s} N=0")
                continue
            print(f"  {lab:34s} N={sel.sum():6d}  min W={W[sel].min():5.2f}  "
                  f"frac(W<3.73)={(W[sel] < 3.73).mean():.3f}")

    sec("4. EFFECTIVE STATISTICS  (POWHEG negative weights, not a bug)")
    for tag, d in (("flux-as-beam-PDF", old), ("binned", new)):
        if d is None:
            continue
        w = F.event_weight(d)
        c = d["n_charm"] >= 1
        print(f"  {tag}")
        for lab, sel in (("inclusive", np.ones(len(w), bool)), ("charm", c)):
            ne = n_eff(w[sel])
            print(f"    {lab:10s} N_raw={sel.sum():7d}  N_eff={ne:9.0f}  "
                  f"retained {ne / max(sel.sum(), 1):.3f}  "
                  f"frac(w<0)={(w[sel] < 0).mean():.3f}")

    sec("5. TAGGED CHARM ASYMMETRY  (error from N_eff, NOT from N_raw)")
    for tag, d in (("flux-as-beam-PDF", old), ("binned", new)):
        if d is None:
            continue
        w = F.event_weight(d)
        tt = X.truth(d)
        s = d["charm_sign_mu"]
        m = np.isfinite(s)
        print(f"  {tag}")
        for lab, sel in (("all semileptonic", m),
                         ("  + p_mu > 5 GeV", m & (d["mu2_p"] > 5)),
                         ("  + x_truth > 0.2", m & (d["mu2_p"] > 5) & (tt["x"] > 0.2))):
            ww, ss = w[sel], s[sel]
            if ww.sum() == 0 or sel.sum() == 0:
                print(f"    {lab:20s} -")
                continue
            A = (ww * ss).sum() / ww.sum()
            ne = n_eff(ww)
            e = 1 / np.sqrt(ne) if ne > 0 else np.inf
            print(f"    {lab:20s} N_raw={sel.sum():6d}  N_eff={ne:7.0f}  "
                  f"A={A:+.4f} +- {e:.4f}  ({abs(A) / e:.1f} sigma)")

    sec("6. x RECONSTRUCTION  (particle sum vs the retired workaround)")
    w = F.event_weight(new) * (new["n_charm"] >= 1)
    res = {s: X.reconstruct(new, np.random.default_rng(3), had_source=s)
           for s in ("particles", "conservation")}
    print(f"  {'method':26s} {'particles':>20s} {'conservation':>20s}   (bias / half-width)")
    for meth in res["particles"]:
        line = f"  {meth:26s}"
        for s in ("particles", "conservation"):
            x = res[s][meth]
            ok = np.isfinite(x) & (t["x"] > 0) & (w > 0)
            rel = (x[ok] - t["x"][ok]) / t["x"][ok]
            ww = w[ok]
            hw = 0.5 * (wquant(rel, ww, 0.84) - wquant(rel, ww, 0.16))
            line += f" {wquant(rel, ww, 0.5):+7.3f} /{hw:7.3f}"
        print(line)

    sec("7. CHAIN YIELDS  (Run 4, 780 /fb, 3DCal only)")
    for tw in (0.1, 0.5, 1.0):
        r = F.chain_scenario(new, tw)
        print(f"  W = {tw * 10:.0f} mm/module:  DIS={r['n_dis']:.4g}  "
              f"charm={r['n_charm']:.4g}  semilep={r['n_semi']:.4g}  "
              f"tagged={r['n_tag']:.4g}  A={r['A_meas']:+.4f}")

    sec("8. PER-BIN ALLOCATION  (where the signal actually is)")
    w = F.event_weight(new)
    s = new["charm_sign_mu"]
    m = np.isfinite(s)
    ib = new["ibin"].astype(int)
    import flux_bins as B
    ec = B.flux_weights()["e_centre"]
    tot = w[m].sum()
    print(f"  {'bin':>4} {'E [GeV]':>9} {'charm yield':>12} {'N_raw':>8} "
          f"{'N_eff':>8} {'share':>7}")
    for i in range(len(ec)):
        sel = m & (ib == i)
        if sel.sum() == 0:
            continue
        print(f"  {i:4d} {ec[i]:9.1f} {w[sel].sum():12.4g} {sel.sum():8d} "
              f"{n_eff(w[sel]):8.0f} {w[sel].sum() / tot:7.3f}")
    print()


if __name__ == "__main__":
    main()
