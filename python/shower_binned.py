#!/usr/bin/env python3
"""
Shower the FIXED-ENERGY-BIN POWHEG samples and write the per-event cache.

This replaces shower_dis.py for the binned production.  Two things change, and
both are improvements:

  1. E_mu IS EXACT.  Each bin is generated with a monochromatic beam
     (`fixed_lepton_beam 1`, x_lepton = 1), so the incoming muon energy is the
     bin energy by construction rather than something sampled inside POWHEG.

  2. THE HADRONIC SUM IS CLEAN.  With x_lepton = 1 the lepton carries the whole
     beam momentum, so there is NO lepton-side beam remnant.  The few-GeV
     forward contamination that broke the particle-level hadronic energy in the
     flux-as-beam-PDF production is simply absent.  This module therefore sums
     the hadronic final state DIRECTLY from the particles, and checks energy
     conservation event by event -- something that failed in 80% of events
     before (see docs/XRECO.md).

     Both are stored: `e_had` from the particle sum, and `e_had_cons` from
     momentum conservation (nu + M), so the two can be compared.  If they agree
     the workaround can be retired.

WEIGHTS.  POWHEG returns a cross section per bin for a monochromatic beam; the
flux enters afterwards as the per-bin integral Int f(x) dx from
python/flux_bins.py.  Summing over bins reproduces the beam-PDF normalisation by
construction, so absolute rates are unchanged.

Run inside the generator's environment (it provides pythia8):
  source .../muondisgenerator/scripts/env.sh
  python3 shower_binned.py --indir <production_binned> --npz out.npz
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pythia8

M_N = 0.9383
SAMPLES = ["mum_proton", "mum_neutron", "mup_proton", "mup_neutron"]

KEYS = ["w_raw", "is_neutron", "e_in", "e_in_bin", "p_out", "theta_mu", "q2",
        "e_had", "px_had", "py_had", "pz_had",
        "e_had_cons", "e_nu_all", "e_nu_charm", "e_em", "e_neuhad", "e_chhad",
        "ebal", "n_charm", "n_semilep_mu", "mu2_q", "mu2_p", "mu2_theta",
        "charm_sign_mu", "ibin"]


def is_charm_hadron(pid):
    a = abs(pid)
    return a >= 100 and ((a // 100) % 10 == 4 or (a // 1000) % 10 == 4)


def charm_ancestor(ev, i):
    seen, st = set(), [i]
    while st:
        j = st.pop()
        if j in seen or j <= 0:
            continue
        seen.add(j)
        p = ev[j]
        if j != i and is_charm_hadron(p.id()):
            return j
        m1, m2 = p.mother1(), p.mother2()
        if m1 > 0:
            st.append(m1)
        if m2 > 0 and m2 != m1:
            st.append(m2)
    return -1


def xsec_from_log(rundir):
    """Total cross section [pb] from the POWHEG stage-4 log."""
    for name in ("s4.log", "s3.log", "s2.log"):
        p = os.path.join(rundir, name)
        if not os.path.exists(p):
            continue
        tot = None
        for line in open(p, errors="ignore"):
            if "total (btilde+remnants) cross section" in line.lower() or \
               "total cross section" in line.lower():
                for tok in line.replace("+-", " ").split():
                    try:
                        tot = float(tok.replace("D", "E"))
                        break
                    except ValueError:
                        continue
        if tot:
            return tot
    # fall back on the btilde + remnant estimates from stage 1
    btl = rem = 0.0
    for line in open(os.path.join(rundir, "s1_3.log"), errors="ignore"):
        if "btilde: estimated absolute value cross section" in line:
            btl = float(line.split(":")[-1].split("+-")[0])
        if "remn: estimated absolute value cross section" in line:
            rem = float(line.split(":")[-1].split("+-")[0])
    return btl + rem


def shower_one(lhe, n_events, is_neutron, e_bin, ibin, w_bin):
    py = pythia8.Pythia("", False)
    for c in ["Beams:frameType = 4", f"Beams:LHEF = {lhe}",
              "Main:numberOfEvents = 0", "Print:quiet = on",
              "Init:showChangedSettings = off", "Next:numberShowEvent = 0",
              "Stat:showProcessLevel = off", "POWHEG:nFinal = 2"]:
        py.readString(c)
    if not py.init():
        print(f"  [warn] init failed: {lhe}", file=sys.stderr)
        return {k: [] for k in KEYS}

    cols = {k: [] for k in KEYS}
    n = 0
    while n < n_events:
        if not py.next():
            if py.infoPython().atEndOfFile():
                break
            continue
        proc = py.process
        in_mu = out_mu = None
        for q in proc:
            if abs(q.id()) != 13:
                continue
            if q.statusAbs() == 21 and in_mu is None:
                in_mu = q
            elif q.status() > 0 and out_mu is None:
                out_mu = q
        if in_mu is None or out_mu is None:
            continue
        q2 = -((in_mu.e() - out_mu.e())**2 - (in_mu.px() - out_mu.px())**2
               - (in_mu.py() - out_mu.py())**2 - (in_mu.pz() - out_mu.pz())**2)
        vi = np.array([in_mu.px(), in_mu.py(), in_mu.pz()])
        vo = np.array([out_mu.px(), out_mu.py(), out_mu.pz()])
        ni, no = np.linalg.norm(vi), np.linalg.norm(vo)
        th = np.arccos(np.clip(np.dot(vi, vo) / (ni * no), -1, 1)) if ni * no else np.nan

        ev = py.event
        charm = {}
        for i in range(ev.size()):
            p = ev[i]
            if not is_charm_hadron(p.id()):
                continue
            dl = [d for d in p.daughterList() if d > 0]
            if not dl or any(is_charm_hadron(ev[d].id()) for d in dl):
                continue
            charm[i] = {"sign": 1 if p.id() > 0 else -1, "mu": None}

        si, bp = -1, -1.0
        for i in range(ev.size()):
            p = ev[i]
            if not p.isFinal() or abs(p.id()) != 13:
                continue
            a = charm_ancestor(ev, i)
            if a >= 0 and a in charm:
                r = charm[a]
                if r["mu"] is None or p.pAbs() > r["mu"].pAbs():
                    r["mu"] = p
            elif p.pAbs() > bp:
                bp, si = p.pAbs(), i

        # --- hadronic system, straight from the particles (now legitimate) ---
        eh = px = py_ = pz = 0.0
        e_nu_all = e_nu_ch = e_em = e_neu = e_ch = 0.0
        for i in range(ev.size()):
            p = ev[i]
            if not p.isFinal() or i == si:
                continue
            ida = abs(p.id())
            if ida in (12, 14, 16):
                e_nu_all += p.e()
                a = charm_ancestor(ev, i)
                if a >= 0 and a in charm:
                    e_nu_ch += p.e()
                continue
            eh += p.e(); px += p.px(); py_ += p.py(); pz += p.pz()
            if ida in (22, 11):
                e_em += p.e()
            elif p.isCharged():
                e_ch += p.e()
            else:
                e_neu += p.e()

        nu = in_mu.e() - out_mu.e()
        n_semi = sum(1 for r in charm.values() if r["mu"] is not None)
        lead, lsign = None, np.nan
        for r in charm.values():
            if r["mu"] is not None and (lead is None or r["mu"].pAbs() > lead.pAbs()):
                lead, lsign = r["mu"], r["sign"]

        row = dict(
            w_raw=w_bin, is_neutron=float(is_neutron),
            e_in=in_mu.e(), e_in_bin=e_bin, p_out=out_mu.pAbs(),
            theta_mu=th, q2=q2,
            e_had=eh, px_had=px, py_had=py_, pz_had=pz,
            e_had_cons=nu + M_N,
            e_nu_all=e_nu_all, e_nu_charm=e_nu_ch,
            e_em=e_em, e_neuhad=e_neu, e_chhad=e_ch,
            ebal=in_mu.e() + M_N - out_mu.e() - eh - e_nu_all,
            n_charm=float(len(charm)), n_semilep_mu=float(n_semi),
            mu2_q=(1.0 if lead is not None and lead.id() < 0
                   else (-1.0 if lead is not None else np.nan)),
            mu2_p=lead.pAbs() if lead is not None else np.nan,
            mu2_theta=(np.arctan2(np.hypot(lead.px(), lead.py()), lead.pz())
                       if lead is not None else np.nan),
            charm_sign_mu=float(lsign), ibin=float(ibin),
        )
        for k in KEYS:
            cols[k].append(row[k])
        n += 1
    py.stat()
    return cols


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", required=True)
    ap.add_argument("--npz", required=True)
    ap.add_argument("--nevents", type=int, default=100000)
    args = ap.parse_args()

    bins = json.load(open(os.path.join(args.indir, "bins.json")))
    ecent = np.array(bins["e_centre"])
    w_mu, w_mubar = np.array(bins["w_mu"]), np.array(bins["w_mubar"])

    cols = {k: [] for k in KEYS}
    for sample in SAMPLES:
        is_n = "neutron" in sample
        flux = w_mubar if sample.startswith("mup") else w_mu
        for ibin in range(len(ecent)):
            d = os.path.join(args.indir, sample, f"bin{ibin:02d}")
            lhe = os.path.join(d, "pwgevents-0001.lhe")
            if not os.path.exists(lhe):
                continue
            sig = xsec_from_log(d)
            nev = sum(1 for _ in open(lhe) if "<event>" in _)
            if nev == 0 or sig is None:
                continue
            # expected events for this bin, divided evenly over the MC events
            w_bin = flux[ibin] * sig / nev
            c = shower_one(lhe, args.nevents, is_n, ecent[ibin], ibin, w_bin)
            for k in KEYS:
                cols[k].extend(c[k])
            print(f"  {sample} bin{ibin:02d}  E={ecent[ibin]:7.1f} GeV  "
                  f"sigma={sig:.4g}  N={len(c['w_raw'])}  w={w_bin:.4g}", flush=True)

    arr = {k: np.array(cols[k], dtype=float) for k in KEYS}
    n = len(arr["w_raw"])
    print(f"\nTotal {n} events, expected yield {arr['w_raw'].sum():.4g}")
    if n:
        bal = arr["ebal"]
        print(f"Energy balance E_in + M - E_mu' - E_had - E_nu:  "
              f"median {np.median(bal):+.3f} GeV, {100*(bal < -1).mean():.1f}% below -1 GeV")
        print("  (the flux-as-beam-PDF production had 80% of events negative)")
        r = arr["e_had"] / np.maximum(arr["e_had_cons"], 1e-9)
        print(f"E_had(particles) / E_had(conservation): median {np.median(r):.4f}")
    os.makedirs(os.path.dirname(args.npz), exist_ok=True)
    np.savez(args.npz, **arr)
    print(f"Cached -> {args.npz}")


if __name__ == "__main__":
    main()
