#!/usr/bin/env python3
"""
Kinematics of BOTH charm hadrons, not just the tagged one.

WHY.  Everything in the analysis so far is driven by the LEADING semileptonic
muon, so the softer charm hadron is silently dropped from the charge tag (see
docs/CHARM_ASYMMETRY.md).  That is defensible for the tag itself, but it leaves
an open question: is the second charm hadron a usable SIGNATURE, or is it
invisible?  In the flavour-excitation / intrinsic-charm topology the struck c is
hard while its partner cbar hadronises with the target remnant, so the two are
expected to be very different -- and the softer one may be below anything the
detector can resolve.

This module records, per event and per weakly-decaying charm hadron:
    p, theta, pid, charge sign, whether it decayed to a muon, that muon's p
sorted by |p|, so the leading and sub-leading charm hadron can be compared.

It also records the proper decay length boosted into the lab,
    L = beta gamma c tau ,
which is what decides whether the decay is a resolvable displaced vertex in the
3DCal's 1 cm voxels rather than a point-like deposit.

Run inside the generator environment (provides pythia8):
  python3 charm_pair.py --indir <production_binned> --npz out.npz [--nevents N]
"""
import argparse
import json
import os
import sys

import numpy as np
import pythia8

M_N = 0.9383
SAMPLES = ["mum_proton", "mum_neutron", "mup_proton", "mup_neutron"]

# c tau in mm, PDG, for the weakly-decaying charm hadrons that matter
CTAU_MM = {411: 0.3118,    # D+
           421: 0.1237,    # D0
           431: 0.1497,    # Ds+
           4122: 0.0602,   # Lambda_c+
           4132: 0.0453,   # Xi_c0
           4232: 0.1320,   # Xi_c+
           4332: 0.0807}   # Omega_c0

KEYS = ["w_raw", "is_neutron", "e_in", "ibin", "q2", "p_out", "theta_mu",
        "n_charm",
        # leading charm hadron (by |p|)
        "c1_p", "c1_theta", "c1_pid", "c1_sign", "c1_mu_p", "c1_len_mm",
        # sub-leading
        "c2_p", "c2_theta", "c2_pid", "c2_sign", "c2_mu_p", "c2_len_mm"]


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


def decay_length_mm(p_gev, m_gev, pid):
    """beta gamma c tau in the lab, in mm.  betagamma = p/m."""
    ct = CTAU_MM.get(abs(int(pid)))
    if ct is None or m_gev <= 0:
        return np.nan
    return (p_gev / m_gev) * ct


def shower_one(lhe, n_events, is_neutron, ibin, flux_i):
    py = pythia8.Pythia("", False)
    for c in ["Beams:frameType = 4", f"Beams:LHEF = {lhe}",
              "Main:numberOfEvents = 0", "Print:quiet = on",
              "Init:showChangedSettings = off", "Next:numberShowEvent = 0",
              "Stat:showProcessLevel = off", "POWHEG:nFinal = 2",
              "LesHouches:matchInOut = off"]:
        py.readString(c)
    cols = {k: [] for k in KEYS}
    if not py.init():
        return cols

    n = 0
    while n < n_events:
        if not py.next():
            if py.infoPython().atEndOfFile():
                break
            continue
        proc, ev = py.process, py.event
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
        th_mu = (np.arccos(np.clip(np.dot(vi, vo) / (ni * no), -1, 1))
                 if ni * no else np.nan)

        # weakly-decaying charm hadrons
        charm = {}
        for i in range(ev.size()):
            p = ev[i]
            if not is_charm_hadron(p.id()):
                continue
            dl = [d for d in p.daughterList() if d > 0]
            if not dl or any(is_charm_hadron(ev[d].id()) for d in dl):
                continue
            charm[i] = dict(p=p.pAbs(),
                            theta=np.arctan2(np.hypot(p.px(), p.py()), p.pz()),
                            pid=p.id(),
                            sign=1 if p.id() > 0 else -1,
                            m=p.m(), mu=np.nan)

        # attach each semileptonic muon to its parent charm hadron
        for i in range(ev.size()):
            p = ev[i]
            if not p.isFinal() or abs(p.id()) != 13:
                continue
            a = charm_ancestor(ev, i)
            if a in charm:
                if not np.isfinite(charm[a]["mu"]) or p.pAbs() > charm[a]["mu"]:
                    charm[a]["mu"] = p.pAbs()

        order = sorted(charm.values(), key=lambda r: -r["p"])
        row = dict(w_raw=flux_i * py.infoPython().weight(),
                   is_neutron=float(is_neutron), e_in=in_mu.e(),
                   ibin=float(ibin), q2=q2, p_out=out_mu.pAbs(),
                   theta_mu=th_mu, n_charm=float(len(charm)))
        for k, tag in ((0, "c1"), (1, "c2")):
            if k < len(order):
                r = order[k]
                row[f"{tag}_p"] = r["p"]
                row[f"{tag}_theta"] = r["theta"]
                row[f"{tag}_pid"] = float(r["pid"])
                row[f"{tag}_sign"] = float(r["sign"])
                row[f"{tag}_mu_p"] = r["mu"]
                row[f"{tag}_len_mm"] = decay_length_mm(r["p"], r["m"], r["pid"])
            else:
                for s in ("p", "theta", "pid", "sign", "mu_p", "len_mm"):
                    row[f"{tag}_{s}"] = np.nan
        for k in KEYS:
            cols[k].append(row[k])
        n += 1
    py.stat()
    return cols


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", required=True, nargs="+")
    ap.add_argument("--npz", required=True)
    ap.add_argument("--nevents", type=int, default=1000000,
                    help="cap per LHE file; weights use the number ACTUALLY "
                         "showered, so a cap only costs statistics")
    a = ap.parse_args()

    bins = json.load(open(os.path.join(a.indir[0], "bins.json")))
    ecent = np.array(bins["e_centre"])
    w_mu, w_mubar = np.array(bins["w_mu"]), np.array(bins["w_mubar"])

    cols = {k: [] for k in KEYS}
    for sample in SAMPLES:
        is_n = "neutron" in sample
        flux = w_mubar if sample.startswith("mup") else w_mu
        for ibin in range(len(ecent)):
            dirs = [os.path.join(r, sample, f"bin{ibin:02d}") for r in a.indir]
            lhes = [os.path.join(x, "pwgevents-0001.lhe") for x in dirs]
            lhes = [x for x in lhes if os.path.exists(x)]
            if not lhes:
                continue
            got = {k: [] for k in KEYS}
            for lhe in lhes:
                c = shower_one(lhe, a.nevents, is_n, ibin, 1.0)
                for k in KEYS:
                    got[k].extend(c[k])
            nsh = len(got["w_raw"])
            if nsh == 0:
                continue
            # weight normalisation uses the events actually showered
            got["w_raw"] = [w * flux[ibin] / nsh for w in got["w_raw"]]
            for k in KEYS:
                cols[k].extend(got[k])
            print(f"  {sample} bin{ibin:02d} E={ecent[ibin]:7.1f} N={nsh}",
                  flush=True)

    arr = {k: np.array(cols[k], dtype=float) for k in KEYS}
    os.makedirs(os.path.dirname(a.npz), exist_ok=True)
    np.savez(a.npz, **arr)
    print(f"\n{len(arr['w_raw'])} events -> {a.npz}")


if __name__ == "__main__":
    main()
