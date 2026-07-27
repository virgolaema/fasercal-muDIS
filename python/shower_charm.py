#!/usr/bin/env python3
"""
Shower the muon-DIS POWHEG LHE events with Pythia8 and record, per charm
hadron, everything the FASERcal charm/anticharm-tagging toy needs.

Physics of the tag
------------------
The observable that separates charm from anticharm is the CHARGE of the muon
from the semileptonic decay of the charmed hadron:

    c  -> s W+  ,  W+ -> mu+ nu     (charm     -> mu+ ,  muon id -13)
    c~ -> s~ W- ,  W- -> mu- nu~    (anticharm -> mu- ,  muon id +13)

so the decay-muon charge == the charm sign (both +1 for c, -1 for c~).  A
3D scintillator like FASERcal has NO magnetic field, so it cannot sign the
muon itself; the charge must come from a downstream magnetised spectrometer.
This step produces the truth inputs (charm sign, decay-muon momentum/angle,
charm-hadron flight length); the toy detector response is applied later in
make_report.py so the detector assumptions stay isolated and easy to vary.

For each charm hadron we store:
  sign      : +1 (contains c) / -1 (contains c~), from the PDG-code sign
  pid       : the charmed-hadron PDG id (411, 421, 431, 4122, ...)
  E, p      : lab energy and |p| of the charm hadron [GeV]
  Ldecay    : flight length production->decay vertex [mm]  (the FASERnu handle)
  has_mu    : whether the decay chain yields a muon
  mu_q      : decay-muon charge (+-1), NaN if has_mu is False
  mu_p      : decay-muon lab |p| [GeV]
  mu_theta  : decay-muon polar angle w.r.t. beam z [rad]
  mu_ptrel  : decay-muon pT relative to the parent charm-hadron direction [GeV]
  w         : tungsten-combined event weight (nominal POWHEG weight)

Tungsten-combined (w_p=74/184, w_n=110/184) over the four beam/target samples,
matching the FASERnu study.  Run inside an LCG view that provides pythia8, e.g.
  source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
  python3 shower_charm.py --nevents 20000
"""
import argparse
import sys
from pathlib import Path

import numpy as np

import pythia8
import config

W_PROTON  = 74.0 / 184.0
W_NEUTRON = 110.0 / 184.0
PROD      = config.get("production_dir")
SAMPLES   = ["mum_proton", "mum_neutron", "mup_proton", "mup_neutron"]


def signal_group(prod_name):
    """(label, LHE path) for a 'production_v1'-layout run (seed1 subdir)."""
    return [(s, f"{PROD}/{prod_name}/{s}/pythia8_seed1/pwgevents.lhe") for s in SAMPLES]


def tung_factor(name):
    return W_NEUTRON if "neutron" in name else W_PROTON


def is_charm_hadron(pid):
    """True if pid is a hadron containing a charm (anti)quark."""
    a = abs(pid)
    if a < 100:
        return False
    return (a // 100) % 10 == 4 or (a // 1000) % 10 == 4


def charm_ancestor(event, i):
    """
    Walk up the mother chain of particle index i; return the index of the
    nearest charmed-hadron ancestor, or -1 if none.  Uses mother1/mother2 so
    both legs of a merge are followed.
    """
    seen = set()
    stack = [i]
    while stack:
        j = stack.pop()
        if j in seen or j <= 0:
            continue
        seen.add(j)
        p = event[j]
        if j != i and is_charm_hadron(p.id()):
            return j
        m1, m2 = p.mother1(), p.mother2()
        if m1 > 0:
            stack.append(m1)
        if m2 > 0 and m2 != m1:
            stack.append(m2)
    return -1


def shower_one(lhe, n_events):
    """Shower up to n_events; return a list of per-charm-hadron dict records."""
    pythia = pythia8.Pythia("", False)
    for cmd in [
        "Beams:frameType = 4",
        f"Beams:LHEF = {lhe}",
        "Main:numberOfEvents = 0",
        "Print:quiet = on",
        "Init:showChangedSettings = off",
        "Init:showChangedParticleData = off",
        "Next:numberShowInfo = 0",
        "Next:numberShowProcess = 0",
        "Next:numberShowEvent = 0",
        "Stat:showProcessLevel = off",
        "POWHEG:nFinal = 2",           # LHE matching (main31.cmnd)
    ]:
        pythia.readString(cmd)
    if not pythia.init():
        print(f"  [warn] Pythia init failed for {lhe}", file=sys.stderr)
        return []

    records = []
    n = 0
    while n < n_events:
        if not pythia.next():
            if pythia.infoPython().atEndOfFile():
                break
            continue
        ev = pythia.event
        w = pythia.infoPython().weight()

        # 1) The weakly-decaying charm hadrons: a charmed hadron whose daughters
        #    are NOT themselves charmed (i.e. the last charm hadron in the chain,
        #    e.g. D0/D+/Ds/Lc after any D* -> D pi cascade).  These are the ones
        #    that decay semileptonically.
        charm_idx = {}      # index -> record
        for i in range(ev.size()):
            p = ev[i]
            if not is_charm_hadron(p.id()):
                continue
            d1, d2 = p.daughterList(), None
            daughters = list(p.daughterList())
            if any(is_charm_hadron(ev[d].id()) for d in daughters if d > 0):
                continue                       # not the final charm hadron
            if not daughters:
                continue                       # no decay recorded (final-state)
            # decay vertex = production vertex of a daughter
            dvtx = ev[daughters[0]]
            Ldecay = np.hypot(np.hypot(dvtx.xProd() - p.xProd(),
                                       dvtx.yProd() - p.yProd()),
                              dvtx.zProd() - p.zProd())       # mm
            charm_idx[i] = dict(
                sign=1 if p.id() > 0 else -1,
                pid=p.id(),
                E=p.e(), p=p.pAbs(),
                Ldecay=Ldecay,
                has_mu=False, mu_q=np.nan, mu_p=np.nan,
                mu_theta=np.nan, mu_ptrel=np.nan,
                w=w,
            )

        # 2) Attach any semileptonic muon to its charm ancestor.  Muons from
        #    charm decay have a charmed-hadron ancestor; the hard scattered muon
        #    does not (its chain goes to the beam), so ancestry cleanly separates.
        for i in range(ev.size()):
            p = ev[i]
            if abs(p.id()) != 13:
                continue
            a = charm_ancestor(ev, i)
            if a < 0 or a not in charm_idx:
                continue
            rec = charm_idx[a]
            if rec["has_mu"]:
                continue                       # keep the first (leading) mu
            ch = ev[a]
            # pT of the muon relative to the parent charm-hadron direction
            pvec = np.array([p.px(), p.py(), p.pz()])
            cvec = np.array([ch.px(), ch.py(), ch.pz()])
            cnorm = np.linalg.norm(cvec)
            if cnorm > 0:
                plong = np.dot(pvec, cvec) / cnorm
                ptrel = np.sqrt(max(0.0, np.dot(pvec, pvec) - plong**2))
            else:
                ptrel = np.nan
            rec.update(
                has_mu=True,
                mu_q=1 if p.id() < 0 else -1,   # mu+ (id -13) -> +1
                mu_p=p.pAbs(),
                mu_theta=np.arctan2(np.hypot(p.px(), p.py()), p.pz()),
                mu_ptrel=ptrel,
            )

        records.extend(charm_idx.values())
        n += 1

    pythia.stat()
    tag = Path(lhe).parent.parent.name + "/" + Path(lhe).parent.parent.parent.name \
          if "pythia8" in lhe else Path(lhe).name
    print(f"  {tag}: {n} events showered, {len(records)} charm hadrons")
    return records


def load_group(group, n_events):
    """Shower all four samples, tungsten-scale the weights, concatenate arrays."""
    keys = ["sign", "pid", "E", "p", "Ldecay",
            "has_mu", "mu_q", "mu_p", "mu_theta", "mu_ptrel", "w"]
    cols = {k: [] for k in keys}
    n_events_total = 0
    for name, lhe in group:
        if not Path(lhe).exists():
            print(f"  [warn] missing {lhe}", file=sys.stderr)
            continue
        recs = shower_one(lhe, n_events)
        if not recs:
            continue
        tf = tung_factor(name)
        for r in recs:
            for k in keys:
                cols[k].append(tf * r[k] if k == "w" else r[k])
        n_ch = len(recs)
        n_mu = sum(r["has_mu"] for r in recs)
        print(f"    {name}: tung={tf:.3f}, charm={n_ch}, "
              f"semileptonic-mu={n_mu} ({100*n_mu/max(n_ch,1):.1f}%)")
    arr = {k: np.array(cols[k], dtype=float) for k in keys}
    return arr


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nevents", type=int, default=20000,
                    help="Max events to shower per sample (default 20000)")
    ap.add_argument("--prod", default="production_v1",
                    help="production_v1-layout run to shower (default: production_v1)")
    ap.add_argument("--npz", default=config.get("charm_npz"),
                    help="Output cache path [default: json settings charm_npz]")
    args = ap.parse_args()

    print(f"Showering {args.prod} ...")
    arr = load_group(signal_group(args.prod), args.nevents)
    n_mu = int(arr["has_mu"].sum()) if len(arr["has_mu"]) else 0
    n_ch = len(arr["sign"])
    print(f"\nTotal: {n_ch} charm hadrons, {n_mu} with a decay muon "
          f"({100*n_mu/max(n_ch,1):.1f}% semileptonic-to-muon)")
    nc  = int((arr["sign"] > 0).sum())
    nca = int((arr["sign"] < 0).sum())
    print(f"Truth split: c = {nc}, c~ = {nca}  (raw counts, unweighted)")

    Path(args.npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.npz, **arr)
    print(f"Cached -> {args.npz}")


if __name__ == "__main__":
    main()
