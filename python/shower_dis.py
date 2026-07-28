#!/usr/bin/env python3
"""
Per-EVENT muon-DIS extraction for the FASERcal charm chain.

Unlike shower_charm.py (one row per charm hadron), this writes one row per DIS
EVENT, which is what the yield chain needs:

    N(muon DIS in FASERcal)
      -> N(events with charm)
        -> N(charm events with a semileptonic muon)
          -> N(those muons reaching the spectrometer)

Per event we record:

  DIS kinematics (from the Pythia HARD PROCESS record = LHE truth):
    w_raw     : nominal POWHEG event weight, NOT target-weighted (see below)
    is_neutron: 1 if the target nucleon was a neutron, else 0
    e_in      : incoming muon energy [GeV]
    p_out     : scattered muon |p| [GeV]
    theta_mu  : scattered-muon polar angle w.r.t. the incoming muon [rad]
    q2        : truth Q^2 [GeV^2]

  Hadronic final state (for the x-reconstruction methods, JB / Sigma / DA):
    e_had     : visible hadronic energy, neutrinos and the scattered muon removed
    px/py/pz_had : hadronic momentum vector [GeV]
    sigma_had : Sum_h (E - pz), the Sigma-method variable

  Charm chain:
    n_charm      : number of weakly-decaying charm hadrons
    n_semilep_mu : how many of them decayed to a muon
    mu2_q/p/theta: charge, |p| [GeV] and polar angle of the LEADING decay muon
    charm_sign_mu: sign of the charm hadron that produced that muon (+1 c, -1 c~)

TARGET COMPOSITION IS DELIBERATELY NOT APPLIED HERE.  The four samples are
per-nucleon (mu-+p, mu-+n, mu++p, mu++n), so storing the raw weight plus an
is_neutron flag lets the analysis recombine for ANY target: tungsten
(w_p=74/184) for the FASERnu cross-check, or scintillator CH (w_p=56/104) for
FASERcal.  Baking tungsten in here, as the FASERnu study did, would be wrong
for a scintillator detector.

Run inside an LCG view providing pythia8, e.g.
  source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
  python3 shower_dis.py --prod production_v1 --seeds 1-10
"""
import argparse
import sys
from pathlib import Path

import numpy as np

import pythia8
import config

PROD    = config.get("production_dir")
SAMPLES = ["mum_proton", "mum_neutron", "mup_proton", "mup_neutron"]


def signal_group(prod_name, seed=1):
    return [(s, f"{PROD}/{prod_name}/{s}/pythia8_seed{seed}/pwgevents.lhe")
            for s in SAMPLES]


def parse_seeds(spec):
    out = set()
    for part in str(spec).split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a), int(b) + 1))
        elif part:
            out.add(int(part))
    return sorted(out)


def is_charm_hadron(pid):
    a = abs(pid)
    if a < 100:
        return False
    return (a // 100) % 10 == 4 or (a // 1000) % 10 == 4


def charm_ancestor(event, i):
    """Index of the nearest charmed-hadron ancestor of particle i, or -1."""
    seen, stack = set(), [i]
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


KEYS = ["w_raw", "is_neutron", "e_in", "p_out", "theta_mu", "q2",
        "e_had", "px_had", "py_had", "pz_had", "sigma_had",
        "n_charm", "n_semilep_mu", "mu2_q", "mu2_p", "mu2_theta",
        "charm_sign_mu"]


def shower_one(lhe, n_events, is_neutron):
    pythia = pythia8.Pythia("", False)
    for cmd in [
        "Beams:frameType = 4", f"Beams:LHEF = {lhe}",
        "Main:numberOfEvents = 0", "Print:quiet = on",
        "Init:showChangedSettings = off", "Init:showChangedParticleData = off",
        "Next:numberShowInfo = 0", "Next:numberShowProcess = 0",
        "Next:numberShowEvent = 0", "Stat:showProcessLevel = off",
        "POWHEG:nFinal = 2",
    ]:
        pythia.readString(cmd)
    if not pythia.init():
        print(f"  [warn] Pythia init failed for {lhe}", file=sys.stderr)
        return {k: [] for k in KEYS}

    cols = {k: [] for k in KEYS}
    n = 0
    while n < n_events:
        if not pythia.next():
            if pythia.infoPython().atEndOfFile():
                break
            continue

        # ---- DIS kinematics from the hard process (LHE truth) ----
        proc = pythia.process
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
        q2 = -((in_mu.e() - out_mu.e())**2
               - (in_mu.px() - out_mu.px())**2
               - (in_mu.py() - out_mu.py())**2
               - (in_mu.pz() - out_mu.pz())**2)
        # scattering angle between the incoming and outgoing muon directions
        vin = np.array([in_mu.px(), in_mu.py(), in_mu.pz()])
        vout = np.array([out_mu.px(), out_mu.py(), out_mu.pz()])
        nin, nout = np.linalg.norm(vin), np.linalg.norm(vout)
        theta_mu = np.arccos(np.clip(np.dot(vin, vout) / (nin * nout), -1, 1)) \
            if nin > 0 and nout > 0 else np.nan

        ev = pythia.event

        # ---- charm hadrons: keep only the weakly-decaying ones ----
        charm_idx = {}
        for i in range(ev.size()):
            p = ev[i]
            if not is_charm_hadron(p.id()):
                continue
            daughters = [d for d in p.daughterList() if d > 0]
            if not daughters or any(is_charm_hadron(ev[d].id()) for d in daughters):
                continue
            charm_idx[i] = {"sign": 1 if p.id() > 0 else -1, "mu": None}

        # ---- muons: separate the scattered muon from charm-decay muons ----
        # A charm-decay muon has a charmed-hadron ancestor; the scattered muon
        # does not (its chain leads back to the beam).  This is what makes the
        # dimuon charm tag possible.
        scattered_idx, best_p = -1, -1.0
        for i in range(ev.size()):
            p = ev[i]
            if not p.isFinal() or abs(p.id()) != 13:
                continue
            a = charm_ancestor(ev, i)
            if a >= 0 and a in charm_idx:
                rec = charm_idx[a]
                if rec["mu"] is None or p.pAbs() > rec["mu"].pAbs():
                    rec["mu"] = p
            else:
                if p.pAbs() > best_p:            # leading non-charm muon
                    best_p, scattered_idx = p.pAbs(), i

        # ---- hadronic final state: everything except neutrinos and the
        #      scattered muon.  The muon-beam remnant photon (status 62, id 22)
        #      MUST be dropped: the flux is used as a beam PDF of a fictitious
        #      7 TeV beam, so Pythia parks the unused balance in a ~5 TeV
        #      photon that is not physical.
        e_had = px = py = pz = sigma = 0.0
        for i in range(ev.size()):
            p = ev[i]
            if not p.isFinal() or i == scattered_idx:
                continue
            if abs(p.id()) in (12, 14, 16):
                continue
            if p.statusAbs() == 62 and p.id() == 22:
                continue
            e_had += p.e(); px += p.px(); py += p.py(); pz += p.pz()
            sigma += p.e() - p.pz()

        # ---- leading semileptonic decay muon ----
        n_semi = sum(1 for r in charm_idx.values() if r["mu"] is not None)
        lead_mu, lead_sign = None, np.nan
        for r in charm_idx.values():
            if r["mu"] is None:
                continue
            if lead_mu is None or r["mu"].pAbs() > lead_mu.pAbs():
                lead_mu, lead_sign = r["mu"], r["sign"]

        row = dict(
            w_raw=pythia.infoPython().weight(), is_neutron=float(is_neutron),
            e_in=in_mu.e(), p_out=out_mu.pAbs(), theta_mu=theta_mu, q2=q2,
            e_had=e_had, px_had=px, py_had=py, pz_had=pz, sigma_had=sigma,
            n_charm=float(len(charm_idx)), n_semilep_mu=float(n_semi),
            mu2_q=(1.0 if lead_mu is not None and lead_mu.id() < 0 else
                   (-1.0 if lead_mu is not None else np.nan)),
            mu2_p=lead_mu.pAbs() if lead_mu is not None else np.nan,
            mu2_theta=(np.arctan2(np.hypot(lead_mu.px(), lead_mu.py()),
                                  lead_mu.pz()) if lead_mu is not None else np.nan),
            charm_sign_mu=float(lead_sign),
        )
        for k in KEYS:
            cols[k].append(row[k])
        n += 1

    pythia.stat()
    print(f"    {Path(lhe).parent.parent.name}: {n} events")
    return cols


def load_group(group, n_events):
    cols = {k: [] for k in KEYS}
    for name, lhe in group:
        if not Path(lhe).exists():
            print(f"  [warn] missing {lhe}", file=sys.stderr)
            continue
        c = shower_one(lhe, n_events, is_neutron=("neutron" in name))
        if not c["w_raw"]:
            continue
        # normalise the weight by the number of generated events in THIS sample
        # so that summing over samples gives the per-nucleon expectation
        nev = len(c["w_raw"])
        c["w_raw"] = [w / nev for w in c["w_raw"]]
        for k in KEYS:
            cols[k].extend(c[k])
    return {k: np.array(cols[k], dtype=float) for k in KEYS}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nevents", type=int, default=20000)
    ap.add_argument("--prod", default="production_v1")
    ap.add_argument("--seeds", default="1-10")
    ap.add_argument("--npz", required=True)
    args = ap.parse_args()

    seeds = parse_seeds(args.seeds)
    print(f"Showering {args.prod}, seeds {seeds}")
    parts = []
    for sd in seeds:
        print(f"--- seed {sd} ---")
        parts.append(load_group(signal_group(args.prod, sd), args.nevents))
    arr = {k: np.concatenate([p[k] for p in parts]) for k in KEYS}
    # each seed is an independent sample of the same luminosity -> average
    arr["w_raw"] = arr["w_raw"] / len(seeds)

    n = len(arr["w_raw"])
    print(f"\nTotal {n} events; "
          f"charm events {int((arr['n_charm'] >= 1).sum())}, "
          f"with semileptonic mu {int((arr['n_semilep_mu'] >= 1).sum())}")
    Path(args.npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.npz, **arr)
    print(f"Cached -> {args.npz}")


if __name__ == "__main__":
    main()
