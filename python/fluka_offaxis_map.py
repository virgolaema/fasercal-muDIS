#!/usr/bin/env python3
"""
Build the transverse muon-fluence map from the FASER FLUKA simulation, and from
it the OFF-AXIS FLUX FACTOR F_flux needed by fasercal_chain.py.

Source: /eos/experiment/fasernu-data0/faser/sim/mc22/fluka/210007/bck/s0010-r0019/
which is the same input the existing muon_fluxes/read_fluka_sim.cpp converter
uses to build the LHAPDF flux grids.  Conventions taken from that converter:

  * positions truth_prod_x/y are in mm, with the ORIGIN on the FASERnu axis
    (the converter's FASERnu window is |x|<140 mm, |y|<200 mm about 0,0);
  * momenta truth_px/py/pz are in MeV;
  * only the FIRST truth particle of each (de-duplicated) event is used;
  * weight = CrossSection * 2.0 / 1.3713  normalises the flux to fb^-1.

F_flux is a RATIO of areal fluence (off-axis window / on-axis FASERnu window),
so the absolute normalisation cancels -- this is exactly the factor the chain
needs, and it is independent of the residual absolute-normalisation puzzle.
"""
import sys, glob
import numpy as np
import ROOT

SRC = "/eos/experiment/fasernu-data0/faser/sim/mc22/fluka/210007/bck/s0010-r0019/"
WIN_X, WIN_Y = 140.0, 200.0          # mm, on-axis FASERnu window (half-widths)
WNORM = 2.0 / 1.3713                 # per-fb^-1 normalisation from the converter


def read(files, max_entries_per_file=None):
    xs, ys, ps, ws, qs = [], [], [], [], []
    for fn in files:
        f = ROOT.TFile.Open(fn)
        if not f or f.IsZombie():
            print(f"  [warn] cannot open {fn}", file=sys.stderr); continue
        t = f.Get("nt")
        n = t.GetEntries()
        if max_entries_per_file:
            n = min(n, max_entries_per_file)
        print(f"  {fn.split('/')[-1]}: {n} entries", flush=True)
        prev = -1
        for i in range(n):
            t.GetEntry(i)
            if t.truth_pdg.size() < 1:
                continue
            pdg = t.truth_pdg[0]
            if abs(pdg) != 13:
                continue
            ev = t.eventID
            if ev == prev:
                continue                     # de-duplicate, as the converter does
            prev = ev
            px, py, pz = t.truth_px[0], t.truth_py[0], t.truth_pz[0]
            xs.append(t.truth_prod_x[0]); ys.append(t.truth_prod_y[0])
            ps.append(np.sqrt(px*px + py*py + pz*pz) * 1e-3)   # GeV
            ws.append(t.CrossSection * WNORM)
            qs.append(1 if pdg < 0 else -1)                    # mu+ -> +1
        f.Close()
    return (np.array(xs), np.array(ys), np.array(ps), np.array(ws), np.array(qs))


def main():
    nfiles = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    maxent = int(sys.argv[2]) if len(sys.argv) > 2 else 400000
    files = sorted(glob.glob(SRC + "*.root"))[:nfiles]
    print(f"Reading {len(files)} FLUKA files (<= {maxent} entries each)...")
    x, y, p, w, q = read(files, maxent)
    print(f"\nmuons: {len(x)}")
    if not len(x):
        return

    onax = (np.abs(x) < WIN_X) & (np.abs(y) < WIN_Y)
    area_on = (2 * WIN_X) * (2 * WIN_Y)          # mm^2
    dens_on = w[onax].sum() / area_on
    print(f"on-axis window |x|<{WIN_X:.0f}, |y|<{WIN_Y:.0f} mm: "
          f"{onax.sum()} muons, areal density {dens_on:.4g} /fb/mm^2")
    print(f"  mu+ / mu- in window: {(q[onax]>0).sum()} / {(q[onax]<0).sum()}")
    print(f"  <p> = {p[onax].mean():.0f} GeV, median {np.median(p[onax]):.0f} GeV")

    # scan a same-sized window across the transverse plane
    print("\n=== F_flux for a same-size window centred at (x0,y0) ===")
    print(f"{'x0[cm]':>7} {'y0[cm]':>7} {'N':>7} {'F_flux':>8} {'<p>[GeV]':>9} {'mu+frac':>8}")
    for x0 in (0, 25, 50, 100, 150, -25, -50, -100, -150):
        for y0 in (0, 50, 100):
            m = (np.abs(x - x0*10) < WIN_X) & (np.abs(y - y0*10) < WIN_Y)
            if m.sum() < 20:
                continue
            F = (w[m].sum() / area_on) / dens_on
            print(f"{x0:7.0f} {y0:7.0f} {m.sum():7d} {F:8.3f} "
                  f"{p[m].mean():9.0f} {(q[m]>0).mean():8.2f}")

    np.savez("/eos/home-e/evilla/faser/fluka_muon_map.npz",
             x=x, y=y, p=p, w=w, q=q, win_x=WIN_X, win_y=WIN_Y)
    print("\nCached -> /eos/home-e/evilla/faser/fluka_muon_map.npz")


if __name__ == "__main__":
    main()
