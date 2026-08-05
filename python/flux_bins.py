"""
Energy binning and flux weights for the BINNED generation strategy.

WHY BINNED GENERATION.  The original samples were produced by handing the muon
flux to POWHEG as a "beam PDF" of a fictitious 7 TeV muon beam
(arXiv:2506.13889 Eq. 2.1).  That is elegant but has two costs:

  * the beam energy is sampled inside POWHEG, so the true E_mu of an event is
    only available through the hard-process record;
  * the fictitious beam has a REMNANT, which hadronises and deposits a few GeV
    of activity that is not part of the DIS vertex.  That contamination broke
    the hadronic-energy sum (energy conservation failed in 80% of events, and
    E_had/nu reached 2.5 at large x), and could not be removed by any
    particle-level cut -- see docs/XRECO.md.  It forced the hadronic four-vector
    to be taken from momentum conservation instead of from the particles.

Generating in fixed-energy bins removes both problems: E_mu is exact and equal to
the bin energy, and with `fixed_lepton_beam 1` the lepton carries the full beam
momentum so there is no lepton-side remnant at all.  The flux is then applied
afterwards as a per-bin weight.

THE REWEIGHTING.  In the beam-PDF approach POWHEG returns
    sigma_tot = Int dx f(x) sigma_hat(x E_beam)
and because f already carries the target factor n_T L_T of Eq. (2.1), sigma_tot
is directly an event count.  Binning that integral gives, for bin i,

    N_i  =  [ Int_{bin i} f(x) dx ]  x  sigma_hat(E_i)

with sigma_hat(E_i) the cross section POWHEG returns for a monochromatic beam at
the bin's energy.  Summing over bins reproduces the original normalisation by
construction, so the absolute scale is unchanged -- this is a change of
generation strategy, not of physics.

The LHAPDF grid stores x*f(x), so f(x) = grid/x.
"""
import numpy as np

# numpy 2 renamed trapz -> trapezoid; keep working on both
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz

# The flux grid used for the original production
FLUX_GRID = ("/afs/cern.ch/work/e/evilla/private/faser/lhapdf_data/"
             "muon_flux_FASERv_Run3_var2/muon_flux_FASERv_Run3_var2_0000.dat")
E_BEAM_REF = 7000.0        # GeV, the fictitious beam the flux x refers to

N_BINS = 20                # 20 bins in log E, as requested
E_MIN, E_MAX = 10.0, 3000.0   # GeV; the flux is zero above ~3.1 TeV


def read_flux(path=FLUX_GRID):
    """Return (x, f_mu(x), f_mubar(x)) from the LHAPDF flux grid."""
    blocks = [b for b in open(path).read().split("---") if b.strip()]
    lines = [l for l in blocks[1].strip().split("\n") if l.strip()]
    xs = np.array([float(v) for v in lines[0].split()])
    qs = np.array([float(v) for v in lines[1].split()])
    flav = lines[2].split()
    data = np.array([[float(v) for v in r.split()] for r in lines[3:] if r.strip()])
    data = data.reshape(len(xs), len(qs), len(flav))
    # grid holds x*f(x); take the lowest Q node (the flux has no Q dependence)
    xf_mu = data[:, 0, flav.index("13")]
    xf_mubar = data[:, 0, flav.index("-13")]
    return xs, xf_mu, xf_mubar


def bin_edges(n=N_BINS, e_min=E_MIN, e_max=E_MAX):
    """Bin edges, equally spaced in log E."""
    return np.logspace(np.log10(e_min), np.log10(e_max), n + 1)


def flux_weights(n=N_BINS, e_min=E_MIN, e_max=E_MAX, path=FLUX_GRID):
    """
    Per-bin flux integrals Int f(x) dx, separately for mu- and mu+.

    Returns a dict with the bin edges, centres (log-mid, used as the generation
    energy), and the two integrals.  Multiply these by the POWHEG cross section
    at the bin energy to get the expected event count.
    """
    xs, xf_mu, xf_mubar = read_flux(path)
    edges = bin_edges(n, e_min, e_max)
    # integrate f(x) = (x f)/x on a fine log grid inside each bin
    out_mu, out_mubar = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        xg = np.logspace(np.log10(lo / E_BEAM_REF), np.log10(hi / E_BEAM_REF), 200)
        f_mu = np.interp(xg, xs, xf_mu) / xg
        f_mb = np.interp(xg, xs, xf_mubar) / xg
        out_mu.append(np.trapezoid(f_mu, xg))
        out_mubar.append(np.trapezoid(f_mb, xg))
    return dict(edges=edges,
                e_centre=np.sqrt(edges[:-1] * edges[1:]),
                w_mu=np.array(out_mu), w_mubar=np.array(out_mubar))


if __name__ == "__main__":
    b = flux_weights()
    print(f"{'bin':>4} {'E_lo':>9} {'E_hi':>9} {'E_gen':>9} "
          f"{'int f dx (mu-)':>16} {'int f dx (mu+)':>16}")
    for i, (lo, hi, ec, wm, wp) in enumerate(zip(b["edges"][:-1], b["edges"][1:],
                                                 b["e_centre"], b["w_mu"],
                                                 b["w_mubar"])):
        print(f"{i:4d} {lo:9.1f} {hi:9.1f} {ec:9.1f} {wm:16.4e} {wp:16.4e}")
    print(f"\ntotal   mu- {b['w_mu'].sum():.4e}   mu+ {b['w_mubar'].sum():.4e}")
