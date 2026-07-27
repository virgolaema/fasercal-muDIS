"""
Toy FASERcal + downstream-spectrometer response for charm/anticharm tagging.

FASERcal is a 3D-granular plastic scintillator: NO magnetic field, so it cannot
measure the sign of a charge on its own.  The charm/anticharm tag therefore has
two independent ingredients, both modelled here as simple, explicit parametric
functions so the assumptions are easy to see and vary:

  1. ACCEPTANCE of the semileptonic decay muon into a downstream magnetised
     spectrometer (FASER2-like): the muon must go forward, above an angular
     acceptance, and above a momentum threshold to be measured / bent.

  2. CHARGE-CONFUSION of that spectrometer vs muon momentum: the higher the
     momentum, the straighter the track, the more likely the sign is wrong.

None of these numbers are FASERcal specifications — they are deliberately
transparent toy defaults (see DEFAULTS) meant to be scanned.  The point of the
study is to see how the c-vs-c~ tagging efficiency and purity depend on them.
"""
import numpy as np

# Plastic-scintillator minimum-ionising energy loss, used to turn a muon
# momentum into an ionisation range (range ~ p / dEdx).  Polystyrene ~2 MeV/cm.
DEDX_MIP_MEV_PER_CM = 2.0

DEFAULTS = dict(
    theta_max_mrad=10.0,   # angular acceptance half-cone of the spectrometer [mrad]
    p_min_gev=5.0,         # muon momentum threshold to be measured [GeV]
    # charge-confusion model: eta(p) = eta0 + (0.5-eta0)*sigmoid((p-p_half)/w)
    eta0=0.005,            # floor charge-confusion at low p
    p_half_gev=1000.0,     # momentum at which confusion reaches half-way to 0.5
    p_width_gev=300.0,     # sigmoid width [GeV]
)


def muon_range_m(p_gev):
    """Ionisation (CSDA-like) range of a muon in plastic scintillator [m]."""
    return (p_gev * 1e3) / DEDX_MIP_MEV_PER_CM / 100.0   # MeV / (MeV/cm) -> cm -> m


def charge_confusion(p_gev, eta0, p_half_gev, p_width_gev):
    """Probability the spectrometer assigns the wrong charge sign at momentum p."""
    sig = 1.0 / (1.0 + np.exp(-(p_gev - p_half_gev) / p_width_gev))
    return eta0 + (0.5 - eta0) * sig


def accepts(mu_p, mu_theta, theta_max_mrad, p_min_gev, **_):
    """Boolean: does the decay muon reach and get measured by the spectrometer?"""
    forward = mu_theta < (theta_max_mrad * 1e-3)
    return forward & (mu_p > p_min_gev)


def tag_charge(true_q, mu_p, rng, eta0, p_half_gev, p_width_gev, **_):
    """
    Reconstructed charge sign: the true sign flipped with probability
    charge_confusion(p).  Returns an int array of +-1.
    """
    eta = charge_confusion(mu_p, eta0, p_half_gev, p_width_gev)
    flip = rng.random(len(mu_p)) < eta
    return np.where(flip, -true_q, true_q)
