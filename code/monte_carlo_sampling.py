"""
Monte Carlo Parameter Scan
==========================

Bayesian MCMC parameter-space exploration for UHDM production in AGN jets
using the emcee affine-invariant ensemble sampler (Foreman-Mackey et al. 2013).

Scanned parameters (all in log10):
    log10(Lambda / GeV)   -- EFT cutoff scale
    log10(m_X    / GeV)   -- UHDM mass
    log10(sigma_jet)      -- jet magnetization
    log10(B_jet  / G)     -- jet magnetic field

Sampler configuration:
    ndim     = 4
    nwalkers = 10   (emcee requires nwalkers >= 2*ndim; 2*4 = 8 minimum)
    nsteps   = 100
    Total    = 10 * 100 = 1000 samples

Random seeds for reproducibility: {42, 137, 271, 314, 628}

Convergence diagnostics:
    - Acceptance fraction
    - Autocorrelation time (emcee.autocorr.integrated_time)
    - Split-R-hat (Gelman-Rubin) computed from walker sub-chains
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import emcee
import corner

# Ensure sibling modules are importable
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from eft_cross_sections import eft_dimension6_xsec
from grmhd_acceleration import gamma_electron_reconnection, gamma_proton_hillas
from observational_constraints import MultiMessengerConstraints, calculate_uhdm_flux
from pdf_integration import CT18Luminosity, hadronic_cross_section

FIGURES_DIR = os.path.join(SCRIPT_DIR, '..', 'figures')

# ---------------------------------------------------------------------------
# Parameter space
# ---------------------------------------------------------------------------
PARAM_NAMES  = [r'$\log_{10}\Lambda$', r'$\log_{10}m_X$',
                r'$\log_{10}\sigma$', r'$\log_{10}B$']
PARAM_LABELS = ['log10(Lambda/GeV)', 'log10(m_X/GeV)',
                'log10(sigma_jet)', 'log10(B_jet/G)']
NDIM = 4

# Log-uniform prior bounds
PRIOR_LO = np.array([8.0, 6.0, 1.0, -1.0])   # Lambda>=1e8, m_X>=1e6, sigma>=10, B>=0.1
PRIOR_HI = np.array([11.0, 10.0, 3.0, 2.0])   # Lambda<=1e11, m_X<=1e10, sigma<=1e3, B<=100


# ---------------------------------------------------------------------------
# Prior, likelihood, posterior
# ---------------------------------------------------------------------------
def log_prior(theta):
    """Flat (log-uniform) prior on all parameters."""
    if np.all((theta >= PRIOR_LO) & (theta <= PRIOR_HI)):
        return 0.0
    return -np.inf


def log_likelihood(theta):
    """
    Log-likelihood from multi-messenger constraints.

    For each observable channel (Fermi, IceCube, Auger) the likelihood
    penalises parameter combinations whose predicted flux exceeds the
    experimental upper limit:

        log L = -0.5 * sum_i [ max(0, F_pred_i - F_lim_i)^2 / sigma_i^2 ]

    This is a one-sided Gaussian penalty (upper-limit likelihood).
    """
    log_Lambda, log_mX, log_sigma, log_B = theta

    Lambda    = 10.0**log_Lambda   # GeV
    m_X       = 10.0**log_mX       # GeV
    sigma_jet = 10.0**log_sigma
    B_jet     = 10.0**log_B        # G

    # -- Particle physics: EFT cross section at AGN energy -------------------
    sqrt_s = 3.16e8   # GeV  (geometric mean of 10^17--10^18 eV)
    s_val  = sqrt_s**2

    xsec = eft_dimension6_xsec(s_val, m_X, Lambda)
    if not np.isfinite(xsec) or xsec <= 0:
        return -1e10   # strongly penalise invalid regions

    # -- Astrophysics: production rate in the jet ----------------------------
    # Lorentz factors (sanity check that sigma supports the energy)
    gamma_e = gamma_electron_reconnection(sigma_jet, eta=0.1)
    R_jet   = 1e15    # cm  (fiducial)
    gamma_p = gamma_proton_hillas(B_jet, R_jet)

    # Effective proton-proton collision rate  N_dot ~ n_jet * sigma_pp * c * V_jet
    # Use simplified estimate: V_jet ~ pi R_jet^2 * R_jet, n ~ 10^3 cm^-3
    n_jet = 1e3  # cm^-3
    V_jet = np.pi * R_jet**3
    c_cgs = 2.998e10  # cm/s

    # Convolve partonic xsec with parton luminosity
    lumi = CT18Luminosity()
    def sigma_hat(s_hat, mx, L=Lambda):
        return eft_dimension6_xsec(s_hat, mx, L)
    sigma_pp = hadronic_cross_section(sigma_hat, s_val, m_X, lumi)
    if sigma_pp <= 0:
        return -1e10

    N_dot = n_jet**2 * sigma_pp * c_cgs * V_jet   # pairs / s

    if N_dot <= 0:
        return -1e10

    # -- Multi-messenger fluxes at Earth (d = 100 Mpc fiducial) ---------------
    fluxes = calculate_uhdm_flux(m_X, N_dot, distance_Mpc=100.0)

    constraints = MultiMessengerConstraints()

    # One-sided Gaussian penalty for each channel
    log_L = 0.0

    # Fermi-LAT
    excess_fermi = max(0.0, fluxes['E2dNdE_gamma'] - constraints.fermi_E2dNdE_limit)
    sigma_fermi  = 0.5 * constraints.fermi_E2dNdE_limit   # 50% uncertainty
    log_L -= 0.5 * (excess_fermi / sigma_fermi)**2

    # IceCube
    excess_ic = max(0.0, fluxes['phi_nu'] - constraints.icecube_phi_limit)
    sigma_ic  = 0.5 * constraints.icecube_phi_limit
    log_L -= 0.5 * (excess_ic / sigma_ic)**2

    # Auger  (very weak constraint at these energies)
    excess_au = max(0.0, fluxes['integral_CR'] - constraints.auger_integral_flux_limit)
    sigma_au  = 0.5 * constraints.auger_integral_flux_limit
    if sigma_au > 0:
        log_L -= 0.5 * (excess_au / sigma_au)**2

    return log_L


def log_probability(theta):
    """Log-posterior = log-prior + log-likelihood."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


# ---------------------------------------------------------------------------
# Convergence diagnostics
# ---------------------------------------------------------------------------
def compute_split_rhat(chain):
    """
    Split-R-hat (Gelman-Rubin) from a single ensemble chain.

    Treats each walker as an independent chain, splits each in half,
    then computes the standard R-hat across all half-chains.

    Parameters
    ----------
    chain : ndarray, shape (nsteps, nwalkers, ndim)

    Returns
    -------
    rhat : ndarray, shape (ndim,)
    """
    nsteps, nwalkers, ndim = chain.shape
    half = nsteps // 2
    # Take the second half (after burn-in) and split again
    second_half = chain[half:]
    n = second_half.shape[0] // 2
    if n < 2:
        return np.full(ndim, np.nan)

    # Create sub-chains: shape (2*nwalkers, n, ndim)
    chains = np.concatenate([second_half[:n], second_half[n:2*n]], axis=1)
    chains = np.transpose(chains, (1, 0, 2))  # (2*nwalkers, n, ndim)

    m = chains.shape[0]   # number of sub-chains
    n = chains.shape[1]   # length of each

    chain_means = chains.mean(axis=1)            # (m, ndim)
    chain_vars  = chains.var(axis=1, ddof=1)     # (m, ndim)

    overall_mean = chain_means.mean(axis=0)      # (ndim,)
    B = n * np.var(chain_means, axis=0, ddof=1)  # between-chain var
    W = np.mean(chain_vars, axis=0)              # within-chain var

    var_hat = (1.0 - 1.0/n) * W + (1.0/n) * B
    rhat = np.sqrt(var_hat / np.maximum(W, 1e-30))

    return rhat


# ---------------------------------------------------------------------------
# Main MCMC runner
# ---------------------------------------------------------------------------
def run_mcmc(nwalkers=10, nsteps=100, seed=42):
    """
    Run the MCMC sampler.

    Parameters
    ----------
    nwalkers : int
        Number of walkers (must be >= 2 * ndim = 8).
    nsteps : int
        Number of steps per walker.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    sampler : emcee.EnsembleSampler
    """
    assert nwalkers >= 2 * NDIM, (
        f"emcee requires nwalkers >= 2*ndim = {2*NDIM}, got {nwalkers}"
    )

    rng = np.random.default_rng(seed)

    # Initialise walkers: uniform within prior volume
    p0 = rng.uniform(PRIOR_LO, PRIOR_HI, size=(nwalkers, NDIM))

    print(f"\nMCMC configuration:")
    print(f"  ndim     = {NDIM}")
    print(f"  nwalkers = {nwalkers}")
    print(f"  nsteps   = {nsteps}")
    print(f"  seed     = {seed}")
    print(f"  total samples = {nwalkers * nsteps}")

    sampler = emcee.EnsembleSampler(nwalkers, NDIM, log_probability)

    print("\nRunning MCMC ...")
    sampler.run_mcmc(p0, nsteps, progress=False)
    print("  ... done.")

    return sampler


def print_diagnostics(sampler):
    """Print convergence diagnostics."""
    chain = sampler.get_chain()   # shape (nsteps, nwalkers, ndim)
    nsteps, nwalkers, ndim = chain.shape

    # Acceptance fraction
    af = sampler.acceptance_fraction
    print(f"\nAcceptance fractions: mean={af.mean():.3f}, "
          f"min={af.min():.3f}, max={af.max():.3f}")

    # Autocorrelation time (may fail with too few samples)
    try:
        tau = emcee.autocorr.integrated_time(
            sampler.get_chain(), quiet=True)
        print(f"Autocorrelation time: {tau}")
        print(f"Effective sample size: "
              f"{nsteps * nwalkers / np.max(tau):.0f}")
    except emcee.autocorr.AutocorrError:
        print("Autocorrelation time: could not estimate "
              "(chain too short; increase nsteps for production runs)")

    # Split R-hat
    rhat = compute_split_rhat(chain)
    print(f"Split R-hat: {rhat}")
    for name, r in zip(PARAM_LABELS, rhat):
        status = 'OK' if (np.isfinite(r) and r < 1.1) else 'WARN'
        print(f"  {name:25s}  R-hat = {r:.3f}  [{status}]")


def plot_corner(sampler, burn_in=20, output_file=None):
    """
    Generate corner plot (Figure S1) from MCMC samples.

    Parameters
    ----------
    sampler : emcee.EnsembleSampler
    burn_in : int
        Number of steps to discard as burn-in.
    output_file : str or None
    """
    if output_file is None:
        output_file = os.path.join(FIGURES_DIR, 'figS1_mcmc_diagnostics.pdf')

    flat_samples = sampler.get_chain(discard=burn_in, flat=True)
    print(f"\nCorner plot: {flat_samples.shape[0]} samples "
          f"(discarded {burn_in} burn-in steps)")

    fig = corner.corner(
        flat_samples,
        labels=PARAM_NAMES,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={'fontsize': 11},
        smooth=1.5,
        bins=20,
    )
    fig.suptitle('MCMC Parameter Scan — Convergence Diagnostics (Figure S1)',
                 fontsize=13, y=1.02)

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    fig.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Figure saved: {output_file}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='MCMC parameter scan for UHDM production in AGN jets')
    parser.add_argument('--nwalkers', type=int, default=10,
                        help='Number of walkers (default: 10, min: 8)')
    parser.add_argument('--nsteps', type=int, default=100,
                        help='Steps per walker (default: 100)')
    parser.add_argument('--seed', type=int, default=42,
                        choices=[42, 137, 271, 314, 628],
                        help='Random seed (default: 42)')
    parser.add_argument('--burn-in', type=int, default=20,
                        help='Burn-in steps to discard (default: 20)')
    args = parser.parse_args()

    print("Monte Carlo Parameter Scan")
    print("=" * 55)

    sampler = run_mcmc(nwalkers=args.nwalkers,
                       nsteps=args.nsteps,
                       seed=args.seed)

    print_diagnostics(sampler)
    plot_corner(sampler, burn_in=args.burn_in)
    print("\nDone.")


if __name__ == '__main__':
    main()
