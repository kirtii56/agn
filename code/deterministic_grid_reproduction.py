"""
Deterministic Grid Reproduction
================================

Reproduces the reported deterministic (log10 Lambda, log10 m_X) inference
result: the 221-point grid characterization of the accepted/rejected region
and the closed-form accepted-region boundaries.

This script does NOT reimplement any physics. It only imports and calls
existing functions:

    log_likelihood      (code/monte_carlo_sampling.py)
        -- which itself calls eft_dimension6_xsec (code/eft_cross_sections.py)
           and hadronic_cross_section (code/pdf_integration.py) and
           calculate_uhdm_flux (code/observational_constraints.py)

It does NOT run the legacy MCMC sampler (`run_mcmc` / `emcee.EnsembleSampler`
are never called here) and does NOT require any package beyond what the
repository's existing `code/requirements.txt` already declares.

sigma_jet and B_jet are not active inference parameters (see
docs/METHODS.md, Section 4) -- any fixed values within their existing
prior range are used here only because log_likelihood's signature
requires a 4-tuple; they do not affect the result (this script's own
output demonstrates that explicitly, see the "sigma_jet/B_jet
invariance check" below).

Usage:
    python code/deterministic_grid_reproduction.py
"""

import os
import sys
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from monte_carlo_sampling import log_likelihood, PRIOR_LO, PRIOR_HI

# Fixed benchmark, EFT-validity, and threshold boundaries, reproduced here
# only as closed-form arithmetic on the same constant already defined
# inside monte_carlo_sampling.log_likelihood (sqrt_s_hat = 3.16e8 GeV) --
# not a re-derivation, a direct restatement for reporting purposes.
SQRT_S_HAT = 3.16e8  # GeV, must match the constant used inside log_likelihood


def accepted_region_boundaries():
    """Closed-form accepted-region boundaries (not fit/estimated)."""
    log_lambda_boundary = np.log10(SQRT_S_HAT)
    log_mx_boundary = np.log10(SQRT_S_HAT / 2.0)
    return log_lambda_boundary, log_mx_boundary


def run_grid(step=0.25, sigma_jet_fixed=1.5, b_jet_fixed=1.0):
    """
    Evaluate log_likelihood on a grid spanning the established prior box
    (PRIOR_LO/PRIOR_HI, imported directly from monte_carlo_sampling.py --
    not redefined here), at 0.25-dex spacing in each of log10(Lambda) and
    log10(m_X), matching the grid resolution used in the reported result.
    """
    log_lambda_vals = np.arange(PRIOR_LO[0], PRIOR_HI[0] + 1e-9, step)
    log_mx_vals = np.arange(PRIOR_LO[1], PRIOR_HI[1] + 1e-9, step)

    results = []
    for log_lambda in log_lambda_vals:
        for log_mx in log_mx_vals:
            ll = log_likelihood((log_lambda, log_mx, sigma_jet_fixed, b_jet_fixed))
            results.append((log_lambda, log_mx, ll))
    return results, log_lambda_vals, log_mx_vals


def sigma_jet_b_jet_invariance_check(n_points=5, seed=0):
    """
    Explicit, executed check that sigma_jet/B_jet do not affect
    log_likelihood: evaluate the same (Lambda, m_X) point at several
    different (sigma_jet, B_jet) values and confirm the result is
    identical every time.
    """
    rng = np.random.default_rng(seed)
    log_lambda, log_mx = 9.5, 8.0  # an arbitrary point inside the accepted region
    values = []
    for _ in range(n_points):
        log_sigma = rng.uniform(PRIOR_LO[2], PRIOR_HI[2])
        log_b = rng.uniform(PRIOR_LO[3], PRIOR_HI[3])
        ll = log_likelihood((log_lambda, log_mx, log_sigma, log_b))
        values.append(ll)
    return values


def main():
    print("=" * 70)
    print("Deterministic Grid Reproduction")
    print("=" * 70)

    log_lambda_boundary, log_mx_boundary = accepted_region_boundaries()
    print(f"\nFixed benchmark sqrt(s_hat) = {SQRT_S_HAT:.3e} GeV")
    print(f"EFT-validity boundary:  log10(Lambda) > {log_lambda_boundary:.4f}")
    print(f"Threshold boundary:     log10(m_X)    < {log_mx_boundary:.4f}")

    print("\n--- sigma_jet/B_jet invariance check ---")
    inv_values = sigma_jet_b_jet_invariance_check()
    print(f"log_likelihood at fixed (Lambda,m_X), varying (sigma_jet,B_jet):")
    print(f"  values = {inv_values}")
    print(f"  all identical: {len(set(inv_values)) == 1}")

    print("\n--- 221-point grid characterization ---")
    results, log_lambda_vals, log_mx_vals = run_grid(step=0.25)
    n_points = len(results)
    print(f"Grid points evaluated: {n_points} "
          f"({len(log_lambda_vals)} x {len(log_mx_vals)})")

    zero_count = sum(1 for (_, _, ll) in results if ll == 0.0)
    invalid_count = sum(1 for (_, _, ll) in results if ll == -1e10)
    other_count = n_points - zero_count - invalid_count

    print(f"  log_L == 0.0  (accepted, no observational penalty): {zero_count}")
    print(f"  log_L == -1e10 (rejected -- EFT-invalid or below threshold): {invalid_count}")
    print(f"  log_L at any other value (a genuine intermediate/soft penalty): {other_count}")

    print("\n--- Summary ---")
    if other_count == 0:
        print("CONFIRMED: no grid point produces an intermediate (soft, data-driven)")
        print("likelihood value. Every point is either fully accepted (log_L=0.0,")
        print("no observational channel exceeds its hardcoded limit) or fully")
        print("rejected by a kinematic/EFT-validity cutoff (log_L=-1e10). This")
        print("matches the reported deterministic characterization exactly.")
    else:
        print("UNEXPECTED: intermediate likelihood values were found. This result")
        print("does NOT match the previously reported deterministic characterization")
        print("and should be investigated before being treated as confirmatory.")

    return {
        "log_lambda_boundary": log_lambda_boundary,
        "log_mx_boundary": log_mx_boundary,
        "n_points": n_points,
        "zero_count": zero_count,
        "invalid_count": invalid_count,
        "other_count": other_count,
        "invariance_values": inv_values,
    }


if __name__ == "__main__":
    main()
