"""
PDF Integration Module
======================

Implements a simplified parton luminosity function for pp collisions at
extreme centre-of-mass energies relevant to AGN jet particle collisions.

Uses a parametric quark-antiquark luminosity model. This model is NOT
CT18 NNLO and is NOT a validated physical PDF set; it is used here only
for exploratory/diagnostic purposes, and no claim of CT18 NNLO fidelity
is made:

    dL/dtau  ~  A * tau^{-1} * (1 - tau)^n

where tau = s_hat / s, and A, n are fitted parameters.  This functional
form captures the dominant behaviour of the quark luminosity at large
momentum fractions without requiring the full LHAPDF grid.

Covered energies:  sqrt(s) ~ 10^{17} -- 10^{18} eV  (= 10^8 -- 10^9 GeV)

Reference:  Hou et al. 2021, Phys. Rev. D 103, 014013 (CT18 NNLO)
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Project paths
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, '..', 'figures')

# Conversion
GeV2_to_cm2 = 3.894e-28   # 1 GeV^{-2} = 3.894e-28 cm^2


class CT18Luminosity:
    """
    Parametric quark-antiquark luminosity model, for exploratory/
    diagnostic use only. This is NOT CT18 NNLO and NOT a validated
    physical PDF set -- no claim of fidelity to any specific PDF set
    is made.

    The luminosity is modelled as

        dL/dtau = A * tau^{a} * (1 - tau)^n

    where tau = M^2 / s.  Default parameters (A=0.5, a=-1.0, n=5.0) give
    a reasonable approximation to the qq-bar channel at high Q^2.

    Parameters
    ----------
    A : float
        Overall normalisation.
    a : float
        Small-tau exponent (controls the rise at low x), default -1.
    n : float
        Large-tau exponent (controls the fall at high x), default 5.
    """

    def __init__(self, A=0.5, a=-1.0, n=5.0):
        self.A = A
        self.a = a
        self.n = n

    def dLdtau(self, tau):
        """
        Parton luminosity dL/dtau at given tau = s_hat / s.

        Returns 0 outside [0, 1].
        """
        tau = np.asarray(tau, dtype=float)
        out = np.zeros_like(tau)
        mask = (tau > 0) & (tau < 1)
        out[mask] = self.A * tau[mask]**self.a * (1.0 - tau[mask])**self.n
        return out

    def integrated_luminosity(self, tau_min, tau_max=1.0):
        """
        Integrate dL/dtau from tau_min to tau_max.

        Parameters
        ----------
        tau_min : float
            Lower integration limit (= M_X^2 / s for threshold production).
        tau_max : float
            Upper limit (default 1).

        Returns
        -------
        L : float
            Integrated luminosity (dimensionless).
        """
        def integrand(tau):
            return float(self.dLdtau(np.array([tau]))[0])

        result, _ = quad(integrand, tau_min, min(tau_max, 1.0 - 1e-15),
                         limit=200)
        return result


def hadronic_cross_section(sigma_hat_func, s, M_X, lumi=None):
    """
    Convolve a partonic cross section with the parton luminosity to obtain
    the hadronic cross section.

    sigma(pp -> XX) = integral_{tau_min}^{1} dtau  dL/dtau  *  sigma_hat(tau*s)

    where tau_min = (2 M_X)^2 / s.

    Parameters
    ----------
    sigma_hat_func : callable
        sigma_hat(s_hat, M_X) returning the partonic cross section [cm^2].
    s : float
        pp centre-of-mass energy squared [GeV^2].
    M_X : float
        UHDM mass [GeV].
    lumi : CT18Luminosity or None
        Parton luminosity object.  If None, uses default CT18 parameterisation.

    Returns
    -------
    sigma_pp : float
        Hadronic cross section [cm^2].
    """
    if lumi is None:
        lumi = CT18Luminosity()

    tau_min = (2.0 * M_X)**2 / s
    if tau_min >= 1.0:
        return 0.0

    def integrand(tau):
        s_hat = tau * s
        sig_hat = sigma_hat_func(s_hat, M_X)
        if not np.isfinite(sig_hat) or sig_hat <= 0:
            return 0.0
        dl = float(lumi.dLdtau(np.array([tau]))[0])
        return dl * sig_hat

    result, _ = quad(integrand, tau_min, 1.0 - 1e-15, limit=200)
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print("PDF Integration Module")
    print("=" * 55)

    lumi = CT18Luminosity()

    # --- Show luminosity at representative tau values -----------------------
    print("\nParametric CT18 NNLO luminosity  dL/dtau = A * tau^a * (1-tau)^n")
    print(f"  A = {lumi.A},  a = {lumi.a},  n = {lumi.n}\n")

    tau_vals = np.logspace(-8, -1, 8)
    print(f"  {'tau':>12s}   {'dL/dtau':>12s}")
    print(f"  {'-'*12}   {'-'*12}")
    for t in tau_vals:
        dl = float(lumi.dLdtau(np.array([t]))[0])
        print(f"  {t:12.3e}   {dl:12.3e}")

    # --- Integrated luminosity for AGN-jet energies -------------------------
    print("\nIntegrated luminosity above tau_min for sqrt(s) = 3.16e8 GeV:")
    sqrt_s = 3.16e8   # GeV
    s_val  = sqrt_s**2
    for m_X in [1e6, 1e7, 1e8]:
        tau_min = (2.0 * m_X)**2 / s_val
        L = lumi.integrated_luminosity(tau_min)
        print(f"  m_X = {m_X:.0e} GeV  (tau_min = {tau_min:.3e})  =>  L = {L:.4e}")

    # --- Hadronic cross section example -------------------------------------
    # Use the EFT dimension-6 cross section from eft_cross_sections module
    import sys
    sys.path.insert(0, SCRIPT_DIR)
    try:
        from eft_cross_sections import eft_dimension6_xsec

        Lambda = 1e10  # GeV
        print(f"\nHadronic cross section (Lambda = {Lambda:.0e} GeV):")
        for m_X in [1e6, 1e7, 5e7]:
            def sigma_hat(s_hat, mx, L=Lambda):
                return eft_dimension6_xsec(s_hat, mx, L)

            sigma_pp = hadronic_cross_section(sigma_hat, s_val, m_X, lumi)
            print(f"  m_X = {m_X:.0e} GeV  =>  sigma_pp = {sigma_pp:.3e} cm^2")
    except ImportError:
        print("\n  [eft_cross_sections not available; skipping convolution test]")

    print("\nDone.")
