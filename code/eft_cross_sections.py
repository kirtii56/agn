"""
EFT Cross Section Module
========================

Computes production cross sections for ultra-heavy dark matter (UHDM)
using an effective-field-theory contact interaction (dimension-6 operator).

    L_eff  =  (1 / Lambda^2)  (q-bar q)(X-bar X)

Partonic cross section:

    sigma_XX(s) = (s / Lambda^4) * (1 / (48 pi)) * beta(s, m_X)**3

where beta = sqrt(1 - 4 m_X^2 / s) is the velocity factor, and the EFT
is valid only when sqrt(s) < Lambda.

Mass range  :  m_X    ~ 10^6  -- 10^10 GeV  (= 10^15 -- 10^19 eV)
Cutoff range:  Lambda ~ 10^8  -- 10^11 GeV  (= 10^17 -- 10^20 eV)
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
GeV2_to_cm2 = 3.894e-28          # 1 GeV^{-2} = 3.894 x 10^{-28} cm^2
M_Planck    = 1.221e19            # Planck mass [GeV]

# Project paths
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, '..', 'figures')


# ---------------------------------------------------------------------------
# Cross-section calculation
# ---------------------------------------------------------------------------
def eft_dimension6_xsec(s, M_X, Lambda):
    """
    EFT contact-interaction cross section for q qbar -> X Xbar.

    sigma = (s / Lambda^4) * (1 / (48 pi)) * beta**3

    Parameters
    ----------
    s : float or array
        Partonic centre-of-mass energy squared [GeV^2].
    M_X : float
        UHDM particle mass [GeV].
    Lambda : float
        EFT cutoff scale [GeV].

    Returns
    -------
    xsec : float or array
        Cross section [cm^2].  Returns 0 below threshold, NaN when EFT
        is invalid (sqrt(s) > Lambda).
    """
    s = np.asarray(s, dtype=float)
    scalar = (s.ndim == 0)
    s = np.atleast_1d(s)

    xsec = np.zeros_like(s)

    # Threshold check: s >= (2 m_X)^2
    above_threshold = s >= 4.0 * M_X**2

    # EFT validity: sqrt(s) < Lambda
    eft_valid = np.sqrt(s) < Lambda

    mask = above_threshold & eft_valid
    if np.any(mask):
        beta = np.sqrt(1.0 - 4.0 * M_X**2 / s[mask])
        xsec[mask] = (s[mask] / Lambda**4) * (1.0 / (48.0 * np.pi)) * beta**3
        xsec[mask] *= GeV2_to_cm2   # convert to cm^2

    # Flag EFT-invalid region
    xsec[above_threshold & ~eft_valid] = np.nan

    if scalar:
        return float(xsec[0])
    return xsec


def eft_validity_check(sqrt_s, Lambda):
    """Return True where the EFT is valid (sqrt(s) < Lambda)."""
    return np.asarray(sqrt_s) < Lambda


# ---------------------------------------------------------------------------
# Figure 1 — EFT validity regions
# ---------------------------------------------------------------------------
def plot_eft_validity(output_file=None):
    """
    Figure 1: EFT validity map in the (m_X, Lambda) plane.

    Current locked criterion: EFT validity requires sqrt(s_hat) < Lambda,
    evaluated at the fixed benchmark sqrt(s_hat) = 3.16e8 GeV used
    throughout the reported deterministic inference. This condition does
    NOT depend on m_X, so the validity boundary is a HORIZONTAL line at
    Lambda = 3.16e8 GeV, not a diagonal Lambda=2*m_X line.

    Also shown, as context only (not additional validity criteria beyond
    what is already used in the reported accepted region): the AGN
    sqrt(s) band from which the fixed benchmark was drawn, the Planck
    mass for scale comparison, reference Lambda values matching the
    scanned prior, and a vertical line at the kinematic (pair-production)
    threshold m_X = sqrt(s_hat)/2 -- the second of the two conditions
    (together with the horizontal EFT-validity line) that define the
    reported accepted region.
    """
    if output_file is None:
        output_file = os.path.join(FIGURES_DIR, 'fig1_eft_validity_corrected.pdf')

    SQRT_S_HAT = 3.16e8  # GeV -- fixed benchmark used throughout the reported inference

    m_X_arr = np.logspace(6, 10, 300)    # GeV

    fig, ax = plt.subplots(figsize=(8, 6))

    # EFT-validity boundary (current, locked criterion): sqrt(s_hat) < Lambda,
    # independent of m_X -- a horizontal line.
    ax.axhline(SQRT_S_HAT, color='k', lw=2,
               label=r'EFT validity boundary: $\Lambda = \sqrt{\hat{s}} = 3.16\times10^8$ GeV')

    # Fill valid region (Lambda above the horizontal boundary)
    ax.axhspan(SQRT_S_HAT, 1e13, alpha=0.12, color='green', label='EFT valid')

    # Kinematic (pair-production) threshold: m_X < sqrt(s_hat)/2 -- the second
    # condition (together with EFT validity above) defining the reported
    # accepted region; shown as a vertical line for direct comparison.
    ax.axvline(SQRT_S_HAT / 2.0, color='k', lw=1.5, ls='-.',
               label=r'Threshold: $m_X = \sqrt{\hat{s}}/2 = 1.58\times10^8$ GeV')

    # AGN jet centre-of-mass energy band (context: sqrt(s_hat) was drawn
    # from this band)
    ax.axhspan(1e8, 1e9, color='orange', alpha=0.15,
               label=r'AGN $\sqrt{s}$ band ($10^{17}$--$10^{18}$ eV)')

    # Planck scale
    ax.axhline(M_Planck, color='grey', ls=':', lw=1.5,
               label=r'$M_\mathrm{Pl} = 1.22\times10^{19}$ GeV')

    # Reference Lambda values matching the scanned prior
    for L, col in [(1e8, 'C0'), (1e9, 'C1'), (1e10, 'C2'), (1e11, 'C3')]:
        ax.axhline(L, color=col, ls='--', lw=1, alpha=0.6,
                   label=rf'$\Lambda = 10^{{{int(np.log10(L))}}}$ GeV')

    ax.set_xlabel(r'UHDM mass $m_X$ [GeV]', fontsize=13)
    ax.set_ylabel(r'EFT cutoff $\Lambda$ [GeV]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('EFT Validity Region at the Fixed Benchmark (Figure 1)', fontsize=13)
    ax.set_xlim(1e6, 1e10)
    ax.set_ylim(1e6, 1e13)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, which='both', alpha=0.2)

    fig.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Figure saved: {output_file}")


# ---------------------------------------------------------------------------
# Figure 3 — Production cross section vs mass
# ---------------------------------------------------------------------------
def plot_production_rate(output_file=None):
    """
    Figure 3: sigma_XX vs m_X at several Lambda values.

    Uses a fixed sqrt(s) representative of AGN jet collisions.
    """
    if output_file is None:
        output_file = os.path.join(FIGURES_DIR, 'fig3_production_rate.pdf')

    m_X_arr = np.logspace(6, 10, 400)  # GeV

    fig, ax = plt.subplots(figsize=(8, 6))

    for Lambda, col, ls in [
            (1e8,  'C0', '-'),
            (1e9,  'C1', '--'),
            (1e10, 'C2', '-.'),
            (1e11, 'C3', ':')]:
        # Fix sqrt(s) at the geometric mean of the AGN band
        sqrt_s = 3.16e8  # GeV  (~10^{17.5} eV)
        s_val  = sqrt_s**2
        xsec = np.array([eft_dimension6_xsec(s_val, mx, Lambda)
                          for mx in m_X_arr])

        # Mask out zero / NaN for clean log-log plot
        valid = np.isfinite(xsec) & (xsec > 0)
        if np.any(valid):
            ax.loglog(m_X_arr[valid], xsec[valid], ls, color=col, lw=2,
                      label=rf'$\Lambda = 10^{{{int(np.log10(Lambda))}}}$ GeV')

    ax.set_xlabel(r'UHDM mass $m_X$ [GeV]', fontsize=13)
    ax.set_ylabel(r'$\sigma_{XX}$ [cm$^{2}$]', fontsize=13)
    ax.set_title(
        r'Production Cross Section at $\sqrt{s}=3.16\times10^{8}$ GeV '
        '(Figure 3)', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, which='both', alpha=0.2)

    fig.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Figure saved: {output_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print("EFT Cross Section Module")
    print("=" * 55)

    # --- Example calculation ------------------------------------------------
    m_X    = 1e8     # GeV  (=10^{17} eV)
    Lambda = 1e10    # GeV  (=10^{19} eV)
    sqrt_s = 3e8     # GeV
    s_val  = sqrt_s**2

    xsec = eft_dimension6_xsec(s_val, m_X, Lambda)
    valid = eft_validity_check(sqrt_s, Lambda)

    print(f"\nExample:")
    print(f"  m_X    = {m_X:.1e} GeV")
    print(f"  Lambda = {Lambda:.1e} GeV")
    print(f"  sqrt(s)= {sqrt_s:.1e} GeV")
    print(f"  EFT valid: {valid}")
    print(f"  sigma_XX  = {xsec:.3e} cm^2")

    # --- Survey over mass range ---------------------------------------------
    print("\nCross-section survey (Lambda = 1e10 GeV, sqrt_s = 3e8 GeV):")
    for mx in [1e6, 1e7, 1e8, 1e9]:
        xs = eft_dimension6_xsec(s_val, mx, Lambda)
        tag = '' if np.isfinite(xs) else '  [EFT invalid or below threshold]'
        print(f"  m_X = {mx:.0e} GeV  =>  sigma = {xs:.3e} cm^2{tag}")

    # --- Generate figures ---------------------------------------------------
    print("\nGenerating Figure 1 ...")
    plot_eft_validity()
    print("Generating Figure 3 ...")
    plot_production_rate()
    print("Done.")
