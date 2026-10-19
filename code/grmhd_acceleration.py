"""
GRMHD Acceleration Module
=========================

Calculates maximum Lorentz factors for electrons and protons accelerated
in magnetically-dominated AGN jets.

Electrons: magnetic reconnection model for sigma >> 1 regime
  gamma_e_max = eta * sigma   (Sironi & Spitkovsky 2014; Werner et al. 2016)

Protons: Hillas criterion
  gamma_p_max = (e * B * R_jet) / (m_p * c^2)

Parameter ranges from GRMHD simulations of MAD jets:
  B     ~ 0.1 -- 100 G       (Zamaninasab et al. 2014)
  sigma ~ 10  -- 1000         (Tchekhovskoy et al. 2011)
  R_jet ~ 10^14 -- 10^16 cm  (sub-pc to pc scale)
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Physical constants (CGS)
# ---------------------------------------------------------------------------
c_cgs = 2.998e10       # speed of light  [cm/s]
m_e   = 9.109e-28      # electron mass   [g]
m_p   = 1.673e-24      # proton mass     [g]
e_cgs = 4.803e-10      # elementary charge [esu]

# ---------------------------------------------------------------------------
# Derived scales
# ---------------------------------------------------------------------------
m_e_c2 = m_e * c_cgs**2   # electron rest energy [erg]
m_p_c2 = m_p * c_cgs**2   # proton rest energy   [erg]

# Project paths
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, '..', 'figures')


def magnetization_parameter(B, n, gamma_bulk):
    """
    Magnetization parameter  sigma = B^2 / (4 pi rho c^2 gamma_bulk^2).

    Parameters
    ----------
    B : float or array
        Magnetic field strength [G].
    n : float or array
        Proton number density [cm^-3].
    gamma_bulk : float
        Bulk Lorentz factor of the jet.

    Returns
    -------
    sigma : float or array
    """
    rho = n * m_p
    return B**2 / (4.0 * np.pi * rho * c_cgs**2 * gamma_bulk**2)


def gamma_electron_reconnection(sigma, eta=0.1):
    """
    Maximum electron Lorentz factor from magnetic reconnection in the
    sigma >> 1 regime.

    gamma_e_max = eta * sigma

    In the high-magnetisation limit the reconnection electric field
    accelerates electrons up to a fraction eta of the available magnetic
    energy per particle.  PIC simulations give eta ~ 0.1--0.3
    (Sironi & Spitkovsky 2014, ApJ 783, L21).

    Parameters
    ----------
    sigma : float or array
        Magnetization parameter (dimensionless).
    eta : float
        Reconnection efficiency, default 0.1.

    Returns
    -------
    gamma_e_max : float or array
    """
    return eta * np.asarray(sigma, dtype=float)


def gamma_proton_hillas(B, R_jet):
    """
    Maximum proton Lorentz factor from the Hillas criterion.

    gamma_p_max = (e * B * R_jet) / (m_p * c^2)

    A charged particle can be confined and accelerated only while its
    Larmor radius fits inside the acceleration region R_jet.

    Parameters
    ----------
    B : float or array
        Magnetic field strength [G].
    R_jet : float or array
        Jet radius (acceleration region size) [cm].

    Returns
    -------
    gamma_p_max : float or array
    """
    return (e_cgs * np.asarray(B, dtype=float) *
            np.asarray(R_jet, dtype=float)) / m_p_c2


def plot_gamma_max_vs_sigma(output_file=None):
    """
    Generate Figure 2: maximum Lorentz factors vs magnetization / B-field.

    Left panel  -- gamma_e_max vs sigma for several eta values.
    Right panel -- gamma_p_max vs B for several R_jet values.

    Saved to figures/fig2_gamma_max.pdf.
    """
    if output_file is None:
        output_file = os.path.join(FIGURES_DIR, 'fig2_gamma_max.pdf')

    sigma_arr = np.logspace(1, 3, 200)         # 10 -- 1000
    B_arr     = np.logspace(-1, 2, 200)         # 0.1 -- 100 G

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # -- Left panel: electron reconnection ---------------------------------
    for eta, ls, col in [(0.05, '--', 'C0'),
                         (0.1,  '-',  'C1'),
                         (0.3,  '-.',  'C2')]:
        gamma_e = gamma_electron_reconnection(sigma_arr, eta=eta)
        ax1.loglog(sigma_arr, gamma_e, ls, color=col, lw=2,
                   label=rf'$\eta = {eta}$')

    ax1.set_xlabel(r'Magnetization $\sigma$', fontsize=13)
    ax1.set_ylabel(r'$\gamma_{e,\mathrm{max}}$', fontsize=13)
    ax1.set_title('Electron reconnection', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.grid(True, which='both', alpha=0.25)

    # -- Right panel: proton Hillas criterion --------------------------------
    for R_jet, ls, col, lbl in [
            (1e14, '--', 'C3', r'$R_\mathrm{jet}=10^{14}\;\mathrm{cm}$'),
            (1e15, '-',  'C4', r'$R_\mathrm{jet}=10^{15}\;\mathrm{cm}$'),
            (1e16, '-.', 'C5', r'$R_\mathrm{jet}=10^{16}\;\mathrm{cm}$')]:
        gamma_p = gamma_proton_hillas(B_arr, R_jet)
        ax2.loglog(B_arr, gamma_p, ls, color=col, lw=2, label=lbl)

    ax2.set_xlabel(r'Magnetic field $B$ [G]', fontsize=13)
    ax2.set_ylabel(r'$\gamma_{p,\mathrm{max}}$', fontsize=13)
    ax2.set_title('Proton Hillas criterion', fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(True, which='both', alpha=0.25)

    fig.suptitle('Maximum Lorentz Factors in AGN Jets (Figure 2)',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Figure saved: {output_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print("GRMHD Acceleration Module")
    print("=" * 55)

    # --- Example calculation ------------------------------------------------
    B_fid     = 10.0    # G
    R_jet_fid = 1e15    # cm
    n_fid     = 1e2     # cm^-3
    gamma_bulk = 10.0

    sigma = magnetization_parameter(B_fid, n_fid, gamma_bulk)
    g_e   = gamma_electron_reconnection(sigma, eta=0.1)
    g_p   = gamma_proton_hillas(B_fid, R_jet_fid)

    print(f"\nFiducial parameters:")
    print(f"  B        = {B_fid} G")
    print(f"  R_jet    = {R_jet_fid:.1e} cm")
    print(f"  n        = {n_fid:.1e} cm^-3")
    print(f"  gamma_bk = {gamma_bulk}")
    print(f"\nDerived:")
    print(f"  sigma          = {sigma:.1f}")
    print(f"  gamma_e_max    = {g_e:.1e}  (reconnection, eta=0.1)")
    print(f"  gamma_p_max    = {g_p:.1e}  (Hillas criterion)")

    # --- Scan over parameter space ------------------------------------------
    print("\nParameter survey:")
    for B in [0.1, 1.0, 10.0, 100.0]:
        for R in [1e14, 1e15, 1e16]:
            gp = gamma_proton_hillas(B, R)
            print(f"  B={B:6.1f} G, R_jet={R:.0e} cm  =>  gamma_p_max = {gp:.2e}")

    # --- Generate Figure 2 --------------------------------------------------
    print("\nGenerating Figure 2 ...")
    plot_gamma_max_vs_sigma()
    print("Done.")
