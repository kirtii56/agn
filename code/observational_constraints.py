"""
Observational Constraints Module
================================

Implements multi-messenger upper limits on signals from UHDM production
in AGN jets.

Experiments and parameterised limits used in this analysis
(see docs/METHODS.md for detailed notes on each):

1. Fermi-LAT (3FHL, 7-year Pass 8)
      E^2 dN/dE  <  10^{-12}  erg cm^{-2} s^{-1}
      Energy range: 10 GeV -- 2 TeV
      Ref: Ajello et al. 2017, ApJS 232, 18

2. IceCube (10-year point-source stacking)
      phi_nu  <  10^{-18}  GeV cm^{-2} s^{-1} sr^{-1}  per AGN
      Energy range: 1 TeV -- 10 PeV
      Ref: IceCube Collaboration 2023, ApJ 954, 75

3. Pierre Auger (UHECR spectrum)
      Upper limits on flux above E > 10^{19.5} eV
      Ref: Pierre Auger Collaboration 2020, PRL 125, 121106
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Constants and paths
# ---------------------------------------------------------------------------
eV_to_GeV = 1.0e-9
GeV_to_erg = 1.602e-3       # 1 GeV = 1.602e-3 erg
Mpc_to_cm  = 3.086e24       # 1 Mpc = 3.086e24 cm

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, '..', 'figures')


class MultiMessengerConstraints:
    """
    Container for multi-messenger upper limits on UHDM-related signals.

    All limits are stored in *instrument-native* units; helper methods
    convert a model prediction to the same units for comparison.
    """

    def __init__(self):
        # Fermi-LAT  --  E^2 dN/dE limit
        self.fermi_E2dNdE_limit = 1e-12          # erg cm^{-2} s^{-1}
        self.fermi_energy_GeV   = (10.0, 2e3)    # 10 GeV -- 2 TeV

        # IceCube  --  per-AGN neutrino flux
        self.icecube_phi_limit  = 1e-18           # GeV cm^{-2} s^{-1} sr^{-1}
        self.icecube_energy_GeV = (1e3, 1e7)      # 1 TeV -- 10 PeV

        # Pierre Auger  --  integral flux above E_threshold
        self.auger_E_threshold_eV = 10**19.5      # 31.6 EeV
        self.auger_E_threshold_GeV = self.auger_E_threshold_eV * eV_to_GeV
        self.auger_energy_GeV = (self.auger_E_threshold_GeV, 1e12)  # up to ~10^{21} eV

        # Approximate Auger integral flux limit above threshold  [cm^{-2} s^{-1} sr^{-1}]
        self.auger_integral_flux_limit = 5e-40    # very rough; see METHODS.md

    # -----------------------------------------------------------------------
    def fermi_excluded(self, E2dNdE_predicted):
        """True if predicted E^2 dN/dE exceeds Fermi-LAT limit."""
        return E2dNdE_predicted > self.fermi_E2dNdE_limit

    def icecube_excluded(self, phi_nu_predicted):
        """True if predicted per-AGN neutrino flux exceeds IceCube limit."""
        return phi_nu_predicted > self.icecube_phi_limit

    def auger_excluded(self, integral_flux_predicted):
        """True if predicted integral UHECR flux exceeds Auger limit."""
        return integral_flux_predicted > self.auger_integral_flux_limit


def calculate_uhdm_flux(M_X_GeV, production_rate, distance_Mpc=100.0):
    """
    Estimate observable fluxes from UHDM decay products at Earth.

    Assumes prompt decay  X -> gamma gamma  (or hadronic cascade) with
    branching fractions  BR_gamma = 0.6, BR_nu = 0.3, BR_CR = 0.1.

    Parameters
    ----------
    M_X_GeV : float
        UHDM mass [GeV].
    production_rate : float
        UHDM production rate [particles s^{-1}].
    distance_Mpc : float
        Luminosity distance to the AGN [Mpc].

    Returns
    -------
    dict with keys  'E2dNdE_gamma' [erg/cm^2/s],
                    'phi_nu'       [GeV/cm^2/s/sr],
                    'integral_CR'  [cm^{-2} s^{-1} sr^{-1}]
    """
    d_cm = distance_Mpc * Mpc_to_cm

    # Total energy release per unit time [erg/s]
    L_total = production_rate * M_X_GeV * GeV_to_erg

    # Isotropic-equivalent flux at Earth [erg/cm^2/s]
    F_total = L_total / (4.0 * np.pi * d_cm**2)

    # Branching fractions (illustrative; model-dependent)
    BR_gamma = 0.6
    BR_nu    = 0.3
    BR_CR    = 0.1

    # E^2 dN/dE for gamma-rays (approximate: most energy near M_X/2)
    E2dNdE_gamma = F_total * BR_gamma

    # Per-AGN neutrino flux  [GeV / cm^2 / s / sr]
    # Convert erg -> GeV, then divide by 4 pi sr
    phi_nu = (F_total * BR_nu / GeV_to_erg) / (4.0 * np.pi)

    # Integral cosmic-ray flux  [cm^{-2} s^{-1} sr^{-1}]
    E_CR_erg = M_X_GeV * GeV_to_erg * BR_CR
    n_CR     = production_rate * BR_CR
    integral_CR = n_CR / (4.0 * np.pi * d_cm**2) / (4.0 * np.pi)

    return {
        'E2dNdE_gamma': E2dNdE_gamma,
        'phi_nu':        phi_nu,
        'integral_CR':   integral_CR,
    }


# ---------------------------------------------------------------------------
# Figure 4 — multi-messenger constraints
# ---------------------------------------------------------------------------
def plot_constraints_summary(output_file=None):
    """
    Figure 4: multi-messenger constraint summary.

    Shows experimental upper limits as horizontal bands and an example
    model prediction curve as a function of UHDM mass.
    """
    if output_file is None:
        output_file = os.path.join(FIGURES_DIR, 'fig4_constraints.pdf')

    constraints = MultiMessengerConstraints()

    # Energy axis in GeV (spans all three experiments)
    E_GeV = np.logspace(1, 12, 500)
    E_eV  = E_GeV / eV_to_GeV   # in eV

    fig, ax = plt.subplots(figsize=(10, 6.5))

    # --- Fermi-LAT excluded region -----------------------------------------
    f_lo, f_hi = constraints.fermi_energy_GeV
    f_mask = (E_GeV >= f_lo) & (E_GeV <= f_hi)
    ax.fill_between(E_eV[f_mask],
                    constraints.fermi_E2dNdE_limit, 1e-8,
                    color='C0', alpha=0.25,
                    label=rf'Fermi-LAT: $E^2 dN/dE < 10^{{-12}}$ erg/cm$^2$/s')
    ax.hlines(constraints.fermi_E2dNdE_limit,
              f_lo / eV_to_GeV, f_hi / eV_to_GeV,
              colors='C0', lw=2)

    # --- IceCube excluded region -------------------------------------------
    ic_lo, ic_hi = constraints.icecube_energy_GeV
    ic_mask = (E_GeV >= ic_lo) & (E_GeV <= ic_hi)
    # Convert IceCube limit to erg/cm^2/s for same y-axis
    ic_erg = constraints.icecube_phi_limit * GeV_to_erg * 4.0 * np.pi
    ax.fill_between(E_eV[ic_mask],
                    ic_erg, 1e-8,
                    color='C2', alpha=0.20,
                    label=rf'IceCube: $\phi_\nu < 10^{{-18}}$ GeV/cm$^2$/s/sr')
    ax.hlines(ic_erg,
              ic_lo / eV_to_GeV, ic_hi / eV_to_GeV,
              colors='C2', lw=2)

    # --- Auger excluded region ---------------------------------------------
    au_lo, au_hi = constraints.auger_energy_GeV
    au_mask = (E_GeV >= au_lo) & (E_GeV <= au_hi)
    # Convert Auger integral flux to erg-equivalent for y-axis
    au_erg = constraints.auger_integral_flux_limit * (constraints.auger_E_threshold_GeV
                                                       * GeV_to_erg) * 4 * np.pi
    if au_erg > 0:
        ax.fill_between(E_eV[au_mask],
                        au_erg, 1e-8,
                        color='C3', alpha=0.18,
                        label=r'Auger: UHECR above $10^{19.5}$ eV')
        ax.hlines(au_erg,
                  au_lo / eV_to_GeV, au_hi / eV_to_GeV,
                  colors='C3', lw=2)

    # --- Example model prediction ------------------------------------------
    # Power-law decline in E^2 dN/dE with mass (illustrative)
    E_ref = 1e6 / eV_to_GeV   # reference energy 10^{15} eV
    model_flux = 1e-13 * (E_eV / E_ref)**(-0.5)
    ax.loglog(E_eV, model_flux, 'k-', lw=2.0, alpha=0.7,
              label=r'Model prediction ($\dot{N}_X = 10^{12}$ s$^{-1}$, $d=100$ Mpc)')

    ax.set_xlabel('Energy [eV]', fontsize=13)
    ax.set_ylabel(r'$E^2\,dN/dE$ [erg cm$^{-2}$ s$^{-1}$]', fontsize=13)
    ax.set_title('Multi-Messenger Constraints on UHDM Production (Figure 4)',
                 fontsize=13)
    ax.set_xlim(1e10, 1e21)
    ax.set_ylim(1e-16, 1e-8)
    ax.legend(fontsize=9, loc='upper right')
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
    print("Observational Constraints Module")
    print("=" * 55)

    constraints = MultiMessengerConstraints()

    print("\nExperimental upper limits:")
    print(f"  Fermi-LAT : E^2 dN/dE < {constraints.fermi_E2dNdE_limit:.0e} "
          f"erg/cm^2/s  ({constraints.fermi_energy_GeV[0]:.0f} -- "
          f"{constraints.fermi_energy_GeV[1]:.0e} GeV)")
    print(f"  IceCube   : phi_nu < {constraints.icecube_phi_limit:.0e} "
          f"GeV/cm^2/s/sr  ({constraints.icecube_energy_GeV[0]:.0e} -- "
          f"{constraints.icecube_energy_GeV[1]:.0e} GeV)")
    print(f"  Auger     : threshold {constraints.auger_E_threshold_eV:.2e} eV "
          f"= {constraints.auger_E_threshold_GeV:.2e} GeV")

    # --- Example model point ------------------------------------------------
    M_X  = 1e8     # GeV
    Ndot = 1e12    # particles/s
    dist = 100.0   # Mpc

    fluxes = calculate_uhdm_flux(M_X, Ndot, dist)
    print(f"\nModel point: M_X = {M_X:.0e} GeV, Ndot = {Ndot:.0e} s^-1, "
          f"d = {dist} Mpc")
    print(f"  E^2 dN/dE (gamma) = {fluxes['E2dNdE_gamma']:.2e} erg/cm^2/s"
          f"  {'EXCLUDED' if constraints.fermi_excluded(fluxes['E2dNdE_gamma']) else 'allowed'}")
    print(f"  phi_nu            = {fluxes['phi_nu']:.2e} GeV/cm^2/s/sr"
          f"  {'EXCLUDED' if constraints.icecube_excluded(fluxes['phi_nu']) else 'allowed'}")
    print(f"  integral CR flux  = {fluxes['integral_CR']:.2e} cm^-2/s/sr"
          f"  {'EXCLUDED' if constraints.auger_excluded(fluxes['integral_CR']) else 'allowed'}")

    # --- Generate Figure 4 --------------------------------------------------
    print("\nGenerating Figure 4 ...")
    plot_constraints_summary()
    print("Done.")
