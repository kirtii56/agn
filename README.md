# Ultra-Heavy Dark Matter Production in AGN Jets

[![arXiv](https://img.shields.io/badge/arXiv-preprint-b31b1b.svg)](https://github.com/kirtii56/agn)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Authors:** Kirti Patidar
**Affiliation:** Independent Researcher
**Contact:** kirtiptdr174@gmail.com
**Version:** 2.0-FINAL
**Last Updated:** November 3, 2025

---

## Abstract

This repository contains the complete research package for **"Ultra-Heavy Dark Matter Production in Active Galactic Nucleus Jets: A Falsifiable Multi-Messenger Framework"**.

We investigate whether AGN jets can produce ultra-heavy dark matter (UHDM) with masses m_X = 10^15–10^19 eV. AGN jets accelerate particles to Lorentz factors γ ~ 10^9, enabling center-of-mass collision energies √s ~ 10^17–10^18 eV. We model UHDM production via the dimension-6 contact interaction L_int = (1/Λ²)(q̄q)(X̄X) for q q̄ → X X̄, with Λ treated as a scanned rather than fixed parameter. Energy-loss-corrected calculations give maximum particle Lorentz factors γ_e,max ~ 10^6–10^7 (electrons, energy-loss-limited) and, from the Hillas confinement limit alone (proton curvature-radiation losses are not modeled), γ_p,max ~ 3.2×10^10–3.2×10^14 (protons).

**Key Finding:** A deterministic, grid-verified analysis of the (Λ, m_X) parameter space shows the region simultaneously above the pair-production threshold and within the EFT-validity condition is exactly 3.16×10^8 < Λ ≤ 10^11 GeV, 10^6 ≤ m_X < 1.58×10^8 GeV. Within this entire region, predicted Fermi-LAT, IceCube, and Pierre Auger fluxes remain 2.6 to 24.8 orders of magnitude below the observational limits implemented in this pipeline — current multi-messenger data, as implemented here, provide no data-driven constraint on this model anywhere in the explored parameter space; the region's boundaries are set by kinematics and EFT self-consistency, not by data. A quantitative cosmic-abundance bound (Ω_X^AGN) is not derived in this analysis — see the manuscript for the reasoning.

---

## Repository Structure

```
UHDM_AGN_Paper_Final/
├── manuscript/
│   ├── manuscript_final.tex          # LaTeX source (RevTeX4-2 format)
│   ├── manuscript_final.pdf          # Compiled PDF
│   ├── references.bib                # BibTeX bibliography
│   └── supplementary_material.pdf    # MCMC diagnostics (Figure S1)
├── figures/
│   ├── fig1_eft_validity_corrected.pdf    # EFT cutoff plot (regenerated)
│   ├── fig2_gamma_max.pdf                 # Max Lorentz factors
│   ├── fig3_production_rate.pdf           # Production rates vs mass
│   ├── fig4_constraints.pdf               # Multi-messenger constraints
│   └── figS1_mcmc_diagnostics.pdf         # Convergence diagnostics
├── code/
│   ├── grmhd_acceleration.py         # Lorentz factor calculations
│   ├── eft_cross_sections.py         # EFT production cross sections
│   ├── monte_carlo_sampling.py       # MCMC parameter scan (emcee)
│   ├── pdf_integration.py            # Parametric/exploratory q-qbar luminosity stand-in (not CT18 NNLO)
│   ├── observational_constraints.py  # Fermi/IceCube/Auger limits
│   └── requirements.txt              # Python dependencies
├── data/
│   ├── README_data.txt               # Data documentation
│   └── [Large data files excluded via .gitignore]
│       # Code uses built-in observational limit parameterizations
├── workflow/
│   ├── Dockerfile                    # Reproducibility container
│   ├── run_analysis.sh               # Master execution script
│   └── environment.yml               # Conda environment
├── docs/
│   ├── README.md                     # Detailed documentation
│   ├── metadata.json                 # Project metadata (arXiv)
│   ├── submission_checklist.txt      # Pre-submission validation
│   └── peer_review_response.txt      # Response to reviewer comments
└── LICENSE                           # MIT License
```

---

## Quick Start

### Installation

**Option 1: Conda (Recommended)**
```bash
conda env create -f workflow/environment.yml
conda activate uhdm-agn
```

**Option 2: pip**
```bash
pip install -r code/requirements.txt
```

**Option 3: Docker**
```bash
docker build -t uhdm-agn -f workflow/Dockerfile .
docker run -it -v $(pwd):/workspace uhdm-agn
```

### Running the Analysis

Execute the complete pipeline:
```bash
bash workflow/run_analysis.sh
```

Or run individual modules:
```bash
# GRMHD acceleration
python code/grmhd_acceleration.py

# EFT cross sections
python code/eft_cross_sections.py

# Monte Carlo scan (~1 min for default 10^3 samples)
python code/monte_carlo_sampling.py --seed 42

# Observational constraints
python code/observational_constraints.py
```

### Reproducing Results

For complete step-by-step reproduction of all results:

```bash
# Run Monte Carlo analysis (~1 min default; use --nsteps 10000 for production)
python code/monte_carlo_sampling.py --seed 42

# Generate all figures
python code/generate_figures.py --output figures/

# Compile manuscript
cd manuscript/
latexmk -pdf manuscript_final.tex
```

See [REPRODUCING.md](REPRODUCING.md) for detailed instructions, troubleshooting, and validation steps.

---

## Key Results

- **Maximum Lorentz factors**: γ_e,max ~ 10^6–10^7 (electrons, energy-loss-limited); γ_p,max ~ 3.2×10^10–3.2×10^14 (protons, Hillas confinement limit only — proton curvature-radiation losses not modeled)
- **Center-of-mass energies**: √s ~ 10^17–10^18 eV; fixed benchmark √ŝ ≈ 3.16×10^8 GeV
- **Allowed (Λ, m_X) region**: 3.16×10^8 < Λ ≤ 10^11 GeV, 10^6 ≤ m_X < 1.58×10^8 GeV — from a deterministic 221-point grid (not MCMC)
- **Multi-messenger constraints**: predicted Fermi-LAT, IceCube, and Pierre Auger fluxes are 2.6–24.8 orders of magnitude below current limits everywhere in the allowed region — current observations provide no data-driven constraint, as implemented in this pipeline
- **Cosmic abundance**: no quantitative Ω_X^AGN bound is derived in this analysis (see manuscript for why)

---

## Results Validation

All figures are regenerated from code on each run.  To validate:

```bash
python code/grmhd_acceleration.py        && echo "GRMHD: OK"
python code/eft_cross_sections.py        && echo "EFT: OK"
python code/pdf_integration.py           && echo "PDF: OK"
python code/observational_constraints.py && echo "Constraints: OK"
python code/monte_carlo_sampling.py --seed 42 && echo "MCMC: OK"
```

Figures verified present after a full run:

| Figure | File | Size |
|--------|------|------|
| 1 | `fig1_eft_validity_corrected.pdf` | 29 KB |
| 2 | `fig2_gamma_max.pdf` | 25 KB |
| 3 | `fig3_production_rate.pdf` | 22 KB |
| 4 | `fig4_constraints.pdf` | 29 KB |
| S1 | `figS1_mcmc_diagnostics.pdf` | 186 KB |

The reported inference is deterministic (see manuscript Section on Deterministic Parameter-Space Inference); `monte_carlo_sampling.py` is a legacy/inactive implementation not used to produce the reported result.

---

## Code Modules

### 1. GRMHD Acceleration (`grmhd_acceleration.py`)
Calculates maximum Lorentz factors from:
- Magnetic reconnection in sigma >> 1 jets (gamma_e = eta * sigma)
- Hillas criterion for protons (gamma_p = eBR / m_p c^2)
- **Generates Figure 2**

### 2. EFT Cross Sections (`eft_cross_sections.py`)
Computes production cross sections using:
- Dimension-6 contact interaction: sigma = s * beta^3 / (48 pi Lambda^4)
- EFT validity check: sqrt(s) < Lambda
- **Generates Figures 1 & 3**

### 3. PDF Integration (`pdf_integration.py`)
Integrates parton luminosities:
- Parametric/exploratory quark-antiquark luminosity stand-in: dL/dtau ~ tau^{-1} (1 - tau)^n (not CT18 NNLO; not benchmarked or validated against CT18 NNLO or any real PDF set)
- Hadronic cross section convolution at sqrt(s) ~ 10^8--10^9 GeV

### 4. Monte Carlo Sampling (`monte_carlo_sampling.py`)
Bayesian parameter inference:
- emcee affine-invariant sampler
- 10^3 samples, 10 walkers x 100 steps (legacy/inactive implementation, not used for the reported result)
- Random seeds: {42, 137, 271, 314, 628} for reproducibility
- **Generates Figure S1** (diagnostics)

### 5. Observational Constraints (`observational_constraints.py`)
Multi-messenger limits:
- Fermi-LAT: E^2 dN/dE < 10^{-12} erg/cm^2/s
- IceCube: phi_nu < 10^{-18} GeV/cm^2/s/sr per AGN
- Pierre Auger: UHECR spectrum above 10^{19.5} eV
- **Generates Figure 4** (illustrative/schematic curve, not a fitted model prediction)

---

## Figures

All figures are generated programmatically from code:

- **Figure 1**: EFT validity regions (corrected)
- **Figure 2**: Maximum Lorentz factors vs magnetization
- **Figure 3**: Production rates vs UHDM mass
- **Figure 4**: Multi-messenger constraints summary (illustrative/schematic curve only, not a calculated model prediction)
- **Figure S1**: MCMC convergence diagnostics

Figures are saved in `figures/` as vector PDFs (300 dpi).

---

## Reproducibility

### Computational Requirements
- **Deterministic (Λ, m_X) grid evaluation**: seconds on a single core
- **Legacy MCMC module** (`monte_carlo_sampling.py`, not used for the reported result): ~1 minute for its 10^3-sample default

### Random Seeds
All random number generators are seeded for bit-exact reproducibility:
```python
np.random.seed(42)  # Default seed (configurable via --seed flag)
```

### Docker Container
Fully containerized environment ensures platform-independent reproducibility:
```bash
docker run -it uhdm-agn bash workflow/run_analysis.sh
```

### Validation Checklist
Before submission, verify:
- [x] All figures regenerate from code
- [x] Manuscript compiles without errors
- [x] Data provenance documented
- [x] Random seeds set for reproducibility
- [x] Docker image builds successfully

See `docs/submission_checklist.txt` for full checklist.

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{patidar2025uhdm,
  author = {Patidar, Kirti},
  title = {Ultra-Heavy Dark Matter Production in Active Galactic Nucleus Jets:
           A Falsifiable Multi-Messenger Framework},
  year = {2025},
  note = {Preprint available at https://github.com/kirtii56/agn}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- Fermi-LAT Collaboration for gamma-ray data
- IceCube Collaboration for neutrino limits
- Pierre Auger Collaboration for UHECR data
- emcee developers for MCMC framework

---

## Contact

**Kirti Patidar**
Independent Researcher
Email: kirtiptdr174@gmail.com

For questions or issues:
- Open an issue on GitHub
- Email the author directly

---

## Documentation

For detailed documentation, see:
- [`docs/README.md`](docs/README.md) - Comprehensive technical documentation
- [`docs/metadata.json`](docs/metadata.json) - Project metadata for arXiv
- [`docs/submission_checklist.txt`](docs/submission_checklist.txt) - Pre-submission validation
- [`manuscript/`](manuscript/) - LaTeX source and compiled PDF

---

**Status**: In preparation for arXiv submission
**Version**: 2.1
**Last Updated**: February 2026
