# Ultra-Heavy Dark Matter Production in Active Galactic Nucleus Jets

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the complete analysis pipeline for our study of ultra-high-density matter (UHDM) production in active galactic nucleus (AGN) jets using multi-messenger astronomy constraints.

**Key Results:**
- GRMHD simulations predict maximum Lorentz factors γ_max ~ 10^6-10^7 in magnetically-dominated jets
- EFT calculations compute production cross sections for a Dirac-fermion UHDM candidate X with masses explored over Λ ∈ [10^8, 10^11] GeV, m_X ∈ [10^6, 10^10] GeV (equivalently ~10^15-10^19 eV)
- Within the theoretically accepted region (3.16×10^8 < Λ ≤ 10^11 GeV, 10^6 ≤ m_X < 1.58×10^8 GeV), no point produces a nonzero observational penalty under the current pipeline's hardcoded Fermi-LAT, IceCube, and Pierre Auger limits — i.e., the current implementation finds no data-driven constraint on (Λ, m_X) over the explored, theoretically valid parameter space. This is a statement about the current pipeline's hardcoded limits, not an independently verified claim about present-day experimental sensitivity.
- The reported inference is a deterministic characterization of the (log10 Λ, log10 m_X) parameter space (a 221-point grid plus the closed-form accepted-region boundaries), not an MCMC posterior

## Repository Structure

```
UHDM_AGN_Paper_Final/
├── manuscript/          # LaTeX manuscript and bibliography
├── figures/            # Publication-quality PDF figures
├── code/               # Python analysis scripts
├── data/               # Observational data and MCMC outputs
├── workflow/           # Docker and automation scripts
├── docs/               # Documentation and metadata
└── LICENSE             # MIT License
```

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

# Monte Carlo scan (WARNING: slow!)
python code/monte_carlo_sampling.py

# Observational constraints
python code/observational_constraints.py
```

## Code Modules

### 1. GRMHD Acceleration (`grmhd_acceleration.py`)
Calculates maximum Lorentz factors from:
- Magnetic reconnection in σ >> 1 jets
- Relativistic shock acceleration
- Generates **Figure 2**

### 2. EFT Cross Sections (`eft_cross_sections.py`)
Computes production cross sections using:
- Contact interaction operators
- EFT validity analysis
- Generates **Figures 1 & 3**

### 3. PDF Integration (`pdf_integration.py`)
Integrates parton luminosities:
- Parametric/exploratory quark-antiquark PDF/luminosity stand-in (NOT CT18 NNLO; not benchmarked or validated against a real PDF set). A validated LHAPDF/CT18 NNLO implementation is deferred to Track B.
- Hadronic cross section convolution

### 4. Monte Carlo Sampling (`monte_carlo_sampling.py`) — LEGACY, NOT the reported inference
This module implements an emcee-based Bayesian sampler. It is a legacy/inactive
implementation and is **not** used to produce the reported result. The reported
inference is deterministic (see below).
- emcee affine-invariant sampler (legacy code path only)
- Generates **Figure S1** (diagnostics) — legacy artifact, not part of the reported result

### 4b. Deterministic Inference (reported result)
The reported inference characterizes the likelihood over `(log10 Λ, log10 m_X)`
directly, without sampling:
- Prior: Λ ∈ [10^8, 10^11] GeV, m_X ∈ [10^6, 10^10] GeV
- Fixed benchmark: sqrt(s_hat) = 3.16×10^8 GeV
- EFT validity: sqrt(s_hat) < Λ
- Production threshold: m_X < sqrt(s_hat)/2
- Accepted region: 3.16×10^8 < Λ ≤ 10^11 GeV, 10^6 ≤ m_X < 1.58×10^8 GeV
- Within this region the implemented likelihood is flat (constant) in the sampled
  log-coordinates; outside it the implementation returns a large negative penalty
  (-1e10), which is numerically negligible but not mathematically exact zero
- Characterized via a 221-point deterministic grid plus the closed-form region
  boundaries above; no MCMC sampler is required to reproduce this result

### 5. Observational Constraints (`observational_constraints.py`)
Multi-messenger limits:
- Fermi-LAT 3FHL catalog
- IceCube neutrino stacking
- Pierre Auger UHECR spectrum
- Generates **Figure 4**

## Data Requirements

Download observational data:

```bash
# Fermi-LAT 3FHL Catalog
wget https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3FHL/gll_psch_v13.fit \
  -O data/fermi_lat_3fhl_catalog.fits

# IceCube and Auger data
# See data/README_data.txt for details
```

Or use built-in parameterizations (less accurate).

## Reproducibility

### Dependencies
- Python 3.9+
- NumPy, SciPy, Matplotlib
- emcee, corner (MCMC)
- h5py, pandas, astropy
- Optional: LHAPDF (for CT18 PDFs)
- LaTeX (for manuscript compilation)

### Computational Requirements
- **Deterministic inference (reported result)**: the (log10 Λ, log10 m_X) grid
  characterization and the individual figure-generating scripts (Figures 1-4) run in
  well under a minute each; no runtime estimate is established by the repository for
  a single consolidated "run everything" command (see `REPRODUCING.md` for the current
  reproduction procedure and its limitations)
- **Legacy MCMC module** (`monte_carlo_sampling.py`, not used for the reported result):
  runtime depends on the `--nwalkers`/`--nsteps` values chosen at invocation; no runtime
  is established by the repository for any specific configuration

### Random Seeds
All random number generators are seeded for reproducibility:
```python
np.random.seed(42)
emcee uses deterministic initialization
```

### Docker Container
Fully containerized environment ensures bit-for-bit reproducibility:
```bash
docker run -it uhdm-agn bash workflow/run_analysis.sh
```

## Figures

All figures are generated programmatically from code:

- **Figure 1**: EFT validity regions (corrected)
- **Figure 2**: Maximum Lorentz factors vs magnetization
- **Figure 3**: Production rates vs UHDM mass
- **Figure 4**: Multi-messenger constraints summary
- **Figure S1**: MCMC convergence diagnostics (legacy artifact — not part of the reported deterministic inference)

Figures are saved in `figures/` as vector PDFs (300 dpi).

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

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.

## Acknowledgments

- Fermi-LAT Collaboration for gamma-ray data
- IceCube Collaboration for neutrino limits
- Pierre Auger Collaboration for UHECR data

## Contact

For questions or issues:
- Open an issue on GitHub
- Email: kirtiptdr174@gmail.com

## Validation Checklist

Before submission, verify:

- [ ] All figures regenerate from code
- [ ] Deterministic (log10 Λ, log10 m_X) grid characterization reproducible
- [ ] Manuscript compiles without errors
- [ ] Data provenance documented
- [ ] Random seeds set for reproducibility
- [ ] Docker image builds successfully

See `docs/submission_checklist.txt` for full checklist.

---

**Status**: Preprint-ready, not yet submitted to arXiv
**Version**: 2.0-FINAL
**Last Updated**: 2026-10-19
