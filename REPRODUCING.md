# Reproducing Results

**Version:** 2.0-FINAL
**Last Updated:** November 3, 2025

The reported scientific result is a **deterministic** characterization of the
`(log10 Lambda, log10 m_X)` parameter space. It does **not** use an MCMC posterior.
This document describes how to reproduce that reported result, and separately
describes the repository's legacy MCMC module, which exists in the codebase but is
**not** used to produce the reported result.

---

## Quick Start

Set up the environment and generate the current figures:

```bash
# 1. Set up environment
conda env create -f workflow/environment.yml
conda activate uhdm-agn

# 2. Generate Figures 1-4 (the figures used in the current pipeline)
python code/generate_figures.py --output figures/ --figures 1,2,3,4
```

This regenerates Figures 1-4 from the current code. It does **not** run the legacy
MCMC module (see "Legacy MCMC Module" below for why, and what it would do if run).

---

## Step-by-Step Reproduction

### 1. Environment Setup

**Option A: Conda (Recommended)**
```bash
conda env create -f workflow/environment.yml
conda activate uhdm-agn
```

**Option B: Docker**
```bash
docker build -t uhdm-agn -f workflow/Dockerfile .
docker run -it -v $(pwd):/workspace uhdm-agn
```

**Option C: pip (Manual)**
```bash
pip install -r code/requirements.txt
```

### 2. Data Acquisition (Optional)

Download observational data (or use the built-in hardcoded parameterizations already
present in `code/observational_constraints.py`):

```bash
# Fermi-LAT 3FHL Catalog (~50 MB)
wget https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3FHL/gll_psch_v13.fit \
  -O data/fermi_lat_3fhl_catalog.fits

# IceCube and Pierre Auger data
# See data/README_data.txt for details
```

### 3. Generate Figures

```bash
# Figure 1: EFT Validity boundary
python code/eft_cross_sections.py

# Figure 2: Maximum Lorentz Factors
python code/grmhd_acceleration.py

# Figure 3: Production Rates vs Mass
# (already generated as part of eft_cross_sections.py / pdf_integration.py)

# Figure 4: Multi-Messenger Constraints (illustrative curve; see caption in the
# generated PDF and figures/README_figures.txt)
python code/observational_constraints.py
```

Or generate all current figures at once:
```bash
python code/generate_figures.py --output figures/ --figures 1,2,3,4
```

**Figure S1 is a legacy artifact** (see below) and is intentionally excluded from
the commands above.

### 4. Reproduce the Deterministic Inference (reported result)

The reported inference characterizes the likelihood over `(log10 Lambda, log10 m_X)`
directly. **No single checked-in script currently performs this characterization
end-to-end** — this is a genuine, current reproducibility gap, stated here rather
than papered over. The characterization was performed using the following existing
functions, called directly (e.g. from a Python shell or a short ad hoc script), not
via a dedicated CLI entry point:

- `eft_dimension6_xsec` (`code/eft_cross_sections.py`) — the partonic cross section
  `sigma_hat = s_hat * beta**3 / (48*pi*Lambda**4)`, `beta = sqrt(1 - 4*m_X**2/s_hat)`.
- `hadronic_cross_section` (`code/pdf_integration.py`) — convolves the partonic
  cross section with the current parametric quark-antiquark PDF/luminosity stand-in
  (not CT18 NNLO; not benchmarked — see below).
- `log_likelihood` (`code/monte_carlo_sampling.py`) — evaluates the likelihood at a
  given `(log10 Lambda, log10 m_X, ...)` point. This function can be called directly,
  at any grid of points, **without invoking the MCMC sampler itself**
  (`sampler.run_mcmc(...)` is a separate call, never made for the reported result).

**Reported methodology, explicitly:**
- Parameters: `(log10 Lambda, log10 m_X)` only.
- Prior: `log10 Lambda` in `[8, 11]` (i.e. `Lambda` in `[1e8, 1e11]` GeV),
  `log10 m_X` in `[6, 10]` (i.e. `m_X` in `[1e6, 1e10]` GeV).
- Fixed benchmark: `sqrt(s_hat) = 3.16e8` GeV, used for every evaluation.
- EFT validity: `sqrt(s_hat) < Lambda`.
- Kinematic (pair-production) threshold: `m_X < sqrt(s_hat)/2`.
- Accepted region: `3.16e8 < Lambda <= 1e11` GeV, `1e6 <= m_X < 1.58e8` GeV.
- Within the accepted region, `log_likelihood` returns a constant value (the flat
  prior with no active Gaussian penalty from any of the three observational
  channels); outside it, `log_likelihood` returns `-1e10` (a large negative
  penalty — numerically negligible when exponentiated, but not mathematically an
  exact zero probability).
- The characterization used in the manuscript evaluates `log_likelihood` on a
  221-point grid spanning the prior box, plus the closed-form region boundaries
  above (derived directly from the fixed benchmark and the two conditions above,
  not fit or estimated).

**To reproduce a grid point yourself**, an example (illustrative, not a checked-in
script):
```python
import numpy as np
from monte_carlo_sampling import log_likelihood

for log_Lambda in np.arange(8.0, 11.01, 0.25):
    for log_mX in np.arange(6.0, 10.01, 0.25):
        # sigma_jet, B_jet are not active inference parameters; any values in their
        # prior range give the same result, since they do not enter log_likelihood's
        # calculation of N_dot/flux (see code/monte_carlo_sampling.py)
        ll = log_likelihood((log_Lambda, log_mX, 1.5, 1.0))
        print(log_Lambda, log_mX, ll)
```

### 5. Compile Manuscript

```bash
cd manuscript/
latexmk -pdf manuscript_final.tex
# or manually:
pdflatex manuscript_final.tex
bibtex manuscript_final
pdflatex manuscript_final.tex
pdflatex manuscript_final.tex
```

**Output:** `manuscript_final.pdf`

---

## Legacy MCMC Module (`monte_carlo_sampling.py`)

This module implements an `emcee`-based Bayesian sampler with parameters
`(log10 Lambda, log10 m_X, log10 sigma_jet, log10 B_jet)`. **It is not used to
produce the reported result.** `sigma_jet` and `B_jet` are computed into standalone
Hillas/reconnection sanity-check values inside `log_likelihood` and do not affect
`N_dot` or the observational-flux comparison — they are not active inference
parameters for the reported deterministic result.

If you choose to run this legacy module anyway (e.g. to explore it, not to
reproduce the reported result):
```bash
python code/monte_carlo_sampling.py
```
No runtime, sample-count, chain-count, or convergence-diagnostic figure is
established by this repository for any specific configuration of this module —
its default settings and any timing depend on the flags you pass; do not assume a
particular runtime without measuring it yourself. Any output figure from this
module (`figS1_mcmc_diagnostics.pdf`) is a legacy artifact, not part of the
reported result.

---

## PDF Treatment — Current Limitation

The current `pdf_integration.py` implements an **exploratory parametric
quark-antiquark PDF/luminosity stand-in** (`dL/dtau ~ tau^a (1-tau)^n`). It is
**not** CT18 NNLO and has **not** been benchmarked or validated against a real PDF
set. It has no factorization/renormalization scale (`Q^2`) dependence. A validated
LHAPDF/CT18 NNLO implementation, with an explicit scale choice, is deferred to
future (Track B) work and is not part of the current reproduction procedure.

---

## Software Dependencies

- Python 3.9+
- NumPy 1.21+, SciPy 1.7+, Matplotlib 3.4+
- emcee 3.1+, corner 2.2+ (required only for the legacy MCMC module; not needed to
  reproduce the reported deterministic result, though currently still listed as a
  core dependency in `code/requirements.txt`)
- h5py 3.3+, pandas 1.3+, astropy 4.3+
- LaTeX (TeXLive 2020+ or MiKTeX), for manuscript compilation
- Optional: LHAPDF 6.5+ (Track B only; not used by the current PDF stand-in)

---

## Troubleshooting

**1. LaTeX Compilation Errors**
```bash
sudo apt-get install texlive-full  # Ubuntu/Debian
brew install --cask mactex          # macOS
```

**2. Memory Errors in PDF Integration**
```bash
# Reduce batch size in PDF integration
# Edit code/pdf_integration.py to use smaller arrays
```

**3. Figure Generation Fails**
```bash
python -c "import matplotlib; print(matplotlib.get_backend())"
export MPLBACKEND=Agg
```

---

## Output Files (current pipeline)

```
figures/
├── fig1_eft_validity_corrected.pdf
├── fig2_gamma_max.pdf
├── fig3_production_rate.pdf
└── fig4_constraints.pdf

manuscript/
└── manuscript_final.pdf
```

`figures/figS1_mcmc_diagnostics.pdf` is a legacy artifact from the (not-used)
MCMC module; it is not part of the current pipeline's output list above unless you
choose to run that legacy module yourself.

---

## Support

If you encounter issues:

1. Check `docs/submission_checklist.txt`
2. Review error messages in terminal output
3. Open an issue on GitHub: https://github.com/kirtii56/agn/issues
4. Email: kirtiptdr174@gmail.com

---

## References

For methodology details, see the manuscript:
- **Section 2**: AGN Jet Acceleration Physics
- **Section 3**: UHDM Production Framework (operator, cross section, PDF treatment)
- Deterministic Parameter-Space Inference section (the reported inference method)
- **Section 4**: Multi-Messenger Constraints

For code documentation:
- `docs/README.md` - detailed technical documentation
- Inline comments and docstrings in the Python scripts in `code/`

---

**Last Updated:** November 3, 2025
**Version:** 2.0-FINAL
