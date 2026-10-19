# Methods Documentation

Technical notes on the physics models and parameterisations used in the
UHDM AGN production analysis.  These notes supplement the manuscript and
explain choices made in the code where the full derivation or data source
is not included in the repository.

---

## 1. EFT Cross Section Formula

### Lagrangian

We use a dimension-6 contact interaction between Standard Model quarks
and the UHDM Dirac fermion X:

    L_eff = (1 / Lambda^2) (q-bar q)(X-bar X)

where Lambda is the EFT cutoff scale.

### Partonic cross section

The tree-level q qbar -> X Xbar cross section from this operator is:

    sigma_hat(s) = (s / Lambda^4) * (1 / (48 pi)) * beta(s, m_X)**3

where

    beta = sqrt(1 - 4 m_X^2 / s)

is the final-state velocity factor, and s is the partonic centre-of-mass
energy squared.

### Derivation sketch

For massless initial quarks q(p1), qbar(p2) and Dirac-fermion final
states X(p3), Xbar(p4), the amplitude from the contact vertex is

    M = (1 / Lambda^2) * [vbar(p2) u(p1)] * [ubar(p3) v(p4)]

Summing the standard fermion completeness relations over final spins
and initial spins gives

    Sum_spins |vbar(p2) u(p1)|^2 = 4 (p1.p2)          (m_q = 0)
    Sum_spins |ubar(p3) v(p4)|^2 = 4 (p3.p4 - m_X^2)

Using p1.p2 = s/2 and p3.p4 - m_X^2 = s*beta^2/2 (with
beta = sqrt(1 - 4 m_X^2/s)), and averaging over the 4 initial spin
states (1/4) and the 9 initial q-qbar color combinations, of which
only the 3 color-singlet-diagonal ones are non-zero (color average
1/3), the spin- and color-averaged squared matrix element is

    |Mbar|^2 = (1/4) * (1/3) * (1/Lambda^4) * [4 (s/2)] * [4 (s*beta^2/2)]
             = s^2 * beta^2 / (3 * Lambda^4)

Combining this with the standard massless-initial-state 2-body phase
space and flux factor, sigma_hat = |Mbar|^2 * beta / (16 pi s), gives

    sigma_hat(s) = (s / Lambda^4) * (1 / (48 pi)) * beta^3

This result has been independently verified by explicit Dirac-spinor
numerical evaluation and by an independent helicity/decay-width
cross-check (see audit_report/ in the repository history for details).

### EFT validity

The effective description is valid only when the partonic energy is below
the cutoff:

    sqrt(s) < Lambda

Above this scale, the contact interaction must be UV-completed by a
mediator, and the cross section formula is no longer reliable.

### Units

The code returns cross sections in cm^2 using the conversion
1 GeV^{-2} = 3.894 x 10^{-28} cm^2.

### Hadronic convolution / PDF treatment

The partonic cross section above is convolved with a quark-antiquark
parton luminosity to obtain the hadronic cross section used downstream.
The current implementation (`code/pdf_integration.py`) uses an
**exploratory parametric stand-in**, `dL/dtau ~ tau^a (1 - tau)^n`. This
is **not** CT18 NNLO and has **not** been benchmarked or validated
against a real CT18 NNLO grid or any other real PDF set; it has no
factorization/renormalization-scale (`Q^2`) dependence. A validated
LHAPDF/CT18 NNLO implementation, with an explicit scale choice, is
deferred to future (Track B) work.

---

## 2. GRMHD Parameter Choices

### Magnetic field strength:  B = 0.1 -- 100 G

This range corresponds to AGN jet environments at distances ~0.1--10 pc
from the central engine.  Observational constraints come from:

- Faraday rotation measure gradients across jets
  (Zavala & Taylor 2005, ApJ 626, L73)
- Jet power arguments equating Poynting flux to observed luminosity
  (Zamaninasab et al. 2014, Nature 510, 126)
- Core-shift measurements in VLBI-resolved jets
  (Lobanov 1998, A&A 330, 79)

Values above ~100 G are relevant only very close to the black hole
(< 10 R_g) and are not representative of the extended acceleration zone.

### Magnetization parameter:  sigma = 10 -- 1000

GRMHD simulations of magnetically arrested disks (MAD) produce jets with
sigma >> 1 near the jet base, decreasing with distance.  Key references:

- Tchekhovskoy et al. 2011, MNRAS 418, L79
  (MAD jets with sigma ~ 10--100 at launch)
- McKinney et al. 2012, MNRAS 423, 3083
  (sigma up to ~1000 in the jet spine)

The sigma >> 1 regime is where magnetic reconnection dominates particle
acceleration over shock acceleration.

### Jet radius:  R_jet = 10^14 -- 10^16 cm

This covers the range from ~0.003 pc to ~0.3 pc:

- 10^14 cm:  inner jet / blazar-zone scale
- 10^15 cm:  typical reconnection / acceleration region
- 10^16 cm:  transition to large-scale jet (~pc)

These scales are set by the jet opening angle (~1/Gamma_bulk) and the
distance from the central engine.

### Electron acceleration:  gamma_e_max = eta * sigma

In the sigma >> 1 regime, PIC simulations show that magnetic reconnection
accelerates electrons to Lorentz factors proportional to the magnetization:

    gamma_e_max = eta * sigma

with eta ~ 0.1--0.3 (Sironi & Spitkovsky 2014, ApJ 783, L21; Werner et
al. 2016, ApJ 816, L8).  The efficiency eta encapsulates the fraction of
magnetic energy converted to the highest-energy particles.

### Proton acceleration:  gamma_p_max = e B R_jet / (m_p c^2)

The Hillas criterion (Hillas 1984, ARA&A 22, 425) gives the maximum
energy a proton can reach while its Larmor radius fits inside the
acceleration region:

    E_max = e B R_jet
    gamma_p_max = E_max / (m_p c^2)

For B = 10 G and R_jet = 10^15 cm this gives gamma_p ~ 3 x 10^9,
consistent with UHECR energies.

---

## 3. Observational Limit Parameterisations

Since the full Fermi-LAT, IceCube, and Auger data files are not included
in this repository, the code uses parameterised upper limits drawn from
published results.  These are sufficient for the order-of-magnitude
predictions in this analysis.

### Fermi-LAT (3FHL catalog)

- **Observable**: point-source sensitivity in E^2 dN/dE
- **Limit used**: E^2 dN/dE < 10^{-12} erg cm^{-2} s^{-1}
- **Energy range**: 10 GeV -- 2 TeV
- **Basis**: Ajello et al. 2017, ApJS 232, 18 (3FHL catalog)
- **Notes**: The 3FHL sensitivity varies across the sky; we use the
  median sensitivity as a representative upper limit.  For specific AGN,
  one should use the per-source upper limit from the catalog.

### IceCube (AGN neutrino stacking)

- **Observable**: per-AGN muon-neutrino flux
- **Limit used**: phi_nu < 10^{-18} GeV cm^{-2} s^{-1} sr^{-1}
- **Energy range**: 1 TeV -- 10 PeV
- **Basis**: IceCube Collaboration (Abbasi et al.) 2023, ApJ 954, 75
- **Notes**: This is the stacking limit from ~1000 gamma-ray-selected
  AGN.  Individual bright sources (e.g. NGC 1068) have tighter
  constraints.  The limit assumes an E^{-2} neutrino spectrum.

### Pierre Auger (UHECR spectrum)

- **Observable**: integral cosmic-ray flux above threshold energy
- **Threshold**: E > 10^{19.5} eV (31.6 EeV)
- **Basis**: Pierre Auger Collaboration 2020, PRL 125, 121106
- **Notes**: The Auger constraint is parameterised as a simple integral
  flux limit above the energy threshold.  In practice, the composition-
  sensitive observables (X_max distributions) provide additional
  constraints on UHDM decay products.  The code uses an approximate
  limit that is conservative (i.e. weaker than the actual data).

---

## 4. Deterministic Parameter-Space Inference (reported methodology)

The reported inference characterizes the likelihood over `(log10 Lambda,
log10 m_X)` directly, without stochastic sampling.

### Parameters and prior

Only `Lambda` and `m_X` are inference parameters. `sigma_jet` and
`B_jet` (see the legacy MCMC subsection below) are **not** active
inference parameters and do not affect the likelihood computed here.

| Parameter | Symbol | Prior range |
|-----------|--------|-------------|
| EFT cutoff | log10(Lambda/GeV) | [8, 11] |
| UHDM mass | log10(m_X/GeV) | [6, 10] |

### Fixed benchmark and validity conditions

A single fixed benchmark collision energy is used for every evaluation:

    sqrt(s_hat) = 3.16 x 10^8 GeV

Two distinct conditions determine whether a given `(Lambda, m_X)` point
is theoretically accepted:

- **EFT validity**: sqrt(s_hat) < Lambda
- **Kinematic (pair-production) threshold**: m_X < sqrt(s_hat)/2

These are physically distinct: the first is a statement about the
self-consistency of the effective-theory expansion; the second is a
statement about whether the fixed benchmark has enough energy to
produce an X Xbar pair at all. Both must hold for a point to be
theoretically accepted.

### Accepted region

    3.16 x 10^8 < Lambda <= 10^11 GeV
    10^6 <= m_X < 1.58 x 10^8 GeV

(The Lambda bound is strict because the code requires
`sqrt(s_hat) < Lambda` strictly. The m_X bound is strict because,
exactly at m_X = sqrt(s_hat)/2, beta = 0 and the cross section is
exactly zero, which the likelihood's `xsec <= 0` check rejects.)

### Likelihood

Within the accepted region, the likelihood is a one-sided Gaussian
penalty:

    log L = -0.5 * sum_i [ max(0, F_pred_i - F_lim_i)^2 / sigma_i^2 ]

where the sum runs over the three messenger channels (gamma-ray,
neutrino, cosmic ray). Outside the accepted region, the implementation
returns a fixed numerical penalty, `log L = -1e10`. This is a large
negative number chosen to make invalid points strongly disfavored; it
is **not** a claim of mathematically exact zero probability — the
corresponding unnormalized density `exp(-1e10)` is a specific, nonzero
positive real number that simply underflows to numerically negligible
values in any practical computation.

### Deterministic characterization

The reported result was obtained by evaluating the likelihood above on
a 221-point grid spanning the prior box (0.25-dex spacing in each of
`log10 Lambda` and `log10 m_X`), together with the closed-form
accepted-region boundaries derived directly from the two validity
conditions above (not fit or estimated from the grid). Across this
grid, within the accepted region, the likelihood is flat (constant) in
the sampled log-coordinates — no grid point produces a nonzero
Gaussian penalty from any of the three observational channels under
the current pipeline's hardcoded limits. This is a statement about the
current pipeline as implemented, not an independently verified claim
about present-day experimental sensitivity.

No MCMC sampler is required to obtain or reproduce this result.

## 5. Legacy MCMC Implementation (`monte_carlo_sampling.py`)

**This subsection describes legacy code. It is not the reported
methodology and was not used to produce the reported result.**

The repository also contains an `emcee`-based Bayesian sampler over a
4-parameter vector `(log10 Lambda, log10 m_X, log10 sigma_jet,
log10 B_jet)`. `sigma_jet` and `B_jet`, though present in this legacy
module's sampled vector, are only used inside standalone Hillas/
reconnection sanity-check calculations and do not affect the predicted
flux or the likelihood value computed by this module — they are not
active inference parameters even within the legacy module itself.

### Sampler

The emcee affine-invariant ensemble sampler (Foreman-Mackey et al.
2013, PASP 125, 306) is used, requiring nwalkers >= 2 * ndim; with
ndim = 4 the module's default is nwalkers = 10.

### Likelihood (identical structure to Section 4 above)

    log L = -0.5 * sum_i [ max(0, F_pred_i - F_lim_i)^2 / sigma_i^2 ]

with sigma_i = 0.5 * F_lim_i (a 50% systematic uncertainty on each
limit).

### Convergence diagnostics — not established for any specific reported run

No specific run of this legacy module (chain count, sample count, or
convergence outcome) is used for or claimed as the reported result.
Diagnostic machinery exists in the module (acceptance fraction,
`emcee.autocorr.integrated_time`, a split Gelman-Rubin `R-hat`
calculation), but no particular values from running it are reported
here. If this module is run, its own default configuration is small
(`nwalkers=10, nsteps=100`); larger configurations have not been
run and characterized in this repository, and no specific runtime,
sample count, or convergence diagnostic value should be assumed
without actually running and measuring it.
