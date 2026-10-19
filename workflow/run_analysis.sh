#!/bin/bash
#
# Master Analysis Execution Script
# =================================
#
# Runs the pipeline that actually produces the reported result:
# 1. GRMHD acceleration calculations
# 2. EFT cross section calculations
# 3. PDF integration
# 4. Observational constraints
# 5. Deterministic (log10 Lambda, log10 m_X) grid reproduction -- this IS the
#    reported scientific result
# 6. Figure generation
# 7. Manuscript compilation
#
# The legacy MCMC module (monte_carlo_sampling.py) is NOT run by this script.
# It is not used to produce the reported result; run it yourself only if you
# specifically want to exercise that legacy code path (see REPRODUCING.md).
#

set -e  # Exit on error
set -u  # Exit on undefined variable

echo "========================================="
echo "UHDM AGN Analysis Pipeline"
echo "========================================="
echo ""

# Navigate to code directory
cd "$(dirname "$0")/.."
CODE_DIR="code"
FIGURES_DIR="figures"
MANUSCRIPT_DIR="manuscript"

# Step 1: GRMHD Acceleration
echo "[1/7] Running GRMHD acceleration calculations..."
python ${CODE_DIR}/grmhd_acceleration.py
echo "✓ Complete"
echo ""

# Step 2: EFT Cross Sections
echo "[2/7] Calculating EFT cross sections..."
python ${CODE_DIR}/eft_cross_sections.py
echo "✓ Complete"
echo ""

# Step 3: PDF Integration
echo "[3/7] Integrating parton distribution functions..."
python ${CODE_DIR}/pdf_integration.py
echo "✓ Complete"
echo ""

# Step 4: Observational Constraints
echo "[4/7] Applying observational constraints..."
python ${CODE_DIR}/observational_constraints.py
echo "✓ Complete"
echo ""

# Step 5: Deterministic grid reproduction -- THE reported result
echo "[5/7] Reproducing the deterministic (log10 Lambda, log10 m_X) grid result..."
python ${CODE_DIR}/deterministic_grid_reproduction.py
echo "✓ Complete"
echo ""

# Step 6: Generate and verify figures
echo "[6/7] Generating figures..."
python ${CODE_DIR}/generate_figures.py --output ${FIGURES_DIR}/
REQUIRED_FIGURES=(
    "fig1_eft_validity_corrected.pdf"
    "fig2_gamma_max.pdf"
    "fig3_production_rate.pdf"
    "fig4_constraints.pdf"
)

for fig in "${REQUIRED_FIGURES[@]}"; do
    if [ -f "${FIGURES_DIR}/${fig}" ]; then
        echo "  ✓ ${fig}"
    else
        echo "  ✗ ${fig} - MISSING!"
    fi
done
echo ""

# Step 7: Compile LaTeX manuscript
echo "[7/7] Compiling LaTeX manuscript..."
cd ${MANUSCRIPT_DIR}

if command -v pdflatex &> /dev/null; then
    pdflatex manuscript_final.tex
    bibtex manuscript_final
    pdflatex manuscript_final.tex
    pdflatex manuscript_final.tex
    echo "✓ Manuscript compiled: manuscript_final.pdf"
else
    echo "⚠ pdflatex not found. Skipping manuscript compilation."
    echo "  Install with: apt-get install texlive-full"
fi

cd ..
echo ""

# Summary
echo "========================================="
echo "Analysis Pipeline Complete!"
echo "========================================="
echo ""
echo "Output files:"
echo "  Figures: ${FIGURES_DIR}/"
echo "  Manuscript: ${MANUSCRIPT_DIR}/manuscript_final.pdf"
echo ""
echo "Next steps:"
echo "  1. Review all figures in ${FIGURES_DIR}/"
echo "  2. Review manuscript PDF"
echo ""
echo "Not run by this script (legacy, not used for the reported result):"
echo "  python ${CODE_DIR}/monte_carlo_sampling.py"
echo ""
