MANUSCRIPT FILES
================

manuscript_final.tex - LaTeX source file in RevTeX4-2 format
references.bib - BibTeX bibliography

To compile:
  pdflatex manuscript_final.tex
  bibtex manuscript_final
  pdflatex manuscript_final.tex
  pdflatex manuscript_final.tex

NOTE: manuscript_final.pdf and supplementary_material.pdf should be generated
by compiling the LaTeX sources. These are binary PDF files not tracked in git
during development, but should be included for final submission.
