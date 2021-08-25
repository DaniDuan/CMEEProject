#!/bin/bash
x=${1%.tex}
pdflatex $x.tex
pdflatex $x.tex
bibtex $x
pdflatex $x.tex
pdflatex $x.tex

## Cleanup
rm *~
rm *.aux
rm *.dvi
rm *.log
rm *.nav
rm *.out
rm *.snm
rm *.toc
rm *.bbl
rm *.blg
