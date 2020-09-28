#!/bin/bash

echo "Performing QC..."

Rscript ./code/R/QC.R
mv Rplots.pdf ./output/figures/qc/

echo "Finished performing QC."
echo
echo