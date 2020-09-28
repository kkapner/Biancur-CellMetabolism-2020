#!/bin/bash

source ~/.bash_profile

cd "$(dirname "${BASH_SOURCE[0]}")"
cd ..

echo
echo "Analyzing Biancur CRISPR Screen..."
echo 

./code/shell/Import-Normalize-Data.sh
./code/shell/QC.sh
source ./code/shell/STARS.sh
Rscript ./code/R/GenerateFigures.R
cd pipeline