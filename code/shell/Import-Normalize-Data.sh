#!/bin/bash

echo "Importing and merging count data..."
Rscript ./code/R/ImportCountData.R
echo "Finished importing and merging count data."
echo
echo

echo "Grouping and normalizing data..."
Rscript ./code/R/NormalizeCountData.R
echo "Finished grouping and normalizing data."
echo
echo

