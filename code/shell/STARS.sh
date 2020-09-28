#!/bin/bash

conda activate STARS

# Running STARS at a threshold of 60 
# Checks to see if Null distribution with a threshold of 60 was created and if it wasn't then it creates it

if ! [[ -f "output/stars/Null_STARSOutput8_60.txt" ]]; then
    echo "Creating null distribution..."
    # Creating Null distribution with a threshold of 60
    python code/STARS/stars_null_v1.3.py --input-file data/meta/Pertubations.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --num-ite 1000 --use-first-pert N

    echo "Finished creating null distribution."
    mkdir -p output/stars
    mv Null_STARSOutput8_60.txt output/stars
    
fi

if ! [[ -d "output/stars/individual/neg" ]]; then
    echo "Running STARS on individual data..."
    # Positive Direction (Individual)
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/LogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir P --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    # Negative Direction (Individual)
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/LogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir N --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    mkdir -p output/stars/individual/neg
    mkdir -p output/stars/individual/pos
    mv *_N_*.txt output/stars/individual/neg
    mv *_P_*.txt output/stars/individual/pos
    echo "Finished running STARS on individual data."
fi

if ! [[ -d "output/stars/comparison/p3/neg" ]]; then
    echo "Running STARS on P3 Comparison data..."
    # Positive Direction (P3 Comparison)
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/P3LogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir P --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    # Negative Direction (P3 Comparison)
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/P3LogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir N --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    mkdir -p output/stars/comparison/p3/neg
    mkdir -p output/stars/comparison/p3/pos
    mv *_N_*.txt output/stars/comparison/p3/neg
    mv *_P_*.txt output/stars/comparison/p3/pos
    echo "Finished running STARS on P3 Comparison data."
fi

if ! [[ -d "output/stars/comparison/p7/neg" ]]; then
    echo "Running STARS on P7 Comparison data..."
    # Positive Direction (P7 Comparison)
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/P7LogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir P --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    # Negative Direction (P7 Comparison)
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/P7LogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir N --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    mkdir -p output/stars/comparison/p7/neg
    mkdir -p output/stars/comparison/p7/pos
    mv *_N_*.txt output/stars/comparison/p7/neg
    mv *_P_*.txt output/stars/comparison/p7/pos
    echo "Finished running STARS on P7 Comparison data..."
fi

if ! [[ -d "output/stars/comparison/nude/neg" ]]; then
    echo "Running STARS on Nude Comparison data..."
    # Positive Direction
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/NudeLogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir P --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    # Negative Direction
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/NudeLogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir N --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    mkdir -p output/stars/comparison/nude/neg
    mkdir -p output/stars/comparison/nude/pos
    mv *_N_*.txt output/stars/comparison/nude/neg
    mv *_P_*.txt output/stars/comparison/nude/pos
    echo "Finished running STARS on Nude Comparison data..."

    # Nude comparison to B6 does not input p-value meta data to output file
    # Extracting header from P3 file
    head -n 1 output/stars/comparison/nude/neg/P3_*.txt > prefix.txt

    # Writing to file for negative and positive directions
    cat prefix.txt output/stars/comparison/nude/neg/B6_*.txt > data.new && mv data.new output/stars/comparison/nude/neg/B6_*.txt

    cat prefix.txt output/stars/comparison/nude/pos/B6_*.txt > data.new && mv data.new output/stars/comparison/nude/pos/B6_*.txt

    rm prefix.txt


fi

if ! [[ -d "output/stars/comparison/3D/neg" ]]; then
    echo "Running STARS on 3D Comparison data..."
    # Positive Direction (P3 Comparison)
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/3DLogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir P --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    # Negative Direction (P3 Comparison)
    python code/STARS/stars_v1.3.py --input-file output/normalized-data/3DLogFoldChangeData.txt --chip-file data/meta/Pertubation-gene_table.txt --thr 60 --dir N --null output/stars/Null_STARSOutput8_60.txt --use-first-pert N

    mkdir -p output/stars/comparison/3D/neg
    mkdir -p output/stars/comparison/3D/pos
    mv *_N_*.txt output/stars/comparison/3D/neg
    mv *_P_*.txt output/stars/comparison/3D/pos
    echo "Finished running STARS on 3D Comparison data."
fi

conda deactivate

#######################################################################
# Merging positive and negative STARS data and adding median LFC values
#######################################################################

echo "Cleaning up STARS data..."
Rscript code/R/CleanSTARSResults.R
echo "Finished cleaning up STARS data."

echo "Merging STARS data..."
conda activate BiancurCollab

python code/py/MergeSTARS.py output/stars/individual/ output/normalized-data/LogFoldChangeData-Gene-Labeled.txt
python code/py/MergeSTARS.py output/stars/comparison/p3 output/normalized-data/P3LogFoldChangeData-GeneLabeled.txt
python code/py/MergeSTARS.py output/stars/comparison/p7 output/normalized-data/P7LogFoldChangeData-GeneLabeled.txt
python code/py/MergeSTARS.py output/stars/comparison/nude output/normalized-data/NudeLogFoldChangeData-GeneLabeled.txt
python code/py/MergeSTARS.py output/stars/comparison/3D output/normalized-data/3DLogFoldChangeData-GeneLabeled.txt

conda deactivate
echo "Finished merging STARS data."
