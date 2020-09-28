Each folder is titled with the condition being held constant in the comparison. 

For example, 2D contains KO1_WT_2D_results.tsv and KO2_WT_2D_results.tsv.
KO1_WT_2D_results.tsv contains the results from the comparison between 2D_KO1 and 2D_WT. The sign of the log fold changes correspond to the the order listed in the file name (e.g. KO1_WT_2D_results.tsv is 2D_KO1 - 2D_WT).

For reference, a map from Ensembl IDs to gene symbols is provided (20200116_EnsemblToSymbols.csv).

Raw sample count data used for DESeq2 analysis:
	STAR-ReadCounts-Sample-Corrected.tsv

Sample FPKM values:
	FPKM-Sample-Corrected.tsv

In original data from Novogene, samples F_12D_1 and F_13D_3 were swapped. All count tables are corrected and analysis was performed on the corrected data.