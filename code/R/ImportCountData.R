source("./code/R/helpers/DataInput.R")

# Using the count data which is pre-log normalized as returned by the Broad
screen.data.paths <- c("./data/readcounts/lognorm/B6_lognorm-FP_GPP414_Aguirre_2018_07_25_P1_SP20180806.txt", 
                    "./data/readcounts/lognorm/B6_lognorm-FP_GPP414_Aguirre_2018_07_25_P2_SP20180806.txt",
                    "./data/readcounts/lognorm/Nude_lognorm-FP_GPP526_20181130_Aguirre_sgRNA_DB1.txt",
                    "./data/readcounts/lognorm/lognorm-FP_GPP817_Biancur_sgRNA_20190813_reseq.txt")

# Merging screen data from multiple data tables
merged.screen.data <- multi_plate_read(screen.data.paths)

# Converting the barcodes to corresponding gene name
gene.labeled.data <- barcode_to_genes(merged.screen.data, "./data/meta/Pertubation-gene_table.xlsx")

# Assigning negative control guide names
gene.labeled.data$GeneID[is.na(gene.labeled.data$GeneID)] <- "Neg_Control"

# Creating merged data table
if (!dir.exists("./data/readcounts/merged")){
    dir.create("./data/readcounts/merged", recursive = TRUE)
}

write.table(gene.labeled.data, "./data/readcounts/merged/merged-lognorm-sgRNA-B6P3P7N3D.txt", sep = "\t",
row.names = FALSE)