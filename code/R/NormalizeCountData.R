source("./code/R/helpers/DataGrouping.R")
source("./code/R/helpers/DataProcess.R")

# Creating globally used varaibles
gene.labeled.data <- read.table("./data/readcounts/merged/merged-lognorm-sgRNA-B6P3P7N3D.txt", header = TRUE)

positive.controls <- c("Rpl10", "Rpl13a", "Rpl14", "Rpl17", "Rpl23",
                        "Rpl27", "Rpl35a", "Rpl7", "Rplp2", "Rps11",
                       "Polr2d", "U2af2", "Eef2", "Psmd1", "Ran",
                       "Tubb5", "Psmb2", "Rptor", "Polr2i", "Psmd6",
                       "Supt5", "Psmc4")

##############################################################################
# Normalizatoin procedure which maintains the biolgoical replicates separatation when calculating the log fold chnage

sample.level.lfc <- data.frame(matrix(NA, nrow = nrow(gene.labeled.data), ncol = 1))[-1]

# Subtracting screen specific full representation from guides
sample.level.lfc[, "GeneID"] <- gene.labeled.data[, "GeneID"]
sample.level.lfc[, c("P7_1", "P7_2", "P7_3")] <- gene.labeled.data[, c("Passage7.1", "Passage7.2", "Passage7.3")] - gene.labeled.data[, c("FR.1", "FR.2", "FR.3")]
sample.level.lfc[, c("P3_1", "P3_2", "P3_3")] <- gene.labeled.data[, c("Passage3.1", "Passage3.2", "Passage3.3")] - gene.labeled.data[, c("FR.1", "FR.2", "FR.3")]
sample.level.lfc[, c("B6_R1C1", "B6_R1C2", "B6_R1C3", "B6_R1C4")] <- gene.labeled.data[, c("R1C1", "R1C2", "R1C3", "R1C4")] - gene.labeled.data[, c("FR.1", "FR.1", "FR.1", "FR.1")]
sample.level.lfc[, c("B6_R2C1", "B6_R2C2", "B6_R2C3", "B6_R2C4")] <- gene.labeled.data[, c("R2C1", "R2C2", "R2C3", "R2C4")] - gene.labeled.data[, c("FR.2", "FR.2", "FR.2", "FR.2")]
sample.level.lfc[, c("B6_R3C1", "B6_R3C2", "B6_R3C3", "B6_R3C4")] <- gene.labeled.data[, c("R3C1", "R3C2", "R3C3", "R3C4")] - gene.labeled.data[, c("FR.3", "FR.3", "FR.3", "FR.3")]
sample.level.lfc[, c("Nude_R1C1", "Nude_R1C2", "Nude_R1C3", "Nude_R1C4")] <- gene.labeled.data[, c("R1C1N", "R1C2N", "R1C3N", "R1C4N")] - gene.labeled.data[, c("FR1.1.x", "FR1.1.x", "FR1.1.x", "FR1.1.x")]
sample.level.lfc[, c("Nude_R2C1", "Nude_R2C2", "Nude_R2C3", "Nude_R2C4")] <- gene.labeled.data[, c("R2C1N", "R2C2N", "R2C3N", "R2C4N")] - gene.labeled.data[, c("FR1.2.x", "FR1.2.x", "FR1.2.x", "FR1.2.x")]
sample.level.lfc[, c("Nude_R3C1", "Nude_R3C2", "Nude_R3C3", "Nude_R3C4")] <- gene.labeled.data[, c("R3C1N", "R3C2N", "R3C3N", "R3C4N")] - gene.labeled.data[, c("FR1.3.x", "FR1.3.x", "FR1.3.x", "FR1.3.x")]
sample.level.lfc[, c("3D_1", "3D_2", "3D_3")] <- gene.labeled.data[, c("X3D_1.1", "X3D_1.2", "X3D_1.3")] - gene.labeled.data[, c("FR1.1.y", "FR1.2.y", "FR1.3.y")]

# Median collapsing replicates
sample.group.regex <- c("P7", "P3", "B6",
                        "Nude", "3D")
collapsed.sample.replicates <- group_by_sample(sample.level.lfc, sample.group.regex)
colnames(collapsed.sample.replicates) <- c("GeneID", "P7", "P3", "B6", "Nude", "ThreeD")

# Median collapsing to gene level
median.collapsed.sample.to.genes <- combine_by_id(collapsed.sample.replicates, "GeneID")

# Normalizing by subtracting median of negative controls and dividing by magnitude of positive controls. We also want to create normalized data which is not collapsed down to the gene level so that it can be run through STARS
normalized.sample.data.gene.level <- normalize_data(median.collapsed.sample.to.genes, positive.controls,
                                  pos_normalization = TRUE, collapsed = TRUE)

normalized.sample.data.guide.level <- normalize_data(collapsed.sample.replicates, positive.controls,
                                  pos_normalization = TRUE, collapsed = FALSE)

sample.data.replicate.gene.level <- combine_by_id(sample.level.lfc, "GeneID")
normalized.sample.data.replicate.gene.level <- normalize_data(sample.data.replicate.gene.level, positive.controls, pos_normalization = TRUE,
collapsed = FALSE)

if (!dir.exists("./output/normalized-data")){
  dir.create("./output/normalized-data", recursive = TRUE)
}

write.table(normalized.sample.data.gene.level, "./output/normalized-data/Replicate-Separate-Gene-Level-Normalized.txt", row.names = FALSE, sep = "\t")
write.table(normalized.sample.data.guide.level, "./output/normalized-data/Replicate-Separate-Guide-Level-Normalized.txt", row.names = FALSE, sep = "\t")
write.table(normalized.sample.data.replicate.gene.level, "./output/normalized-data/Normalized-Replicates-Gene-Level.txt", row.names = FALSE, sep = "\t")

##############################################################################

##############################################################################
# Normalizatoin procedure which median collapses the replicates and then proceeds with the normalization procedure

sample.group.regex <- c("FR.\\d$", "Passage7", "Passage3",
                        "R\\dC\\d$", "FR\\d.\\d.x$", "R\\dC\\dN",
                        "FR\\d.\\d.y$", "3D")
collapsed.to.replicates <- group_by_sample(gene.labeled.data,
                                                      sample.group.regex)
colnames(collapsed.to.replicates) <- c("GeneID", "FR_2DB6", "P7", 
                                       "P3", "B6", "FR_Nude", "Nude",
                                       "FR_3D", "ThreeD")

log.fold.change.gene.labeled <- data.frame(matrix(
  NA, nrow = nrow(collapsed.to.replicates), ncol = 1))[-1]
log.fold.change.gene.labeled[, "GeneID"] <- collapsed.to.replicates[, "GeneID"]
log.fold.change.gene.labeled[, "P7"] <- collapsed.to.replicates[, "P7"] -
  collapsed.to.replicates[, "FR_2DB6"]
log.fold.change.gene.labeled[, "P3"] <- collapsed.to.replicates[, "P3"]-
  collapsed.to.replicates[, "FR_2DB6"]
log.fold.change.gene.labeled[, "B6"] <- collapsed.to.replicates[, "B6"] -
  collapsed.to.replicates[, "FR_2DB6"]
log.fold.change.gene.labeled[, "Nude"] <- collapsed.to.replicates[, "Nude"] -
  collapsed.to.replicates[, "FR_Nude"]
log.fold.change.gene.labeled[, "ThreeD"] <- collapsed.to.replicates[, "ThreeD"] -
  collapsed.to.replicates[, "FR_3D"]

median.collapsed.to.genes <- combine_by_id(log.fold.change.gene.labeled, "GeneID")
normalized.data <- normalize_data(median.collapsed.to.genes, positive.controls,
                                  pos_normalization = TRUE, collapsed = TRUE)
normalized.lfc.data <- normalize_data(log.fold.change.gene.labeled, positive.controls,
                                      pos_normalization = TRUE, collapsed = FALSE)

write.table(collapsed.to.replicates, "./output/normalized-data/Collapsed-Replicates.txt", sep = "\t", row.names = FALSE)
write.table(normalized.data, "./output/normalized-data/NormalizedData.txt", sep = "\t", row.names = FALSE)
write.table(normalized.lfc.data, "./output/normalized-data/LogFoldChangeData-Gene-Labeled.txt", row.names = FALSE, sep = "\t")

lfc.no.genes <- normalized.lfc.data[, c("P7", "P3", "B6", "Nude", "ThreeD")]
lfc.no.genes <- cbind(gene.labeled.data[, "Construct.IDs"], lfc.no.genes)
colnames(lfc.no.genes) <- c("Barcodes", "P7", "P3", "B6", "Nude", "ThreeD")

write.table(lfc.no.genes, file = "./output/normalized-data/LogFoldChangeData.txt", sep = "\t", row.names = FALSE)

##############################################################################

# Creating data sets of expression relative to P3, P7, Nude, and 3D

# P3 
p3.compared.lfc <- normalized.lfc.data[, c("P7", "P3", "B6", "Nude", "ThreeD")] - 
  normalized.lfc.data[, "P3"]

# Adding in geneID columns
p3.compared.lfc <- cbind(gene.labeled.data[, "Construct.IDs"], p3.compared.lfc)

# Removing P3 from data (it is all zeros) and renaming columns
p3.compared.lfc <- within(p3.compared.lfc, rm("P3"))
colnames(p3.compared.lfc) <- c("Barcodes","P7", "B6", "Nude", "ThreeD")

write.table(p3.compared.lfc, file = "./output/normalized-data/P3LogFoldChangeData.txt", sep = "\t", row.names = FALSE)

p3.compared.lfc.gene.labeled <- cbind(gene.labeled.data[, "GeneID"], p3.compared.lfc)
colnames(p3.compared.lfc.gene.labeled) <- c("GeneID", "Barcodes","P7", "B6", "Nude", "ThreeD")
write.table(p3.compared.lfc.gene.labeled, file = "./output/normalized-data/P3LogFoldChangeData-GeneLabeled.txt", row.names = FALSE, sep = "\t")

# P7

p7.compared.lfc <- normalized.lfc.data[, c("P7", "P3", "B6", "Nude", "ThreeD")] - 
  normalized.lfc.data[, "P7"]

# Adding in geneID columns
p7.compared.lfc <- cbind(gene.labeled.data[, "Construct.IDs"], p7.compared.lfc)
# Removing P3 from data (it is all zeros) and renaming columns
p7.compared.lfc <- within(p7.compared.lfc, rm("P7"))
colnames(p7.compared.lfc) <- c("Barcodes","P3", "B6", "Nude", "ThreeD")
write.table(p7.compared.lfc, file = "./output/normalized-data/P7LogFoldChangeData.txt", sep = "\t", row.names = FALSE)
p7.compared.lfc.gene.labeled <- cbind(gene.labeled.data[, "GeneID"], p7.compared.lfc)
colnames(p7.compared.lfc.gene.labeled) <- c("GeneID", "Barcodes","P3", "B6", "Nude", "ThreeD")
write.table(p7.compared.lfc.gene.labeled, file = "./output/normalized-data/P7LogFoldChangeData-GeneLabeled.txt", row.names = FALSE, sep = "\t")

# Nude
nude.compared.lfc <- normalized.lfc.data[, c("P7", "P3", "B6", "Nude", "ThreeD")] - 
  normalized.lfc.data[, "Nude"]

# Adding in geneID columns
nude.compared.lfc <- cbind(gene.labeled.data[, "Construct.IDs"], nude.compared.lfc)

# Removing P3 from data (it is all zeros) and renaming columns
nude.compared.lfc <- within(nude.compared.lfc, rm("Nude"))
colnames(nude.compared.lfc) <- c("Barcodes", "P7", "P3", "B6", "ThreeD")
write.table(nude.compared.lfc, file = "./output/normalized-data/NudeLogFoldChangeData.txt", sep = "\t", row.names = FALSE)
nude.compared.lfc.gene.labeled <- cbind(gene.labeled.data[, "GeneID"], nude.compared.lfc)
colnames(nude.compared.lfc.gene.labeled) <- c("GeneID", "Barcodes", "P7", "P3", "B6", "ThreeD")
write.table(nude.compared.lfc.gene.labeled, file = "./output/normalized-data/NudeLogFoldChangeData-GeneLabeled.txt", row.names = FALSE, sep = "\t")

# 3D
threeD.compared.lfc <- normalized.lfc.data[, c("P7", "P3", "B6", "Nude", "ThreeD")] - 
  normalized.lfc.data[, "ThreeD"]

# Adding in geneID columns
threeD.compared.lfc <- cbind(gene.labeled.data[, "Construct.IDs"], threeD.compared.lfc)

# Removing P3 from data (it is all zeros) and renaming columns
threeD.compared.lfc <- within(threeD.compared.lfc, rm("ThreeD"))
colnames(threeD.compared.lfc) <- c("Barcodes", "P7", "P3", "B6", "Nude")
write.table(threeD.compared.lfc, file = "./output/normalized-data/3DLogFoldChangeData.txt", sep = "\t", row.names = FALSE)
threeD.compared.lfc.gene.labeled <- cbind(gene.labeled.data[, "GeneID"], threeD.compared.lfc)
colnames(threeD.compared.lfc.gene.labeled) <- c("GeneID", "Barcodes", "P7", "P3", "B6", "Nude")
write.table(threeD.compared.lfc.gene.labeled, file = "./output/normalized-data/3DLogFoldChangeData-GeneLabeled.txt", row.names = FALSE, sep = "\t")

##############################################################################

# Creating data sets of expression relative to P3 and P7
p3.compared.lfc.sample <- normalized.sample.data.guide.level[, c("P7", "P3", "B6", "Nude", "ThreeD")] - 
  normalized.sample.data.guide.level[, "P3"]

# Adding in geneID columns
p3.compared.lfc.sample <- cbind(gene.labeled.data[, "Construct.IDs"], p3.compared.lfc.sample)

# Removing P3 from data (it is all zeros) and renaming columns
p3.compared.lfc.sample  <- within(p3.compared.lfc.sample , rm("P3"))
colnames(p3.compared.lfc.sample ) <- c("Barcodes","P7", "B6", "Nude", "ThreeD")

write.table(p3.compared.lfc.sample , file = "./output/normalized-data/P3LogFoldChangeDataSample.txt", sep = "\t", row.names = FALSE)

p7.compared.lfc.sample <- normalized.sample.data.guide.level[, c("P7", "P3", "B6", "Nude", "ThreeD")] - 
  normalized.sample.data.guide.level[, "P7"]

# Adding in geneID columns
p7.compared.lfc.sample<- cbind(gene.labeled.data[, "Construct.IDs"], p7.compared.lfc.sample)
# Removing P3 from data (it is all zeros) and renaming columns
p7.compared.lfc.sample<- within(p7.compared.lfc.sample, rm("P7"))
colnames(p7.compared.lfc.sample) <- c("Barcodes","P3", "B6", "Nude", "ThreeD")

write.table(p7.compared.lfc.sample, file = "./output/normalized-data/P7LogFoldChangeDataSample.txt", sep = "\t", row.names = FALSE)