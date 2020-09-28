source("./code/R/helpers/DataPlots.R")
library(ggfortify)
library(GGally)

raw.log.normalized <- read.table("./data/readcounts/merged/merged-lognorm-sgRNA-B6P3P7N3D.txt", header = TRUE)
collapsed.to.replicates <- read.table("./output/normalized-data/Collapsed-Replicates.txt", header = TRUE)
normalized.data <- read.table("./output/normalized-data/NormalizedData.txt", header = TRUE)
sample.level.data <- read.table("./output/normalized-data/Normalized-Replicates-Gene-Level.txt", header = TRUE)
normalized.data.guide.level <- read.table("./output/normalized-data/Replicate-Separate-Guide-Level-Normalized.txt", header = TRUE)

if (!dir.exists("./output/figures/qc/qc-data")){
    dir.create("./output/figures/qc/qc-data", recursive = TRUE)
}

# Creating sgRNA QC figures (distribution, ECDF)

suppressMessages(sgRNA_distribution(collapsed.to.replicates, x_steps = 0.3,
                   plot.file.path = "./output/figures/qc/sgRNADistrubtion.jpg"))

suppressMessages(sgRNA_cumulative_distribution(collapsed.to.replicates, 
                              plot.file.path = "./output/figures/qc/sgRNAecdf.jpg"))

suppressMessages(sgRNA_cumulative_distribution(collapsed.to.replicates[, c("P3", "P7", "B6", "FR_2DB6")], 
                              plot.file.path = "./output/figures/qc/sgRNAecdf_P3P7B6.jpg"))
                              
suppressMessages(sgRNA_distribution(collapsed.to.replicates[, c("P3", "P7", "B6", "FR_2DB6")], x_steps = 0.3,
                   plot.file.path = "./output/figures/qc/sgRNADistrubtion_P3P7B6.jpg"))             
suppressMessages(sgRNA_cumulative_distribution(collapsed.to.replicates[, c("Nude", "FR_Nude")], 
                              plot.file.path = "./output/figures/qc/sgRNAecdf_Nude.jpg"))
                              
suppressMessages(sgRNA_distribution(collapsed.to.replicates[, c("Nude", "FR_Nude")], x_steps = 0.3,
                   plot.file.path = "./output/figures/qc/sgRNADistrubtion_Nude.jpg")) 

suppressMessages(sgRNA_cumulative_distribution(collapsed.to.replicates[, c("ThreeD", "FR_3D")], 
                              plot.file.path = "./output/figures/qc/sgRNAecdf_3D.jpg"))
                              
suppressMessages(sgRNA_distribution(collapsed.to.replicates[, c("ThreeD", "FR_3D")], x_steps = 0.3,
                   plot.file.path = "./output/figures/qc/sgRNADistrubtion_3D.jpg"))                   

###############################################################################

# Creating PCA plots 

# The sequence of steps is:
#   1. Transpose the data so genes are columns and rows are conditions, label the columns with teh gen
#   2. Rename columns with GeneIDs
#   3. Remove any genes with 0 expression
#   4. Create a condition column with the row names (conditoins) as entries
#   5. Create pca plot and plot data

gene.level.normalized <- as.data.frame(t(normalized.data[, 2:length(colnames(normalized.data))]))
colnames(gene.level.normalized) <- normalized.data$GeneID
gene.level.normalized <- gene.level.normalized[, which(apply(gene.level.normalized, 2, var) != 0)]
gene.level.normalized$condition <- row.names(gene.level.normalized)
gene.level.normalized.plot <- prcomp(gene.level.normalized[, 1:(length(colnames(gene.level.normalized))-1)], scale = TRUE)

autoplot(gene.level.normalized.plot, data = gene.level.normalized, colour = "condition") + theme_bw()
suppressMessages(ggsave("./output/figures/qc/Normalized-Gene-Level-PCA.jpg"))

write.table(gene.level.normalized.plot$rotation[order(-gene.level.normalized.plot$rotation[, "PC1"]), ], "./output/figures/qc/qc-data/Normalized-Gene-Level-PCA.txt", row.names = TRUE, sep = "\t")
write.table(gene.level.normalized.plot$x, "./output/figures/qc/qc-data/Normalized-Gene-Level-PCA-Coordinates.txt", row.names = TRUE, sep = "\t")

# Normalization procedure applied prior to gene collapse
guide.level.normalized <- as.data.frame(t(normalized.data.guide.level[, 2:length(colnames(normalized.data.guide.level))]))
colnames(guide.level.normalized) <- normalized.data.guide.level$GeneID
guide.level.normalized <- guide.level.normalized[, which(apply(guide.level.normalized, 2, var) != 0)]
guide.level.normalized$condition <- row.names(guide.level.normalized)
guide.level.normalized.plot <- prcomp(guide.level.normalized[, 1:(length(colnames(guide.level.normalized))-1)], scale = TRUE)
autoplot(guide.level.normalized.plot, data = guide.level.normalized, colour = "condition") + theme_bw()
suppressMessages(ggsave("./output/figures/qc/Normalized-Guide-Level-PCA.jpg"))
write.table(guide.level.normalized.plot$rotation[order(-guide.level.normalized.plot$rotation[, "PC1"]), ], "./output/figures/qc/qc-data/Normalized-Guide-Level-PCA.txt", row.names = TRUE, sep = "\t")
write.table(guide.level.normalized.plot$x, "./output/figures/qc/qc-data/Normalized-Guide-Level-PCA-Coordinates.txt", row.names = TRUE, sep = "\t")

# No 3D condition in PCA plot
gene.level.normalized.no3d <- as.data.frame(t(normalized.data[, 2:(length(colnames(normalized.data))-1)]))
colnames(gene.level.normalized.no3d) <- normalized.data$GeneID
gene.level.normalized.no3d <- gene.level.normalized.no3d[, which(apply(gene.level.normalized.no3d, 2, var) != 0)]
gene.level.normalized.no3d$condition <- row.names(gene.level.normalized.no3d)
gene.level.normalized.no3d.plot <- prcomp(gene.level.normalized.no3d[, 1:(length(colnames(gene.level.normalized.no3d))- 1)], scale = TRUE)
autoplot(gene.level.normalized.no3d.plot, data = gene.level.normalized.no3d, colour = "condition") + theme_bw()
suppressMessages(ggsave("./output/figures/qc/Normalized-Gene-Level-No-3D-PCA.jpg"))

write.table(gene.level.normalized.no3d.plot$rotation[order(-gene.level.normalized.no3d.plot$rotation[, "PC1"]), ], "./output/figures/qc/qc-data/Normalized-Gene-Level-No-3D-PCA.txt", row.names = TRUE, sep = "\t")
write.table(gene.level.normalized.no3d.plot$x, "./output/figures/qc/qc-data/Normalized-Gene-Level-No-3D-PCA-Coordinates.txt", row.names = TRUE, sep = "\t")

# Normalization procedure after collapse to genes
sample.level.normalized <- as.data.frame(t(sample.level.data[, 2:length(colnames(sample.level.data))]))
sample.level.normalized <- sample.level.normalized[, which(apply(sample.level.normalized, 2, var) != 0)]
sample.level.normalized$Group <- c(rep("P7", 3), rep("P3", 3), rep("B6", 12), rep("Nude", 12), rep("3D", 3))
sample.level.normalized.plot <- prcomp(sample.level.normalized[, 1:(length(colnames(sample.level.normalized))- 1)], scale = TRUE)
autoplot(sample.level.normalized.plot, data = sample.level.normalized, colour = "Group") + theme_bw()
suppressMessages(ggsave("./output/figures/qc/Normalized-Sample-Gene-Level-PCA.jpg"))
write.table(sample.level.normalized.plot$rotation[order(-sample.level.normalized.plot$rotation[, "PC1"]), ], "./output/figures/qc/qc-data/Normalized-Sample-Gene-Level-PCA.txt", row.names = TRUE, sep = "\t")
write.table(sample.level.normalized.plot$x, "./output/figures/qc/qc-data/Normalized-Sample-Gene-Level-PCA-Coordinates.txt", row.names = TRUE, sep = "\t")

###############################################################################

###############################################################################
# Correlation Plots

suppressPackageStartupMessages(library(corrplot))
normalized.correlation.data <- normalized.data[, -1]
colnames(normalized.correlation.data) <- c("P7", "P3", "B6", "Nude", "3D")
normalized.correlations <- cor(normalized.correlation.data)
pdf("./output/figures/qc/Correlation-Normalized-Data.pdf")
corrplot.mixed(normalized.correlations, lower = "number", upper = "ellipse", tl.col = "black")
invisible(capture.output(dev.off()))

sample.correlations.data <- sample.level.data[, -1]
sample.correlations <- cor(sample.correlations.data)
pdf("./output/figures/qc/Correlation-Normalized-Sample-Data.pdf")
corrplot(sample.correlations, method = "color", order = "hclust", tl.col = "black")
invisible(capture.output(dev.off()))

# Guide Level Data
normalized.guide.correlation.data <- normalized.data.guide.level[, -1]
colnames(normalized.guide.correlation.data) <- c("P7", "P3", "B6", "Nude", "3D")
normalized.guide.correlations <- cor(normalized.guide.correlation.data)
pdf("./output/figures/qc/Correlation-Normalized-Data-Guide-Level.pdf")
corrplot.mixed(normalized.guide.correlations, lower = "number", upper = "ellipse", tl.col = "black")
invisible(capture.output(dev.off()))

# Pool Distribution Clustering (Condition Level)
pool.distribution.correlations.data <- dplyr::select_if(collapsed.to.replicates, is.numeric)
pool.correlations <- cor(pool.distribution.correlations.data)
pdf("./output/figures/qc/Correlation-Pool-Distributions-Condition-Level.pdf")
corrplot(pool.correlations, method = "color", order = "hclust", tl.col = "black")
invisible(capture.output(dev.off()))

# Pool Distribution Clustering (Sample Level)
sample.pool.distribution.correlations.data <- dplyr::select_if(raw.log.normalized, is.numeric)
sample.pool.correlations <- cor(sample.pool.distribution.correlations.data)
pdf("./output/figures/qc/Correlation-Pool-Distributions-Sample-Level.pdf")
corrplot(sample.pool.correlations, method = "color", order = "hclust", tl.col = "black")
invisible(capture.output(dev.off()))


# Correlation plot with added scatter plots

# Extracting numerical data and collapsing by GeneID
log_fold_change <- read.table("./output/normalized-data/NormalizedData.txt", header = TRUE)
log_fold_change.data <- dplyr::select_if(log_fold_change, is.numeric)
log_fold_change.collapsed <- aggregate(x = log_fold_change.data,
                                       by = list(GeneID = log_fold_change$GeneID),
                                       FUN = median)

correlation.data <- dplyr::select_if(log_fold_change.collapsed, is.numeric)
colnames(correlation.data) <- c("P7", "P3", "B6", "Nude", "3D")


color_correlaiton <- function(data, mapping, method="p", use="pairwise", ...){
  # Taken from stack overflow: 
  # https://stackoverflow.com/questions/45873483/ggpairs-plot-with-heatmap-of-correlation-values
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value  
  colFn <- colorRampPalette(c("#600f21", "#8ea6d5", "#0b1f3c"), interpolate = 'spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
  
  ggally_cor(data = data, mapping = mapping, color = "white", ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill)) 
}

pdf("./output/figures/qc/Correlation_Normalized_LogFoldChange_GeneLevel.pdf")
ggpairs(correlation.data, 
        upper = list(continuous = color_correlaiton), 
        progress = FALSE)
dev.off()
