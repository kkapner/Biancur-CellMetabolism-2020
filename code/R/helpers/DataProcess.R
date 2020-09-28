suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(matrixStats))

z_score <- function(normalized.data){
  #' This function returns the z-scored normalized data. It computes the z-score
  #' for each sample (column) and returns a data frame of the same dimensions
  #' as normalized.data
  #' 
  #' @param normalized.data a data frame with samples as columns and log2
  #' fold change data as values
  #' 
  #' @return a data frame of same shape as normalized.data with the z-score
  #' for each log2 fold change value
  
  # Initializing data frame
  z.scored.data <- data.frame(matrix(NA,
                                     nrow = nrow(normalized.data),
                                     ncol = 1))[-1]
  
  # Extracting non numeric columns
  z.scored.data[, c("Construct.Barcode", "Construct.IDs", "GeneID")] <-
    dplyr::select_if(normalized.data, negate(is.numeric))
  
  # Calculating z-score
  numeric.data <- dplyr::select_if(normalized.data, is.numeric)
  numeric.data <- scale(numeric.data)  
  
  # Appending to z.scored.data
  z.scored.data <- cbind(z.scored.data, numeric.data)
  return(z.scored.data)
}

combine_by_id <- function(normalized.data, id){
  #' This function collapses the data in normalized.data to a per gene level 
  #' (if GeneID is supplied for value) by taking the median of each of the values 
  #' in the columns for a given gene.
  #' 
  #' @param normalized.data a data frame to be collapsed by the indicated id
  #' @param id the id to group and the collapse the data by
  #' 
  #' @return the collapsed data frame 
  
  collapsed.data <- normalized.data %>% dplyr::group_by_at(id) %>% 
    dplyr::summarise_all("median")
  
  if (id == "GeneID"){
    collapsed.data$Construct.Barcode <- NULL
    collapsed.data$Construct.IDs <- NULL
    collapsed.data$`"GeneID"` <- NULL
  }
  
  return(collapsed.data)
}

analyze_by_gene <- function(normalized.data, group1, group2, method, 
                            adjustment.method = "BH"){
  #' This function performs the statistical test indicated by method (t, mann)
  #' for each gene between group1 and group2. The p-values are then adjusted
  #' using the indicated method.
  #' 
  #' @param normalized.data a data frame of z-score normalized log fold change 
  #' data that is not collapsed by gene (multiple guides point to the same gene)
  #' @param group1 the name of the first group to perform comparison. This name 
  #' must match one of the columns of normalized. 
  #' @param group2 the name of the second group to perform comparison. This name 
  #' must match one of the columns of normalized. 
  #' @param method the statistical method to perform to compare the two groups.
  #' Current options are t (for the two-tailed t-test with Welch's correction
  #' for unequal variances) or mann for the Mann-Whitney U test.
  #' @param adjustment.method the p-value multiple comparisons adjustment to 
  #' be used after all comparisons are performed. The default is the 
  #' Benjamini-Hochberg correction, but any of the methods available to
  #' p.adjust are available.
  #' 
  #' @return a new data frame with the columns GeneID, MeanDifference, 
  #' p_value, p_value_adj, neg_log10_p_value, neg_log10_p_value_adj)
  
  # Getting unique genes for iterating 
  unique.genes <- unique(normalized.data$GeneID)
  
  # Initializing empty data frame 
  analyzed.data <- data.frame(matrix(NA,
                                     nrow = length(unique.genes),
                                     ncol = 1))[-1]
  
  # Iterating over each gene and performing test
  for (gene in unique.genes){
    gene.subset <- dplyr::filter(normalized.data, GeneID == gene)
    group1.data <- dplyr::select(gene.subset, group1)
    group2.data <- dplyr::select(gene.subset, group2)
    
    # Performign the test
    if (method == 't'){
      test.results <- t.test(group1.data, group2.data, 
                             alternative = c("two.sided"),
                             mu = 0, var.equal = FALSE,
                             conf.level = 0.95)
    } else if (method == "mann"){
      test.results <- wilcox.test(group1.data[[1]], group2.data[[1]],
                                  alternative = c("two.sided"),
                                  mu = 0, paired = FALSE,
                                  conf.level = 0.95)
    }

    mean_difference <- mean(group1.data[[1]]) - mean(group2.data[[1]])
    analyzed.data <- rbind(analyzed.data, 
                           data.frame("GeneID" = gene, 
                                      "MeanDifference" = mean_difference, 
                                      "p_value" = test.results$p.value))
  }
  
  # Performing p-value adjustment 
  analyzed.data$p_value_adj <- p.adjust(analyzed.data$p_value, 
                                        method = adjustment.method)
  
  # Taking negative logs of p_values
  analyzed.data$neg_log10_p_value <- -log10(analyzed.data$p_value)
  analyzed.data$neg_log10_p_value_adj <- -log10(analyzed.data$p_value_adj)
  
  return(analyzed.data)
  
}

normalize_data <- function(data, positive.controls,
                           pos_normalization, collapsed) {
  #' Given data that is median collapsed to gene level and a array of 
  #' positive controls, this function normalizes the data in each column
  #' using the following procedure:
  #' 1) Subtract the median of the negative controls from each value
  #' 2) Divide by the magnitude of the median of the positive controls
  #' This function returns a data frame of the same dimensions as the 
  #' median.collapsed.to.genes data frame
  #' 
  #' @param data data frame of gene level data
  #' @param positive.controls array of genes to use as positive controls
  #' @param pos_normalization a boolean value indicating if division by
  #' the absolute value of the median magniutde of the positive controls 
  #' should occur. If the data should just be centered on the
  #' negative controls, then this value should be FALSE
  #' @param collapsed a boolean value indicating if the data is already 
  #' collapsed to gene or if it is at the guide level
  #' 
  #'  @return a data frame of the same dimension as data with the normalzied
  #'  data that can be used for further analysis
  
  # Getting indices of negative and positive control medians
  negative.median.ind <- which(grepl("Control", data$GeneID))
  positive.median.ind <- which(data$GeneID %in%
                                 positive.controls)

  if (collapsed){
  # Subtracting negative control value (data is already collapsed) 
  negative.subtracted <- mutate_if(data, is.numeric,
                          function(x){x - x[negative.median.ind]})
  } else {
    # Subtracting median of negative controls (data is not collapsed)
    negative.subtracted <- mutate_if(data, is.numeric,
                                     function(x){x - 
                                         median(x[negative.median.ind])})
  }
  
  if (pos_normalization){
  # Dividing by the magnitude of the median of the positive controls
  normalized <- mutate_if(negative.subtracted, is.numeric,
                          function(x){x / abs(median(x[positive.median.ind]))})
  return(normalized)
  } else {
    return(negative.subtracted)
  }
}

match_stars <- function(stars.input, stars.results, condition){
  #' Given data frames consisting of the input to STARS and the resultant
  #' output from STARS, this function combines the results together into a 
  #' data frame, which can the be used to create volcano plots and be used 
  #' for further analysis. It adds in a -log base 10 column of the q-values
  #' and makes sure that the raw data is matched up to the correct genes.
  #' 
  #' @param stars.input the input data used to genereate STARS
  #' @param stars.results the results data from STARS
  #' @param condition the condition the STARS results are for
  #' 
  #' @return a data frame representing the combined data between the stars input
  #' and stars.results
  
  # First have to find differences between groups (in the event of thresholding
  # for STARS)
  
  set.diffs <- union(setdiff(stars.input[, "GeneID"],
                             stars.results[, "Gene.Symbol"]),
                     setdiff(stars.results[, "Gene.Symbol"],
                             stars.input[, "GeneID"]))
  
  # Filtering to make sure they have the same genes
  condition.filtered <- dplyr::filter(stars.results,
                                      !(Gene.Symbol %in% set.diffs))
  input.filtered <- dplyr::filter(stars.input,
                                  !(GeneID %in% set.diffs))
  
  # Median collapsing down to gene level
  input.filtered <- combine_by_id(input.filtered, "GeneID")
  
  # Ordering condition.filtered by gene name (input.filtered is already sorted)
  condition.filtered <- 
    condition.filtered[order(condition.filtered[, "Gene.Symbol"]), ]
  
  # Combining data and taking -log10 of q-values
  combined.data <- cbind(condition.filtered, input.filtered[, c("GeneID", 
                                                                condition)])
  
  combined.data[, "LogQ"] <- -log10(combined.data[, "q.value"])
  
  # If the q-value is 0, then this means that the value is < 3.4258307635e-07,
  # which corresponds to a LogQ of 6.465234
  combined.data[, "LogQ"][combined.data[, "LogQ"] == Inf] <- 6.465234
  
  return(combined.data)
  
}
