library(matrixStats)

group_by_sample <- function(normalized.read.counts, group.names){
  #' Produces a median reading for the given group names. It looks
  #' for just the main group name provided and averages all of the replicates.
  #' For example, FR.1, FR.2, and FR.3 will all be averaged if FR is provided
  #' in group.names
  #' 
  #' @param normalized.read.counts a normalized read count matrix with sample
  #' names for column names
  #' @param group.names the general group name for each sample (such as 
  #' Passage7 or FR or RC1), provided as an array
  #' 
  #' @return data frame with a single median value for each group. 
  
  if (length(group.names) == 0){
    stop("Group names not provided.")
  }
  if (nrow(normalized.read.counts) == 0){
    stop("Normalized read counts has no data.")
  }
  
  # Creating data frame for grouped reads and labeling by gene ID
  grouped.read.counts <- data.frame(matrix(NA, 
                                           nrow = nrow(normalized.read.counts),
                                           ncol = 1))[-1]
  
  grouped.read.counts["GeneID"] <- normalized.read.counts["GeneID"]
  
  # Iterate over each group name
  for (group in 1:length(group.names)){
    name.group <- group.names[group]
    
    # Extracting columns where grepl finds a pattern match
    group.data <- normalized.read.counts[, grepl(name.group, 
                                                 colnames(
                                                   normalized.read.counts))]
    # Taking row medians (if more than one column present)
    if (!is.null(dim(group.data))){
      group.data.medians <- rowMedians(as.matrix(group.data))
    } else {
      group.data.medians <- group.data
    }
    # Appending to grouped.read.counts
    grouped.read.counts[name.group] <- group.data.medians
  }
  return(grouped.read.counts)
}

