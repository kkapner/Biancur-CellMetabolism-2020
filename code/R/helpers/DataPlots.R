# Library loading
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(UpSetR))

sgRNA_distribution <- function(normalized.data, x_steps, plot.file.path){
  #' This function produces overlaid histograms that display the breakdown
  #' of the normalized read counts. The plots are made for each sample
  #' provided (indicated by the column names in normalized.data). 
  #' This can be used t o see how much dropout occurs and how tight the 
  #' original distribution of sgRNA is.
  #' 
  #' @param normalized.data data frame of normalized read count data
  #' @param x_steps the increment for the x axis
  #' @param plot.file.path the full path to where the plot should be saved, 
  #' including the desired file name and type
  #' 
  #' @return Saves the created plot in the specified location with the 
  #' specified file name
  
  # Pulling out only numeric columns
  normalized.data.numeric <- dplyr::select_if(normalized.data, is.numeric)
  
  # Melting data for plotting
  plot.data <- melt(normalized.data.numeric, variable.name = "Sample", 
                    value.name = "Data")
  
  # Creating ggplot object
  sgRNAdist <- ggplot(plot.data, aes(x = Data, fill = Sample)) +
    geom_histogram(binwidth = 0.05, position="identity", alpha = 0.3) +
    xlab(expression(Log[2]~(average~normalized~read~count))) +
    ylab("Number of sgRNA guides") +
    ggtitle("Average sgRNA Distrubtion") +
    scale_x_continuous(breaks = seq(min(plot.data$Data),
                                    max(plot.data$Data), by = x_steps)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 75, hjust = 1),
          legend.position = "right")
  
  print(sgRNAdist)
  ggsave(plot.file.path)
}

sgRNA_cumulative_distribution <- function(normalized.data, plot.file.path){
  #' This function produces the empirical cumulative distribution function for 
  #' the normalized data. The plots are made for each sample
  #' provided (indicated by the column names in normalized.data). 
  #' This can be used t o see how much dropout occurs and how tight the 
  #' distribution of sgRNA is.
  #' 
  #' @param normalized.data data frame of normalized read count data
  #' @param plot.file.path the full path to where the plot should be saved, 
  #' including the desired file name and type
  #' 
  #' @return Saves the created plot in the specified location with the 
  #' specified file name
  
  if (nrow(normalized.data) == 0) {
    stop("Normalized.data is empty.")
  }
  
  # Pulling out only numeric columns
  normalized.data.numeric <- dplyr::select_if(normalized.data, is.numeric)
  
  # Melting data for plotting
  plot.data <- melt(normalized.data.numeric, variable.name = "Sample", 
                    value.name = "Data")
  
  # colnames(plot.data) <- "Data"
  plot.data$Data <- 2 ^ plot.data$Data
  
  
  # Creating ggplot object
  sgRNAdist <- ggplot(plot.data,aes(x = Data, color = Sample)) +
    stat_ecdf(geom = "step", position = "identity") +
    
    # Relabeling y axis and log scaling x
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    scale_x_log10() +
    
    # Labeling graph and setting theme
    ggtitle("Cumulative Fraction of Samples") +
    xlab("Read Count") +
    ylab("Cumulative Fraction") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(sgRNAdist)
  ggsave(plot.file.path)
}

gene_dependencies <- function(scored.data, condition, experiment.name,
                              plot.file.path){
  #' This function produces a plot of scores vs gene rank and labels the top
  #' and bottom points to illustrate gene dependencies. Can only be used
  #' when scored.data only has two columns (LogFoldChange and Score). 
  #' Cannot be used when keep.samples in fold_change or multisample in 
  #' z_score_normalize are set to TRUE.
  #' 
  #' @param scored.data data frame of z scored log fold change read count data
  #' @param condition name of condition from scored.data that you want to
  #' plot the gene dependencies for
  #' @param experiment.name condition being studied for naming the plot. 
  #' Title will read "expermient.name Dependencies"
  #' @param plot.file.path the full path to where the plot should be saved, 
  #' including the desired file name and type
  #' 
  #' @return Saves the created plot in the specified location with the 
  #' specified file name
  
  # Sorting in ascending order and labeling each with a Dependent Gene Rank,
  # which corresponds to its position from the lowest value in the list
  ordered.scored.data <- scored.data[order(pull(scored.data, condition)),]
  ordered.scored.data$Rank <- seq.int(nrow(ordered.scored.data))
  ordered.scored.data$Label <- 0
  
  # Creating top and bottom 5 genes
  ordered.scored.data["Label"][1:5, ] <- -1
  ordered.scored.data["Label"][(nrow(ordered.scored.data) - 5):(nrow(ordered.scored.data)), ] <- 1
  ordered.scored.data$Label[ordered.scored.data$GeneID == "Fdft1"] <- 1

  # Labeling control guides
  ordered.scored.data$ControlGuide <- is.na(ordered.scored.data$GeneID)
  
  # Creating plot
  # Creating ggplot object
  gene_rank <- ggplot(ordered.scored.data, aes_string(x = "Rank", 
                                                        y = condition)) +
    # Creating points
    geom_point(alpha = 0.5, color = "gray57") +
    geom_point(data = ordered.scored.data[ordered.scored.data$ControlGuide, ],
               alpha = 0.5, color = "red") +
    geom_text_repel(aes(label = ifelse(Label == 1, as.character(GeneID),'')),
              hjust = -1.5,vjust = 0, size = 3, nudge_y = 0.1) +
    geom_text_repel(aes(label = ifelse(Label == -1, as.character(GeneID),'')),
                    hjust = -1.5,vjust = 0, size = 3, nudge_x = 0.1, 
                    direction = "y") +
    
    # Labeling graph and axes
    ggtitle(paste0(experiment.name, " Depdendencies")) +
    xlab("Dependent Gene Rank") +
    ylab(expression("Score")) +
    
    # Setting theme
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(gene_rank)
  ggsave(plot.file.path)
}

volcano_plot <- function(raw.plot.data, x_data, y_data, x_label, y_label,
                         genes_of_interest, plot.title, plot.file.path) {
  #' Given a data frame with fold changes as one column and p-values as the
  #' second column and GeneIDs as the third column, this function plots the 
  #' fold change on the x-axis and the  p-values on the y-axis to create a 
  #' volcano plot. The function will automatically rename the columns to xlabel 
  #' and ylabel, so no standardized naming for the columns is needed. The top 
  #' 30 genes above 1.031 are highlighed  (as well as any genes in 
  #' genes_of_interest)
  #' 
  #' @param raw.plot.data data frame with two columns. The first column is the 
  #' measured change and the second column are the p values for plotting. NOTE:
  #' This function does not perform any p-value adjustments for multiple 
  #' comaprisons, it assumes this is already applied.
  #' @param x_data the column name correpsonding to the x axis data
  #' @param y_data the column name corresponding to the y axis data
  #' @param x_label the label for the x-axis of the graph (character or 
  #' expression)
  #' @param y_label the label for the y-axis of the graph (character or 
  #' expression)
  #' @param genes_of_interest character array of genes whose labels should be 
  #' plotted regardless of the significance
  #' @param plot.title the title for the graph (character or expression)
  #' @param plot.file.path the full file path, including the desired file name 
  #' and type to be used for saving the created plot.
  
  # Extracting columns for plotting
  plot.data <- raw.plot.data[, c(x_data, y_data, "GeneID")]
  
  # Creating column for significance
  plot.data <- plot.data[order(plot.data[, y_data], decreasing = TRUE), ]
  plot.data[, "Significant"] <- plot.data[, y_data]  > 1.031
  plot.data[, "LabelName"] <- plot.data[, "Significant"] & 
    (plot.data[, "GeneID"] %in% plot.data[1:30, "GeneID"])
  plot.data[, "LabelName"] <- plot.data[, "LabelName"] | 
    (plot.data[, "GeneID"] %in% genes_of_interest)
  
  volcano <- ggplot(data = plot.data, aes_string(x = x_data, y = y_data)) +
    geom_point() +
    geom_text_repel(aes(label = ifelse(plot.data[, "LabelName"], 
                                 as.character(GeneID),
                                 "")), 
              hjust = 0, vjust = 0, size = 3) +
    xlab(x_label) +
    ylab(y_label) +
    theme_bw() +
    ggtitle(plot.title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(volcano)
  ggsave(plot.file.path)
}

intersection_plots <- function(scored.data, score.condition, threshold, below,
                               category.column, save.path){
  #' This function uses the UpSetR package to create intersection plots
  #' for the given categories and given score conditions. It thresholds the 
  #' score.condition at the provided threshold value and then creates 
  #' intersections of the given category.column for the data
  #' 
  #' @param scored.data list of data frames containing the data to be used for 
  #' finding the intersection and for thresholding
  #' @param score.condition the condition to threshold
  #' @param threshold the value to threshold (either above or below)
  #' @param below logcial value indicating if values below threshold should be
  #' taken
  #' @param category.column the column for the data to be used for measuring 
  #' intersections
  #' @param save.path path to where image file should be saved, including
  #' file exntension. For jpg files, use the extension jpeg NOT jpg.
  #' 
  #' @return None
  
  if (below) {
    filtered.data <- lapply(scored.data, function(x){
      dplyr::filter(x, get(score.condition) < threshold)
    })} else {
      filtered.data <- lapply(scored.data, function(x){
        dplyr::filter(x, get(score.condition) > threshold)
      })}
  
  filtered.categories <- lapply(filtered.data, function(x){
    as.character(x[, category.column])
  })
  
  upset_plot <- upset(fromList(filtered.categories), order.by = "freq")
  print(upset_plot)
  dev.copy(get(tools::file_ext(save.path)), save.path)
  dev.off()
}

threshold_volcano <- function(scored.data, score.condition, threshold, below,
                              x_axis, y_axis, xlabel, ylabel, point.labels,
                              interest.genes, plot.title, save.path){
  #' This function returns a volcano plot of the scored data, only coloring and
  #' labeling the points that meet the thresholding requirements for the 
  #' score.condition. This does not have to be the same data as the provided
  #' y_axis. Additionally, genes of interest can be supplied and these will be 
  #' labeled (but only colored if significant).
  #' 
  #' @param scored.data the data frame containing all of the columns that are to
  #' be used in plotting (it can have additional columns)
  #' @param score.condition the condition to apply thresholding to
  #' @param threshold a numerical value for comparison in the score.condition
  #' @param below should data below or above the threshold be selected
  #' @param x_axis the data to be used as the x axis of the plot
  #' @param y_axis the data to be used as the y axis of the plot
  #' @param xlabel the label for the x axis
  #' @param ylabel the label for the y axis
  #' @param point.labels the column name containing the labels for the points
  #' @param interest.genes character array of additional genes of interest to
  #' label on the plot
  #' @param plot.title the title for the plot
  #' @param save.path the full file path (including file name and extension) of 
  #' where to save the plot
  
  scored.data[, c("Color", "Label")] <- FALSE

  if (below) {
    scored.data <- scored.data[order(scored.data[, score.condition]), ]
    scored.data$Color[scored.data[, score.condition] < threshold] <- TRUE
    scored.data$Label[scored.data[, score.condition] < threshold] <- TRUE
    scored.data$Label[26:nrow(scored.data)] <- FALSE
  } else {
    scored.data <- scored.data[order(-scored.data[, score.condition]), ]
    scored.data$Color[scored.data[, score.condition] > threshold] <- TRUE
    scored.data$Label[scored.data[, score.condition] > threshold] <- TRUE
    scored.data$Label[26:nrow(scored.data)] <- FALSE
  }
  
  
  scored.data$Label[scored.data[, point.labels] %in% interest.genes] <- TRUE
  
  volcano <- ggplot(scored.data, aes_string(x = x_axis, y = y_axis)) +
    geom_point(color = ifelse(scored.data[, "Color"], "blue", "gray50")) +
    geom_text_repel(aes(label = ifelse(scored.data[, "Label"], 
                                     as.character(get(point.labels)),
                                     "")), 
                  hjust = 0, vjust = 0, size = 3) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(plot.title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
    
  print(volcano)
  ggsave(save.path)
}

comparison_scatter <- function(plot.data, x_axis, y_axis, plot.title, 
                               plot.file.path){
  #' Produces a generic plot where the data corresponding to the provided 
  #' paramter (column name) y_axis gets plotted against the provided parameter
  #' x_axis. The linear regression line between the two variables is determined
  #' and plotted.
  #' 
  #' @param plot.data the data frame containing the two variables of interest
  #' @param x_axis the column name corresponding to the variable to be placed
  #' on the x_axis
  #' @param y_axis the column name corresponding to the variable to be placed
  #' on the y_axis
  #' @param plot.title the title for the plot
  #' @param plot.file.path the full path to where the plot should be saved, 
  #' including the desired file name and type
  #' 
  #' @return Saves the created plot in the specified location with the 
  #' specified file name
  
  # Determining plot bounds
  total_min <- min(min(plot.data[, y_axis]), min(plot.data[, x_axis])) - 0.1
  total_max <- max(max(plot.data[, y_axis]), max(plot.data[, x_axis])) + 0.1
  
  ggplotRegression <- function(fit) {
    
    # Getting r^2 value
    r_squared = signif(summary(fit)$adj.r.squared, 3)
    
    # Plotting
    ggplot(fit$model, aes_string(x = names(fit$model)[2], 
                                 y = names(fit$model)[1])) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue") +
      annotate(geom = "text",
               x = -1, y = 0.5,
               label = paste("R^2 == ", r_squared), parse = TRUE) +
      # Setting bounds for graph so xlim = ylim
      xlim(total_min, total_max) +
      ylim(total_min, total_max) +
      xlab(x_axis) +
      ylab(y_axis) +
      ggtitle(plot.title) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  # Creating plot object
  formla <- as.formula(paste(y_axis, "~", x_axis))
  print(ggplotRegression(lm(formla, data = plot.data)))
  ggsave(plot.file.path)
}

control_validation_plots <- function(plot.data, condition, experiment.name, 
                                     positive.controls, plot.file.path){
  #' Given a data frame of gene and sample collapsed log fold changes which are
  #' un-normalized, this function will create plots for the gene dependencies
  #' for each sample, highlighting positive controls in green and negative
  #' controls in red. Positive controls are expected to drop out (be at the 
  #' lower end of the gene ranking), while negative controls should be centered
  #' around 0. This function also places a symbol on the plot representing
  #' the median value for the controls that will be used in the normalization
  #' procedure.
  #' 
  #' @param plot.data data frame of z scored log fold change read count data
  #' @param condition name of condition from scored.data that you want to
  #' plot the gene dependencies for
  #' @param experiment.name condition being studied for naming the plot. 
  #' Title will read "expermient.name Dependencies"
  #' @param positive.controls an array of GeneIDs to use as positive controls
  #' @param plot.file.path the full path to where the plot should be saved, 
  #' including the desired file name and type
  #' 
  #' @return Saves the created plot in the specified location with the 
  #' specified file name
  
  # Sorting in ascending order and labeling each with a Dependent Gene Rank,
  # which corresponds to its position from the lowest value in the list
  ordered.plot.data <- plot.data[order(pull(plot.data, condition)),]
  ordered.plot.data$Rank <- seq.int(nrow(ordered.plot.data))
  ordered.plot.data$Label <- "gray75"
  
  # Creating coloring for positive (green) and negative (red) control guides
  ordered.plot.data$Label[ordered.plot.data$GeneID %in% positive.controls] <- "#4dac26"
  ordered.plot.data$Label[is.na(ordered.plot.data$GeneID)] <- "#d01c8b"
  
  # Finding median of positive controls
  positive.control.data <- pull(ordered.plot.data, condition)
  positive.median <- median(positive.control.data[ordered.plot.data$GeneID 
                                                  %in% positive.controls])
  
  # Getting negative control value
  negative.control.data <- pull(ordered.plot.data, condition)
  negative.control <- negative.control.data[is.na(ordered.plot.data$GeneID)]
  
  # Creating ggplot object
  gene_rank <- ggplot(ordered.plot.data, aes_string(x = "Rank", 
                                                    y = condition, 
                                                    colour = "Label")) +
    # Creating points
    geom_point(alpha = 0.5) +
    
    scale_color_manual(values = c("#4dac26","#d01c8b","gray75"),
                       labels = c("Positive Control", "Negative Control",
                                  "Non-Control")) +
    
    geom_point(data = subset(ordered.plot.data, Label == "#4dac26"),
               aes_string(x = "Rank", y = condition), colour = "#4dac26") +
    geom_point(data = subset(ordered.plot.data, Label == "#d01c8b"),
               aes_string(x = "Rank", y = condition), colour = "#d01c8b") +
    geom_hline(yintercept = positive.median, linetype = "dashed",
               color = "#4dac26") +
    geom_hline(yintercept = negative.control, linetype = "dashed",
               color = "#d01c8b") +
  
    # Labeling graph and axes
    ggtitle(paste0(experiment.name, " Control Validation")) +
    xlab("Dependent Gene Rank") +
    ylab(expression("LogFoldChange")) +
    
    # Setting theme
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right", legend.title = element_blank())

  print(gene_rank)
  ggsave(plot.file.path)
}


