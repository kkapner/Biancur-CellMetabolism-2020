suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(plyr))


#### Basic Read Functions ####
read_tab_scores <- function(reads.file.path, remove.blank = TRUE){
  #' Reads in RNA count data from a tab-delimited file with RNA constructs as 
  #' rows and samples as columns. This function removes the blank column by 
  #' default, but can be set to retain the value. 
  #' 
  #' @param reads.file.path path to tab-delimited txt file containing read 
  #' count data
  #' @param remove.blank removes blank sample from data. default value is TRUE.
  #' 
  #' @return a tibble of the data for further processing, normalization, and
  #' grouping
  
  # Loading in data as tibble
  read.data <- dplyr::as_tibble(
    read.delim(reads.file.path, header = TRUE, sep = "\t", dec = ".", 
               stringsAsFactors = FALSE))
  
  # Removing blanks if desired
  if (remove.blank){
    read.data.no.blanks <- suppressWarnings(select(read.data, 
                                                   -one_of(c("blank", "BLANK",
                                                       "Blank"))))
    return (read.data.no.blanks)
  } else {
    return (read.data)
  }
}

barcode_to_genes <- function(read.counts, ids.to.genes.path){
  #'  Converts construct IDs to gene names in the read.counts tibble using an
  #'  xlsx file that matches the Clone ID of a given target sequence to the 
  #'  targeted gene name.
  #'  
  #'  @param read.counts a tibble of read counts
  #'  @param ids.to.genes.path full file path to the xlsx file containing the
  #'  mapping from constructs to gene names.
  #'  
  #'  @return a modified read.counts tibble with an appended column representing
  #'  the targeted gene.
  
  # Reading in map data from xlsx file
  mapping.data <- read_xlsx(ids.to.genes.path)
  
  # Extracting Clone ID (Construct.IDs) and gene names
  clone.ids <- pull(mapping.data["Clone ID"])
  gene.names <- pull(mapping.data["Orig. Target Gene Symbol"])
  
  # Mapping values and returning updated data
  mapped.construct.ids <- plyr::mapvalues(read.counts$Construct.IDs, 
                                          clone.ids,
                                          gene.names)
  gene.mapped.read.counts <- read.counts
  gene.mapped.read.counts["GeneID"] <- mapped.construct.ids
  
  return (gene.mapped.read.counts)
}

### Multi-plate Read Functions ###

multi_plate_read <- function(multiple.reads.file.paths, remove.blank = TRUE){
  #' Given a character array of file paths to normalized read data from the
  #' same experiment (all guides are the same), this function reads in all the
  #'  files and merges them together (ensures data stays matched to the barcodes)
  #' If column names match between given plates, they are all renamed with a
  #' plate identifier put at the end of the name.
  #' 
  #' @param multiple.reads.file.paths character array of full file paths to
  #'  tab-delimited txt files containing read count data
  #' @param remove.blank removes blank sample from data. default value is TRUE.
  #' 
  #' @return a tibble of the data for further processing, normalization, and
  #' grouping
  
  # Creating list of data frames
  merge.data.frames <- list()
  for (reads_i in 1:length(multiple.reads.file.paths)) {
    merge.data.frames[[reads_i]] <- read_tab_scores(
      multiple.reads.file.paths[reads_i], remove.blank)
  }
  
  # Merging data frames
  merged.data <- Reduce(function(...) merge(..., 
                                            by = c("Construct.Barcode", 
                                                   "Construct.IDs"),
                                            all = TRUE), merge.data.frames)
  
  return(merged.data)
  
}

