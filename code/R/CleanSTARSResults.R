#!/usr/bin/env Rscript

#Threshold 60
#Merging Individual Data
pos.ind.stars.files <- list.files("./output/stars/individual/pos", full.names = TRUE)
names(pos.ind.stars.files) <- lapply(pos.ind.stars.files, function(x){sub("\\_.*", "", basename(x))})
neg.ind.stars.files <- list.files("./output/stars/individual/neg", full.names = TRUE)
names(neg.ind.stars.files) <- lapply(neg.ind.stars.files, function(x){sub("\\_.*", "", basename(x))})

 pos.ind.stars.data <- lapply(pos.ind.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})
 neg.ind.stars.data <- lapply(neg.ind.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})

 merged_frames_complete <- list()

 for (cond in names(pos.ind.stars.data)){
   # Merging data frames just by rows
   pos_cond <- pos.ind.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
   neg_cond <- neg.ind.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
   neg_cond[, "STARS.Score"] <- -1*neg_cond[, "STARS.Score"]
   merged <- rbind(pos_cond, neg_cond)

   merged_frames_complete[[cond]] <- merged
   write.table(merged_frames_complete[[cond]], paste0("./output/stars/individual/STARS-Complete-", cond, ".csv"), row.names = FALSE, sep = "\t")
}

 # Merging P3 Comparison Data
 pos.comp3.stars.files <- list.files("./output/stars/comparison/p3/pos", full.names = TRUE)
 names(pos.comp3.stars.files) <- lapply(pos.comp3.stars.files, function(x){sub("\\_.*", "", basename(x))})
 neg.comp3.stars.files <- list.files("./output/stars/comparison/p3/neg", full.names = TRUE)
 names(neg.comp3.stars.files) <- lapply(neg.comp3.stars.files, function(x){sub("\\_.*", "", basename(x))})

 pos.comp3.stars.data <- lapply(pos.comp3.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})
 neg.comp3.stars.data <- lapply(neg.comp3.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})

 merged_frames_p3_full <- list()
 for (cond in names(pos.comp3.stars.data)){
   # Merging data frames just by rows
   pos_cond <- pos.comp3.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
   neg_cond <- neg.comp3.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
   neg_cond[, "STARS.Score"] <- -1*neg_cond[, "STARS.Score"]
   merged <- rbind(pos_cond, neg_cond)
   merged_frames_p3_full[[cond]] <- merged

   write.table(merged_frames_p3_full[[cond]], paste0("./output/stars/comparison/p3/STARS-Complete-", cond, ".csv"), row.names = FALSE, sep = "\t")
 }

 # Merging P7 Comparison Data
 pos.comp7.stars.files <- list.files("./output/stars/comparison/p7/pos", full.names = TRUE)
 names(pos.comp7.stars.files) <- lapply(pos.comp7.stars.files, function(x){sub("\\_.*", "", basename(x))})
 neg.comp7.stars.files <- list.files("./output/stars/comparison/p7/neg", full.names = TRUE)
 names(neg.comp7.stars.files) <- lapply(neg.comp7.stars.files, function(x){sub("\\_.*", "", basename(x))})

 pos.comp7.stars.data <- lapply(pos.comp7.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})
 neg.comp7.stars.data <- lapply(neg.comp7.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})

 merged_frames_p7_full <- list()

 for (cond in names(pos.comp7.stars.data)){
   # Merging data frames just by rows
   pos_cond <- pos.comp7.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
   neg_cond <- neg.comp7.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
   neg_cond[, "STARS.Score"] <- -1*neg_cond[, "STARS.Score"]
   merged <- rbind(pos_cond, neg_cond)
   merged_frames_p7_full[[cond]] <- merged

   write.table(merged_frames_p7_full[[cond]], paste0("./output/stars/comparison/p7/STARS-Complete-", cond, ".csv"), row.names = FALSE, sep = "\t")
 }

 # Merging Nude Comparison Data
 pos.nude.stars.files <- list.files("./output/stars/comparison/nude/pos", full.names = TRUE)
 names(pos.nude.stars.files) <- lapply(pos.nude.stars.files, function(x){sub("\\_.*", "", basename(x))})
 neg.nude.stars.files <- list.files("./output/stars/comparison/nude/neg", full.names = TRUE)
 names(neg.nude.stars.files) <- lapply(neg.nude.stars.files, function(x){sub("\\_.*", "", basename(x))})

 pos.nude.stars.data <- lapply(pos.nude.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})
 neg.nude.stars.data <- lapply(neg.nude.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})

 merged_frames_nude_full <- list()

 for (cond in names(pos.nude.stars.data)){
  # Merging data frames just by rows
   pos_cond <- pos.nude.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
   neg_cond <- neg.nude.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
   neg_cond[, "STARS.Score"] <- -1*neg_cond[, "STARS.Score"]
   merged <- rbind(pos_cond, neg_cond)
   merged_frames_nude_full[[cond]] <- merged

   write.table(merged_frames_nude_full[[cond]], paste0("./output/stars/comparison/nude/STARS-Complete-", cond, ".csv"), row.names = FALSE, sep = "\t")
}

pos.3d.stars.files <- list.files("./output/stars/comparison/3D/pos", full.names = TRUE)
names(pos.3d.stars.files) <- lapply(pos.3d.stars.files, function(x){sub("\\_.*", "", basename(x))})
neg.3d.stars.files <- list.files("./output/stars/comparison/3D/neg", full.names = TRUE)
names(neg.3d.stars.files) <- lapply(neg.3d.stars.files, function(x){sub("\\_.*", "", basename(x))})

pos.3d.stars.data <- lapply(pos.3d.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})
neg.3d.stars.data <- lapply(neg.3d.stars.files, function(x){read.delim(x, header = TRUE, skip = 1)})

merged_frames_3d_full <- list()
for (cond in names(pos.3d.stars.data)){
  # Merging data frames just by rows
  pos_cond <- pos.3d.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
  neg_cond <- neg.3d.stars.data[[cond]][c("Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value")]
  neg_cond[, "STARS.Score"] <- -1*neg_cond[, "STARS.Score"]
  merged <- rbind(pos_cond, neg_cond)
  merged_frames_3d_full[[cond]] <- merged
  
  write.table(merged_frames_3d_full[[cond]], paste0("./output/stars/comparison/3D/STARS-Complete-", cond, ".csv"), row.names = FALSE, sep = "\t")
}

