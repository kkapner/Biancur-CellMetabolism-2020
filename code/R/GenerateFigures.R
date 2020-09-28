library(ggplot2)
library(ggrepel)

create.volcano <- function(plot.data, title, save.path){
  column.ids <- colnames(plot.data)
  colnames(plot.data) <- c("Gene", "LFC", "logp")
  
  volcano <- ggplot(plot.data, aes(x = LFC, y = logp)) +
    geom_point() +
    xlab(column.ids[2]) +
    ylab(column.ids[3]) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  volcano + geom_text_repel(data = head(plot.data, 50), aes(label = Gene), parse=TRUE)
  
  ggsave(save.path)
}

create.delta.delta <- function(plot.data, title, save.path, x.axis.group, comparison.group){
  
  positive.controls <- c("Rpl10", "Rpl13a", "Rpl14", "Rpl17", "Rpl23",
                        "Rpl27", "Rpl35a", "Rpl7", "Rplp2", "Rps11",
                       "Polr2d", "U2af2", "Eef2", "Psmd1", "Ran",
                       "Tubb5", "Psmb2", "Rptor", "Polr2i", "Psmd6",
                       "Supt5", "Psmc4")
  positive.controls <- paste0("italic('", positive.controls, "')")
  column.ids <- colnames(plot.data)
  colnames(plot.data) <- c("Gene", "DD", "logp", "Comparison")
  pos.control.median.data <- dplyr::filter(plot.data, Gene %in% positive.controls)
  pos.control.median.x <- -1
  pos.control.median.y <- median(pos.control.median.data[, "DD"])
  pos.control.mad.y <- mad(pos.control.median.data[, "DD"])
  pos.control.mad.x <- mad(pos.control.median.data[, "Comparison"])
  negative.control.data <- dplyr::filter(plot.data, Gene == "italic('Control-target')")
  negative.control.x <- negative.control.data[, "Comparison"]
  negative.control.y <- negative.control.data[, "DD"]

  
  delta_delta <- ggplot(plot.data, aes(x = Comparison, y = DD)) +
    geom_point(data = dplyr::filter(plot.data, !(Gene %in% positive.controls)), color = "#d1d1d1") +
    xlab(bquote(.(parse(text=column.ids[4])))) +
    ylab(bquote(.(parse(text=column.ids[2])))) +
    ggtitle(title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept=0, colour = "#b8b8b8", linetype = "longdash") +
    geom_vline(xintercept=0, colour = "#b8b8b8", linetype = "longdash")

  
  delta_delta + geom_point(data = dplyr::filter(plot.data, logp > -log10(0.05) & Comparison < 0 & DD < 0), 
               aes(x = Comparison, y = DD), 
               colour = '#2c7bb6') +
    geom_point(data = dplyr::filter(plot.data, logp > -log10(0.05) & 
    Comparison < 0 & DD > 0), 
            aes(x = Comparison, y = DD), 
            colour = '#75573f') +
    geom_point(data = dplyr::filter(plot.data, logp > -log10(0.05) & Comparison > 0), aes(x = Comparison, y = DD), colour = "#fdae61") +
    geom_point(data = dplyr::filter(plot.data, Gene == "Control-target"), 
               aes(x = Comparison, y = DD), colour = '#ff00ff', alpha = 0.5, size = 0.8) + 
    geom_text_repel(data = dplyr::filter(plot.data, logp > -log10(0.05)), aes(label = Gene), parse=TRUE) +
    
    # Labeling points around positive control
    # geom_text_repel(data = dplyr::filter(plot.data, abs(DD) < 0.05 & Comparison < -.8 & !(Gene %in% positive.controls)), aes(label = Gene), parse=TRUE) +
    geom_point(x = pos.control.median.x, y = pos.control.median.y, colour = "#05ba02", alpha = 0.5, size = 0.8) +
    
    # Lighly labeling positive controls
    geom_point(data = dplyr::filter(plot.data, Gene %in% positive.controls),
    colour = "#05ba02", alpha = 0.4) +
    
    # Positive Control Error
    geom_segment(aes(x = -1, y = pos.control.median.y-pos.control.mad.y,
                        xend = -1, yend = pos.control.median.y + pos.control.mad.y), colour = "#05ba02",
                        size = 0.2) +
    geom_segment(aes(x = pos.control.median.x - pos.control.mad.x, 
    y = pos.control.median.y, xend = pos.control.median.x + pos.control.mad.x,
                     yend = pos.control.median.y), colour = "#05ba02", size=0.2) +
    
    
    # Negative Control Error (Delta Delta mad: 0.1565666, LFC mad: 0.2075198)
    geom_segment(aes(x = negative.control.x, y = negative.control.y - 0.1565666,
                     xend = negative.control.x, yend = negative.control.y + 0.1565666), 
                 colour = "#ff00ff", size = 0.15) +
    geom_segment(aes(x = negative.control.x - 0.2075198, y = negative.control.y,
                     xend = negative.control.x + 0.2075198, yend = negative.control.y),
                 colour = "#ff00ff", size = 0.15)



  ggsave(save.path)
}


if (!dir.exists("./output/figures/stars")){
  dir.create("./output/figures/stars", recursive = TRUE)
}

# comparison.data.b6.to.p3.data <- read.table("./output/stars/comparison/p3/STARS-with-Median-LFC-B6.csv",
#                                        header = TRUE, sep = ",")

# 0.00010472300764474696 is the lower limit for the q-value as calculated by STARS

b6.data <- read.table("./output/stars/individual/STARS-with-Median-LFC-B6.csv", header = TRUE, sep = ",")
b6.data[, "Gene.Symbol"] <- paste0("italic('", b6.data[, "Gene.Symbol"], "')")

comparison.data.b6.to.p7.data <- read.table("./output/stars/comparison/p7/STARS-with-Median-LFC-B6.csv",
                                       header = TRUE, sep = ",")
p7.data <- read.table("./output/stars/individual/STARS-with-Median-LFC-P7.csv",
                      header = TRUE, sep = ",")

comparison.data.b6.to.p7.data[, "Gene.Symbol"] <- paste0("italic('", comparison.data.b6.to.p7.data[, "Gene.Symbol"], "')")
p7.data[, "Gene.Symbol"] <- paste0("italic('", p7.data[, "Gene.Symbol"], "')")

comparison.data.b6.to.p7.data.plot <- comparison.data.b6.to.p7.data[, c("Gene.Symbol", "B6", "q.value")]
comparison.data.b6.to.p7.data.plot[,"q.value"][comparison.data.b6.to.p7.data.plot[,"q.value"] == 0] <- 0.00010472300764474696
comparison.data.b6.to.p7.data.plot[, "q.value"] <- -log10(comparison.data.b6.to.p7.data.plot[, "q.value"])
colnames(comparison.data.b6.to.p7.data.plot) <- c("Gene", "Median Normalized LFC Difference", "-log(q-value)")
delta.delta.plot <- comparison.data.b6.to.p7.data.plot[, c("Gene", "Median Normalized LFC Difference", "-log(q-value)")]
delta.delta.plot <- delta.delta.plot[order(delta.delta.plot[, "Gene"]), ]
p7.data <- p7.data[order(p7.data[, "Gene.Symbol"]), ]
delta.delta.plot[, "B6"] <- b6.data[, "B6"]
delta.delta.plot <- delta.delta.plot[order(-delta.delta.plot[, "-log(q-value)"]), ]
colnames(delta.delta.plot) <- c("Gene", "B6[LFC[Normalized]] - P7[LFC[Normalized]]", "-log(q-value)", "B6[LFC[Normalized]]")

delta.delta.plot <- delta.delta.plot[order(-delta.delta.plot[, "-log(q-value)"]), ]

dir.create("./output/figures/stars/experimental/")
create.delta.delta(delta.delta.plot, "In Vivo vs. 2D", "./output/figures/stars/experimental/InVivo2DComparison-Point-Change.jpg", b6.data, p7.data)

nude.data <- read.table("./output/stars/individual/STARS-with-Median-LFC-Nude.csv",
                      header = TRUE, sep = ",")
comparison.data.b6.to.nude.data <- read.table("./output/stars/comparison/nude/STARS-with-Median-LFC-B6.csv",header = TRUE, sep = ",") 
comparison.data.b6.to.nude.data[, "Gene.Symbol"] <- paste0("italic('", comparison.data.b6.to.nude.data[, "Gene.Symbol"], "')")
nude.data[, "Gene.Symbol"] <- paste0("italic('", nude.data[, "Gene.Symbol"], "')")

comparison.data.b6.to.nude.data.plot <- comparison.data.b6.to.nude.data[, c("Gene.Symbol", "B6", "q.value")]
comparison.data.b6.to.nude.data.plot[,"q.value"][comparison.data.b6.to.nude.data.plot[,"q.value"] == 0] <- 0.00010472300764474696
comparison.data.b6.to.nude.data.plot[, "q.value"] <- -log10(comparison.data.b6.to.nude.data.plot[, "q.value"])
colnames(comparison.data.b6.to.nude.data.plot) <- c("Gene", "Median Normalized LFC Difference", "-log(q-value)")
delta.delta.plot <- comparison.data.b6.to.nude.data.plot[, c("Gene", "Median Normalized LFC Difference", "-log(q-value)")]
delta.delta.plot <- delta.delta.plot[order(delta.delta.plot[, "Gene"]), ]
nude.data <- nude.data[order(nude.data[, "Gene.Symbol"]), ]
delta.delta.plot[, "B6"] <- b6.data[, "B6"]
delta.delta.plot <- delta.delta.plot[order(-delta.delta.plot[, "-log(q-value)"]), ]
colnames(delta.delta.plot) <- c("Gene", "B6[LFC[Normalized]] - Nude[LFC[Normalized]]", "-log(q-value)", "B6[LFC[Normalized]]")

delta.delta.plot <- delta.delta.plot[order(-delta.delta.plot[, "-log(q-value)"]), ]

create.delta.delta(delta.delta.plot, "B6 vs. Nude", "./output/figures/stars/B6vsNude-Comparison-Point-Change.jpg", b6.data, nude.data) 
          


comparison.data.nude.to.p7.data <- read.table("./output/stars/comparison/p7/STARS-with-Median-LFC-Nude.csv",
                                       header = TRUE, sep = ",")

comparison.data.nude.to.p7.data[, "Gene.Symbol"] <- paste0("italic('", comparison.data.nude.to.p7.data[, "Gene.Symbol"], "')")
comparison.data.nude.to.p7.data.plot <- comparison.data.nude.to.p7.data[, c("Gene.Symbol", "Nude", "q.value")]
comparison.data.nude.to.p7.data.plot[,"q.value"][comparison.data.nude.to.p7.data.plot[,"q.value"] == 0] <- 0.00010472300764474696
comparison.data.nude.to.p7.data.plot[, "q.value"] <- -log10(comparison.data.nude.to.p7.data.plot[, "q.value"])
colnames(comparison.data.nude.to.p7.data.plot) <- c("Gene", "Median Normalized LFC Difference", "-log(q-value)")
delta.delta.plot <- comparison.data.nude.to.p7.data.plot[, c("Gene", "Median Normalized LFC Difference", "-log(q-value)")]
delta.delta.plot <- delta.delta.plot[order(delta.delta.plot[, "Gene"]), ]
p7.data <- p7.data[order(p7.data[, "Gene.Symbol"]), ]
delta.delta.plot[, "Nude"] <- nude.data[, "Nude"]
delta.delta.plot <- delta.delta.plot[order(-delta.delta.plot[, "-log(q-value)"]), ]
colnames(delta.delta.plot) <- c("Gene", "Nude[LFC[Normalized]] - P7[LFC[Normalized]]", "-log(q-value)", "Nude[LFC[Normalized]]")

delta.delta.plot <- delta.delta.plot[order(-delta.delta.plot[, "-log(q-value)"]), ]

create.delta.delta(delta.delta.plot, "In Vivo (Nude) vs. 2D (P7)", "./output/figures/stars/Nude-2D-Comparison-Point-Change.jpg", nude.data, p7.data)


p3.data <- read.table("./output/stars/individual/STARS-with-Median-LFC-P3.csv",
                      header = TRUE, sep = ",")
threeD.data <- read.table("./output/stars/individual/STARS-with-Median-LFC-ThreeD.csv",
header = TRUE, sep = ",")
comparison.data.3d.to.p3.data <- read.table("./output/stars/comparison/p3/STARS-with-Median-LFC-ThreeD.csv",
                                       header = TRUE, sep = ",")
comparison.data.3d.to.p3.data[, "Gene.Symbol"] <- paste0("italic('", comparison.data.3d.to.p3.data[, "Gene.Symbol"], "')")
comparison.data.3d.to.p3.data.plot <- comparison.data.3d.to.p3.data[, c("Gene.Symbol", "ThreeD", "q.value")]
comparison.data.3d.to.p3.data.plot[,"q.value"][comparison.data.3d.to.p3.data.plot[,"q.value"] == 0] <- 0.00010472300764474696
comparison.data.3d.to.p3.data.plot[, "q.value"] <- -log10(comparison.data.3d.to.p3.data.plot[, "q.value"])
colnames(comparison.data.3d.to.p3.data.plot) <- c("Gene", "Median Normalized LFC Difference", "-log(q-value)")
delta.delta.plot <- comparison.data.3d.to.p3.data.plot[, c("Gene", "Median Normalized LFC Difference", "-log(q-value)")]
delta.delta.plot <- delta.delta.plot[order(delta.delta.plot[, "Gene"]), ]
p3.data <- p3.data[order(p3.data[, "Gene.Symbol"]), ]
delta.delta.plot[, "ThreeD"] <- threeD.data[, "ThreeD"]
delta.delta.plot <- delta.delta.plot[order(-delta.delta.plot[, "-log(q-value)"]), ]
colnames(delta.delta.plot) <- c("Gene", "Three-D[LFC[Normalized]] - P3[LFC[Normalized]]", "-log(q-value)", " Three-D[LFC[Normalized]]")

delta.delta.plot <- delta.delta.plot[order(-delta.delta.plot[, "-log(q-value)"]), ]

create.delta.delta(delta.delta.plot, "3D vs. 2D (P3)", "./output/figures/stars/3D-2D-Comparison-Point-Change.jpg", threeD.data, p3.data)

