# This script make plots the key outlier region on chromosome8
# This refers to POS = 20550000 - 20750000
#
# Library -----------------------------------------------------------------
library(data.table)
library(tidyverse)
library(cowplot)
library(egg)
# library(patchwork)
# library(ggpubr)
# Define input ------------------------------------------------------------
XtXfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_summary_pi_xtx_SNPs.out"
C2file <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_IS_C2_GEA_summary_contrast_CON_001_SNPs.out"
C2CON004file <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_IS_C2_GEA_summary_contrast_CON_004_SNPs.out"

## Output plots -----------------------------------------------------
outpng <- paste("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/manhattanplot_XtX_C2.png", sep = "")

chr_order <- c("1", "1A", "2", "3", "4", "4A", "5a", "5b", "5c", "6", "7", "8", "9", 
               "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", 
               "20", "21", "22", "23", "24", "25", "26", "27", "29", 
               "30", "31", "32")
# Read in input files ---------------------------------------------------------
## Read XtX file --------------------------------------------------------------
dt_XtX <- fread(file = XtXfile, header = T) |>
  rename(value = M_XtX)

dt_XtX$MRK <- dt_XtX$MRKALL
dt_XtX$MRKALL <- NULL
dt_XtX$chr <- gsub(pattern = "^Superscaffold_chr", replacement = "", x = dt_XtX$chr)
dt_XtX$chr <- factor(x = dt_XtX$chr, levels = chr_order)
dt_meta_SNP <- dt_XtX[,1:4]

## Create cumulative bp length table --------------------------------------
WGS_cum <- dt_meta_SNP |>
  group_by(chr) |>
  summarise(max_bp = max(POS)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(chr, bp_add)

## Merge SNP meta table with cumulated bp length table --------------------
dt_meta_SNP <- dt_meta_SNP |>
  inner_join(WGS_cum, by = c("chr")) |>
  mutate(bp_cum = POS + bp_add)

## Create a table for labels ----------------------------------------------
axis_set <- dt_meta_SNP |>
  group_by(chr) |>
  summarize(center = mean(bp_cum))

## Merge table with XtX with SNP meta table -------------------------------
dt_XtX <- merge(dt_XtX, dt_meta_SNP[,c("ID", "bp_cum")], by = "ID")
dt_XtX$data <- "XtX"

## Create table with colours ----------------------------------------------
axis_set$col = rep(c("#276FBF", "#183059"),
                   length.out = unique(length(axis_set$chr)))

dt_XtX <- merge(dt_XtX, axis_set[,c("chr", "col")], by = "chr") 

## Define outliers --------------------------------------------------------
quant_thres <- 0.99999
path2XtXscript <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_XtX_threshold.R"
PODxtx <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/subset_001/G.btapods_summary_pi_xtx.out"
commandArgs <- function(...) c(PODxtx,quant_thres)
XtXthresh <-  as.numeric(capture.output(source(path2XtXscript)))
dt_XtX$col[dt_XtX$value > XtXthresh] <- "red"

## Read in C2 (CON001) file ------------------------------------------------------------
dt_C2 <- fread(file = C2file, header = T) |>
  rename(value = M_C2)

dt_C2$chr <- gsub(pattern = "^Superscaffold_chr", replacement = "", x = dt_C2$chr)
dt_C2$chr <- factor(x = dt_C2$chr, levels = chr_order)
dt_C2 <- merge(dt_C2, dt_meta_SNP[,c("ID", "bp_cum")], by = "ID")
dt_C2 <- merge(dt_C2, axis_set[,c("chr", "col")], by = "chr") 

path2C2script <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_C2_threshold.R"
PODC2 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/subset_001/G.btapods_IS_C2_GEA_summary_contrast.out"
quant_thres <- 0.99999
commandArgs <- function(...) c(PODC2, "M_C2", 1, quant_thres)
C2thresh <-  as.numeric(capture.output(source(path2C2script)))
dt_C2$col[dt_C2$value > C2thresh] <- "red"

dt_XtX$data <- "XtX"
dt_C2$data <- "C2"



## Read in C2 (CON004) file ------------------------------------------------------------
dt_C2CON004 <- fread(file = C2CON004file, header = T) |>
  rename(value = M_C2)

dt_C2CON004$chr <- gsub(pattern = "^Superscaffold_chr", replacement = "", x = dt_C2CON004$chr)
dt_C2CON004$chr <- factor(x = dt_C2CON004$chr, levels = chr_order)
dt_C2CON004 <- merge(dt_C2CON004, dt_meta_SNP[,c("ID", "bp_cum")], by = "ID")
dt_C2CON004 <- merge(dt_C2CON004, axis_set[,c("chr", "col")], by = "chr") 

path2C2script <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_C2_threshold.R"
PODC2CON004 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/subset_001/G.btapods_IS_C2_GEA_summary_contrast.out"
quant_thres <- 0.99999
commandArgs <- function(...) c(PODC2CON004, "M_C2", 4, quant_thres)
C2thresh <-  as.numeric(capture.output(source(path2C2script)))
dt_C2CON004$col[dt_C2CON004$value > C2thresh] <- "red"

# Make plots ------------------------------------------------------------------
theme_plots <- theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(face = "bold", colour = "black", size = 15),
        axis.title.x = element_text(face = "bold", colour = "black"),
        panel.grid = element_blank())

theme_empty <- theme_plots + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))


plot_val <- function(dt, axis_set, ylab){
  p <- ggplot() +
    geom_point(dt, mapping = aes(
      x = bp_cum, y = value), alpha = 0.75,
      size = 0.25, col = dt$col) +
    scale_x_continuous(expand = c(0,0),
                       label = axis_set$chr,
                       breaks = axis_set$center
    ) +
    scale_y_continuous(expand = c(0, 0)) + 
    labs(
      x = "Chromosome",
      y = ylab
    ) + 
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
          plot.margin = unit(c(0, 5.5, 0, 5.5), "pt"),
          axis.title.y = element_text(size = 13))
  return(p)
}
## XtX --------------------------------------------------------------------
# p2 <- plot_val(sample_n(dt_XtX, 10000), axis_set, ylab = "XtX")
p2 <- plot_val(dt_XtX, axis_set, ylab = "XtX")

# p2

## C2 ---------------------------------------------------------------------
# p3 <- plot_val(sample_n(dt_C2, 10000), axis_set, ylab = expression(C[2]))
p3 <- plot_val(dt_C2, axis_set, ylab = expression(C[2]))
## CON004
p3 <- plot_val(dt_C2CON004, axis_set, ylab = expression(C[2]))

png(outpng, width = 8.3, height = 11.7/2, units = "in", res = 600)
ggarrange(p2 + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank()), 
          p3 + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank()), 
          p4, ncol = 1, labels = c("A)", "B)", "C)"),
          label.args = list(gp=grid::gpar(font=2)))
dev.off()


