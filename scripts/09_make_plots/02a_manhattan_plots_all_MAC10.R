
# Define input arguments --------------------------------------------#####

# XtX
XtX <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_summary_pi_xtx_SNPs.out"
PODxtx <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/subset_001/G.btapods_summary_pi_xtx.out"
path2XtXscript <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_XtX_threshold.R"

# Library ---------------------------------------------------------------------
library(tidyverse)
library(gridExtra)
library(data.table)

# Define functions --------------------------------------------------------
chr_order <- c("1", "1A", "2", "3", "4", "4A", "5a", "5b", "5c", "6", "7", "8", "9", 
               "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", 
               "20", "21", "22", "23", "24", "25", "26", "27", "29", 
               "30", "31", "32")

# XtX data ----------------------------------------------------------------
## Read XtX tables --------------------------------------------------------
dt_XtX <- fread(file = XtX, header = T) |> 
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
commandArgs <- function(...) c(PODxtx,quant_thres)
XtXthresh <-  as.numeric(capture.output(source(path2XtXscript)))
dt_XtX$col[dt_XtX$value > XtXthresh] <- "red"
# dt_XtX$col[is.na(dt_XtX$col)] <- "black"

# C2 ----------------------------------------------------------------------
path2C2script <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_C2_threshold.R"
PODC2 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/subset_001/G.btapods_IS_C2_GEA_summary_contrast.out"
for (i in 1:7){
  C2 <- paste("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_IS_C2_GEA_summary_contrast_CON_00", i, "_SNPs.out", sep = "")
  dt_C2 <- fread(file = C2, header = T) |> 
    rename(value = M_C2)
  dt_C2$chr <- gsub(pattern = "^Superscaffold_chr", replacement = "", x = dt_C2$chr)
  dt_C2$chr <- factor(x = dt_C2$chr, levels = chr_order)
  dt_C2 <- merge(dt_C2, dt_meta_SNP[,c("ID", "bp_cum")], by = "ID")
  dt_C2$data <- paste("CON00", i, sep = "")
  dt_C2 <- merge(dt_C2, axis_set[,c("chr", "col")], by = "chr") 
  
  quant_thres <- 0.99999
  commandArgs <- function(...) c(PODC2, "M_C2", i, quant_thres)
  C2thresh <-  as.numeric(capture.output(source(path2C2script)))
  dt_C2$col[dt_C2$value > C2thresh] <- "red"
  # dt_C2$col[is.na(dt_C2$col)] <- "black"
  if (i == 1) {
    dt_C2_comb <- dt_C2
  } else {
    dt_C2_comb <- rbind(dt_C2_comb, dt_C2)
  }
}

# iHS ---------------------------------------------------------------------
IHSfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/WGS_IHS.txt"
dt_IHS <- fread(file = IHSfile, header = T) |>
  rename(value = LOGPVALUE) |>
  rename(POS = POSITION) |>
  rename(chr = CHR)
dt_IHS$chr <- gsub(pattern = "^Superscaffold_chr", replacement = "", x = dt_IHS$chr)
dt_IHS$chr <- factor(x = dt_IHS$chr, levels = chr_order)
dt_IHS <- cbind(dt_IHS, dt_meta_SNP[,c("bp_cum")]) # Only possible because they are in the same order. If not, see code below
# dt_IHS <- merge(dt_IHS, dt_meta_SNP[,c("chr", "POS", "bp_cum")], by = c("chr", "POS"))
dt_IHS$data <- "iHS"
dt_IHS <- merge(dt_IHS, axis_set[,c("chr", "col")], by = "chr") 
dt_IHS$col[dt_IHS$value > 6] <- "red"

dt_outlier_region <- fread(input = "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/combined_outliers/outliers_n5_summary.csv", header = T)
dt_outlier_region <- dt_outlier_region %>%
  rename(data = outlier_statistics)
dt_outlier_region$data <- factor(dt_outlier_region$data, levels = c("XtX", "iHS", "CON001", "CON002", "CON003", "CON004", "CON005", "CON006", "CON007"))
dt_outlier_region$chr <- gsub(pattern = "^Superscaffold_chr", replacement = "", x = dt_outlier_region$chr)
dt_outlier_region <- merge(dt_outlier_region, dt_meta_SNP[,c("chr", "POS", "bp_cum")], by.x = c("chr", "start_pos"), by.y = c("chr", "POS"))
dt_outlier_region <- dt_outlier_region %>% 
  rename(startbpcum = bp_cum)
dt_outlier_region <- merge(dt_outlier_region, dt_meta_SNP[,c("chr", "POS", "bp_cum")], by.x = c("chr", "end_pos"), by.y = c("chr", "POS"))
dt_outlier_region <- dt_outlier_region %>% 
  rename(endbpcum = bp_cum)

rm(dt_meta_SNP)
gc()
# Combine all plots -------------------------------------------------------

col_common <- intersect(intersect(colnames(dt_C2_comb), colnames(dt_XtX)), colnames(dt_IHS))
dt_merge <- rbindlist(list(dt_XtX[, ..col_common],
                       dt_C2_comb[, ..col_common],
                       dt_IHS[, ..col_common]))

dt_merge$data <- factor(dt_merge$data, levels = c("XtX", "iHS", "CON001", "CON002", "CON003", "CON004", "CON005", "CON006", "CON007"))

theme_manh <- theme_bw() + 
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

manhplot <- ggplot() +
  geom_point(dt_merge, mapping = aes(
    x = bp_cum, y = value), alpha = 0.75,
    size = 0.25, col = dt_merge$col) +
  scale_x_continuous(expand = c(0,0),
    label = axis_set$chr,
    breaks = axis_set$center
  ) +
  scale_y_continuous(expand = c(0, 0)) + #, limits = c(0, ylim)) +
  # geom_hline(dt.pod.comb.thresh, mapping = aes(yintercept = POD_quantiles, col = lab_comp)) +
  # scale_color_manual(name = "POD quantiles", values = pod.col) +
  # geom_point(dt_XtX_outliers, mapping = aes(x = bp_cum, y = M_XtX, fill = label), col = NA, alpha = 0.75, size = 0.25, shape = 21) +
  # scale_fill_manual(name = "POD threshold", values = quant_col) +
  labs(
    x = "Chromosome",
    y = NULL
  ) + 
  # ggtitle("XtX (M_XtX); DArT_ALL vs DArT_WGSsamp vs WGS") + 
  theme_manh +
  facet_grid(rows = vars(data), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.placement = "outside")

outpng <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/manhattanplot_allv1.png"
png(outpng, width = 8.3, height = 11.7, units = "in", res = 600)
manhplot
dev.off()

outpng <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/manhattanplot_allv2.png"
png(outpng, width = 8.3, height = 11.7, units = "in", res = 600)
manhplot +
  geom_rect(dt_outlier_region, mappin = aes(xmin=startbpcum/1000000, xmax=endbpcum/1000000,
                                            ymin= -Inf, ymax= Inf), fill = "red", alpha = 0.4, col = NA)
dev.off()


# Just XtX, iHS and CON001 ------------------------------------------------
dt_merge <- dt_merge %>% # Overwriting the object to save some memory 
  filter(data %in% c("XtX", "iHS", "CON001"))
dt_merge$data <- factor(dt_merge$data, levels = c("XtX", "iHS", "CON001"))

dt_outlier_region <- dt_outlier_region %>% # Overwriting the object to save some memory 
  filter(data %in% c("XtX", "iHS", "CON001"))
dt_outlier_region$data <- factor(dt_outlier_region$data, levels = c("XtX", "iHS", "CON001"))

manhplot <- ggplot() +
  geom_point(dt_merge, mapping = aes(
    x = bp_cum, y = value), alpha = 0.75,
    size = 0.25, col = dt_merge$col) +
  scale_x_continuous(expand = c(0,0),
                     label = axis_set$chr,
                     breaks = axis_set$center
  ) +
  scale_y_continuous(expand = c(0, 0)) + #, limits = c(0, ylim)) +
  labs(
    x = "Chromosome",
    y = NULL
  ) + 
  # ggtitle("XtX (M_XtX); DArT_ALL vs DArT_WGSsamp vs WGS") + 
  theme_manh +
  facet_grid(rows = vars(data), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.placement = "outside")

outpng <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/manhattanplot_XtX_iHS_CON001_v1.png"
png(outpng, width = 8.3, height = 5, units = "in", res = 600)
manhplot
dev.off()

outpng <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/manhattanplot_XtX_iHS_CON001_v2.png"
png(outpng, width = 8.3, height = 5, units = "in", res = 600)
manhplot +
  geom_rect(dt_outlier_region, mappin = aes(xmin=startbpcum/1000000, xmax=endbpcum/1000000,
                                            ymin= -Inf, ymax= Inf), fill = "red", alpha = 0.4, col = NA)
dev.off()
