# This script make plots the key outlier region on chromosome8
# This refers to POS = 20550000 - 20750000
#
# Library -----------------------------------------------------------------
library(data.table)
library(tidyverse)
library(cowplot)
# library(patchwork)
# library(ggpubr)
# Define input ------------------------------------------------------------
XtXfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_summary_pi_xtx_SNPs.out"
C2file <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_IS_C2_GEA_summary_contrast_CON_001_SNPs.out"
IHSfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/WGS_IHS.txt"
depthfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/vcftools_sum_stats/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.ldepth.mean"
AFfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/MAF_per_pop_chr8.out"
## Annotation files -------------------------------------------------------
Genefile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/outlier_region/chr8_20500000_20800000_genes.txt"
genenamefile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/outlier_region/chr8_20500000_20800000_genes_idwithnames.txt"
## Position parameters ----------------------------------------------------
chrno <- "chr8"
chrval <- paste("Superscaffold_", chrno, sep = "")
start_pos <- 20580000
end_pos <- 20750000

x_axis_brks <- seq(20600000, end_pos-50000, by = 50000)
x_axis_brks2 <- sort(c(x_axis_brks, 20678727, 20687524))


## Output plots -----------------------------------------------------
outpng <- paste("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/combined_", chrno, "_", start_pos, "_", end_pos, "_XtX_C2_IHS_AF_withgenetracksv2.png", sep = "")

# Read in input files ---------------------------------------------------------
## Read XtX file --------------------------------------------------------------
dt_XtX <- fread(file = XtXfile, header = T) |>
  filter(chr == chrval) |>
  filter(POS >= start_pos & POS <= end_pos) |> 
  rename(value = M_XtX)

path2script <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_XtX_threshold.R"
PODxtx <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/subset_001/G.btapods_summary_pi_xtx.out"
quant_thres <- 0.99999
commandArgs <- function(...) c(PODxtx,quant_thres)
XtXthresh <-  as.numeric(capture.output(source(path2script)))
dt_XtX$colour[dt_XtX$value > XtXthresh] <- "red"
dt_XtX$colour[is.na(dt_XtX$colour)] <- "#183059"

## Read in C2 file ------------------------------------------------------------
dt_C2 <- fread(file = C2file, header = T) |>
  filter(chr == chrval) |>
  filter(POS >= start_pos & POS <= end_pos) |>
  rename(value = M_C2)

path2script <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_C2_threshold.R"
PODC2 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/subset_001/G.btapods_IS_C2_GEA_summary_contrast.out"
quant_thres <- 0.99999
commandArgs <- function(...) c(PODC2, "M_C2", 1, quant_thres)
C2thresh <-  as.numeric(capture.output(source(path2script)))
dt_C2$colour[dt_C2$value > C2thresh] <- "red"
dt_C2$colour[is.na(dt_C2$colour)] <- "#183059"
## Read IHS file --------------------------------------------------------------
dt_IHS <- fread(file = IHSfile, header = T) |>
  filter(CHR == chrval) |>
  filter(POSITION >= start_pos & POSITION <= end_pos) |>
  rename(value = LOGPVALUE) |>
  rename(POS = POSITION) |>
  rename(chr = CHR)

dt_IHS$colour[dt_IHS$value > 6] <- "red"
dt_IHS$colour[is.na(dt_IHS$colour)] <- "#183059"
dt_XtX$data <- "XtX"
dt_C2$data <- "C2"
dt_IHS$data <- "iHS"
## Creat dummy gene track table for plotting ------------------------------
dt_genetracks <- data.table(chr = chrval,
                            POS = NA,
                            value = c(0.9,2.4),
                            data = "",
                            colour = NA)

## Read genes -------------------------------------------------------------
dt_genes <- fread(Genefile, 
                  header = F)

dt_genes <- data.table(CHR = dt_genes$V1,
                       startPOS = dt_genes$V4,
                       endPOS = dt_genes$V5,
                       ID = gsub(pattern = "ID=", replacement = "", dt_genes$V9))

dt_genes <- dt_genes |> 
  filter(endPOS > start_pos & endPOS < end_pos)
dt_genes$starty <- c(1, 1, 1, 1, 1, 1, 1)
dt_genes$endy <- dt_genes$starty + 1
dt_genes$label <- 1:7
dt_genes$fill <- "chartreuse4"
dt_genes$fill[dt_genes$label == 6] <- "orange"
dt_genes$ylabpos <- 2
# rbindlist(list(dt_genes, c("Superscaffold_chr8", 20678727, 20687524, "repeat_482137", 1, 2, "grey", "", 2)))
# dt_rep <- data.table(CHR = "Superscaffold_chr8", 
#                      startPOS = 20678727, 
#                      endPOS = 20687524, 
#                      ID = "repeat_482137", 
#                      starty = 1, 
#                      endy = 2, 
#                      label = "", 
#                      fill = "grey", 
#                      ylabpos = 2)
# dt_rep_comb <- rbind(dt_rep |> mutate(data = ""))
# dt_rep_comb$data <- factor(dt_rep_comb$data, levels = c("", "XtX", "C2", "iHS"))
# 
# # dt_genes <- rbindlist(list(dt_genes, dt_rep))
# dt_genes_comb <- rbind(dt_genes |> mutate(data = ""))
# dt_genes_comb$data <- factor(dt_genes_comb$data, levels = c("", "XtX", "C2", "iHS"))


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
## Gene tracks -------------------------------------------------------------
p1 <- ggplot() + 
  geom_rect(dt_rep, mappin = aes(xmin=startPOS/1000000, xmax=endPOS/1000000,
                                      ymin= starty, ymax= endy), fill = dt_rep$fill, col = NA) +
  geom_rect(dt_genes, mappin = aes(xmin=startPOS/1000000, xmax=endPOS/1000000,
                                        ymin= starty, ymax= endy), fill = dt_genes$fill, col = NA) +
  geom_text(dt_genes, mapping = aes(x = (endPOS/1000000)+0.0015, y = ylabpos, label = label)) + 
  xlab("Position (Mb)") +
  scale_y_continuous(limits = c(0.9, 2.3), expand = c(0,0)) +
  scale_x_continuous(limits = c(start_pos/1000000, end_pos/1000000), expand = c(0,0), position = "top",  breaks = x_axis_brks/1000000) +
  theme_plots + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 5.5, 2.5, 5.5), "points"))

# p1

plot_val <- function(dt, ylab){
  p <- ggplot(data = dt) +
    geom_point(mapping = aes(x = POS/1000000, y = value), colour = dt$colour, size = 0.5) + 
    scale_x_continuous(expand = c(0,0)) +
    ylab(ylab) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.title.y = element_text(size = 13))
  return(p)

}
## XtX --------------------------------------------------------------------
p2 <- plot_val(dt_XtX, ylab = "XtX")
# p2

## C2 ---------------------------------------------------------------------
p3 <- plot_val(dt_C2, ylab = expression(C[2]))
# p3


## iHS --------------------------------------------------------------------
iHSlab <- expression(paste(-log[10](italic("p")), " iHS", sep = ""))
p4 <- ggplot(data = dt_IHS) +
  geom_point(mapping = aes(x = POS/1000000, y = value), colour = dt_IHS$colour, size = 0.5) + 
  scale_x_continuous(expand = c(0,0), breaks = x_axis_brks/1000000, position = "top", 
                     sec.axis =  dup_axis(name = NULL,
                                          breaks = x_axis_brks2/1000000, labels = NULL)) +
  ylab(iHSlab) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.y = element_text(size = 13))
  

pcomb <- cowplot::plot_grid(p1,
                   p2,
                   p3, 
                   p4,
                   align = "v", 
                   ncol = 1, 
                   # axis = "lr", 
                   rel_heights = c(1.7, 4,4,4), labels = c("A)", "B)", "C)", "D)"))

## Read in AF file ------------------------------------------------------------
col_ls <- c("CHR", "SNP", "CLST", "A1", "A2", "MAF", "MAC", "NCHROBS")
CLSTorder <- c("Madhya_Pradesh", "Tamil_Nadu", "Maharashtra", "Maharashtra_subpopulation_A", "Fiji", "Melbourne", "Cairns", "Napier", "Leigh", "Great_Barrier_Island", "South_Africa")
dt_AF <- fread(input = AFfile, header = F, 
               col.names = col_ls)

dt_AF <- dt_AF[dt_AF$SNP %in% dt_XtX$ID, ]
dt_AF$NAT_INV <- "Invasive"
dt_AF$NAT_INV[dt_AF$CLST %in% c("Madhya_Pradesh", "Maharashtra", "Maharashtra_subpopulation_A", "Tamil_Nadu")] <- "Native"
dt_AF$CLST <- factor(dt_AF$CLST, levels = CLSTorder)
dt_AF$POS <- sapply(strsplit(x = dt_AF$SNP, split = "_"), "[[" , 3)
dt_AF$POS <- as.numeric(dt_AF$POS)
dt_AF <- dt_AF %>% 
  filter(POS >= start_pos & POS <= end_pos)
dt_AF$ID <- dt_AF$SNP

ls_SNP_convert <- dt_AF$SNP[dt_AF$MAF > 0.5 & dt_AF$CLST == "Madhya_Pradesh"]


dt_AF <- dt_AF %>%
  mutate(MAF_STD = case_when(SNP %in% ls_SNP_convert ~ 1 - MAF, 
                             TRUE ~ MAF))

dt_AF$SNPPOS <- as.character(dt_AF$POS)
poplab <- c("Madhya\n Pradesh", "Tamil Nadu", "Maharashtra", "Maharashtra\nsubpop. A", "Fiji", "Melbourne", "Cairns", "Napier", "Leigh", "Great\nBarrier Isl.", "South\nAfrica")
brk <- c("Madhya_Pradesh", "Tamil_Nadu", "Maharashtra", "Maharashtra_subpopulation_A", "Fiji", "Melbourne", "Cairns", "Napier", "Leigh", "Great_Barrier_Island", "South_Africa")
dt_POS <- data.table(POS = sort(unique(dt_AF$POS)))
dt_POS$Rank <- rank(dt_POS$POS, ties.method = "min")
dt_AF <- merge(dt_POS, dt_AF, by = c("POS"))
dt_poplab <- data.table(CLST = brk, label = poplab)
dt_poplab$order <- c(seq(1, dim(dt_poplab)[1]))
dt_AF <- merge(dt_poplab, dt_AF, by = c("CLST"))

# Calculate key positions in ranks ----------------------------------------
val <- c()
for (i in x_axis_brks2){
  val <- c(val, mean(c(unique(dt_AF$Rank[dt_AF$POS == min(dt_AF$POS[dt_AF$POS >= i])]),
                       unique(dt_AF$Rank[dt_AF$POS == max(dt_AF$POS[dt_AF$POS <= i])]))))
}

labx <- min(dt_AF$Rank) + 0.02*diff(range(dt_AF$Rank))

phm <- ggplot(dt_AF, aes(x = Rank, y = order, fill = MAF_STD)) +
  geom_tile() +
  scale_x_discrete(expand = c(0,0)) +
  # geom_text(y = "South Africa", x = labx, label = "D") +
  # scale_y_discrete() +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_fill_gradient2(name = "Allele Frequency", low = "red", mid = "white", high = "blue", midpoint = 0.5, breaks = c(0, 0.5, 1)) + 
  scale_y_continuous(breaks = dt_poplab$order, labels = dt_poplab$label, expand = c(0,0)) + 
  # ggtitle("Allele frequency vs population (standardised to Madhya Pradesh)") + 
  theme(legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.background = element_rect(colour = "black"),
        legend.direction = "horizontal",
        legend.margin = margin(5, 12, 5, 12),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(colour = "black", face = "bold", size = 13),
        axis.text.x = element_blank(),
        plot.margin = unit(c(9, 7, 7, 7), "points")) + 
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# phm
phm <- phm + geom_vline(xintercept = val,
                        linetype = "dashed", colour = "black", size = 0.5)


# Combine plots
allset2plots <- cowplot::plot_grid(pcomb,
                   phm, 
                   # align = "v", 
                   ncol = 1, 
                   rel_heights = c(13, 9),
                   labels = c("", "E)"))


png(outpng, width = 8.3, height = 11, units = "in", res = 600)
allset2plots
dev.off()


