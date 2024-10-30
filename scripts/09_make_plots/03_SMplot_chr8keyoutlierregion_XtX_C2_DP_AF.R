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
IHSfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/chr8_IHS.txt"
depthfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/vcftools_sum_stats/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.ldepth.mean"
AFfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/MAF_per_pop_chr8.out"
## Annotation files -------------------------------------------------------
TEfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/TE/TE_chr8_20500000_20800000.txt"
SVfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/SV/SV_chr8_20500000_20800000.txt"
VEPfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/output/VEPout_chr8_20500000_20800000_missense_start_stop.txt" 
Genefile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/outlier_region/chr8_20500000_20800000_genes.txt"
Repeatfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/RepeatMasker/Repeats_chr8_20500000_20800000.txt"
## Position parameters ----------------------------------------------------
chrno <- "chr8"
chrval <- paste("Superscaffold_", chrno, sep = "")
start_pos <- 20580000
end_pos <- 20750000

x_axis_brks <- seq(20600000, end_pos, by = 50000)
x_axis_brks2 <- sort(c(x_axis_brks, 20678727, 20687524))


## Output plots -----------------------------------------------------
outpng <- paste("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/combined_", chrno, "_", start_pos, "_", end_pos, "_XtX_C2_IHS_DP_AF.png", sep = "")

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

## Read SNP depth file --------------------------------------------------------
dt_DP <- fread(depthfile, header = T) |>
  filter(CHROM == chrval) |> 
  filter(POS >= start_pos & POS <= end_pos) |>
  rename(value = MEAN_DEPTH) |>
  rename(chr = CHROM)

dt_DP$colour <- "black"

dt_XtX$data <- "XtX"
dt_C2$data <- "C2"
dt_IHS$data <- "IHS"
dt_DP$data <- "DP"

## Combine data tables --------------------------------------------------------
col2keep <- c("chr", "POS", "value", "data", "colour")
dt_comb <- rbind(dt_XtX[, ..col2keep], 
                 dt_C2[, ..col2keep], 
                 dt_IHS[, ..col2keep], 
                 dt_DP[, ..col2keep])

dt_comb$data <- factor(dt_comb$data, levels = c("XtX", "C2", "IHS", "DP"))

## Create table for adding labels ---------------------------------------------
# labrat <- 0.04
# dt_lab <- dt_comb |> 
#   group_by(data) |>
#   summarise(x = min(POS) + labrat*diff(range(POS)),
#             y = max(value) - 0.1*diff(range(value))) |>
#   mutate(label = paste0(c("i", "ii", "iii", "iv"),")"))

# Read in genes/TEs/etc -------------------------------------------------------
## Read genes ----------------------------------------------------------------- 
dt_genes <- fread(Genefile, 
                  header = F)

dt_genes <- data.table(CHR = dt_genes$V1,
                       startPOS = dt_genes$V4,
                       endPOS = dt_genes$V5,
                       ID = gsub(pattern = "ID=", replacement = "", dt_genes$V9))
dt_genes <- dt_genes |> 
  filter(endPOS > start_pos & endPOS < end_pos)
dt_genes_comb <- rbind(dt_genes |> mutate(data = "XtX"),
                       dt_genes |> mutate(data = "C2"),
                       dt_genes |> mutate(data = "IHS"))
dt_genes_comb$data <- factor(dt_genes_comb$data, levels = c("XtX", "C2", "IHS", "DP"))
dt_genes_comb$dtype <- "Gene"

## Read TEs -------------------------------------------------------------------
dt_TE <- fread(TEfile, sep = "\t")

dt_TE <- data.table(CHR = dt_TE$V1,
                    POS = dt_TE$V2) |> 
  filter(POS >= start_pos & POS <= end_pos)

## Read SVs -------------------------------------------------------------------
dt_SV <- fread(SVfile,
               header = T)
colnames(dt_SV)[1] <- "CHROM"
# colnames(dt_SV) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")

inls <- c("<DEL>", "<DUP>", "<DUP:TANDEM>", "<INV>")
dt_SV <- dt_SV |>
  mutate(dtype = case_when(ALT == "<DEL>" ~ "SV: Deletion",
         ALT == "<DUP>" ~ "SV: Duplicate",
         ALT == "<DUP:TANDEM>" ~ "SV: Tandem duplicate",
         ALT == "<INV>" ~ "SV: Inversion",
         TRUE ~ "SV: Insertion/Deletion?"))


## Read VEP missense mutations ------------------------------------------------
dt_missense <- fread(VEPfile, header = F, col.names = c("ID", "gene_ID", "mutation_type"))

dt_missense <- dt_XtX[dt_XtX$ID %in% dt_missense$ID,]

## Combine TE and missense files ----------------------------------------------
dt_missense$dtype <- "missense mutation"
dt_TE$dtype <- "TE variant"
# dt_SV$dtype <- "SV variant"
dt_vline <- rbind(dt_TE[,c("POS", "dtype")], 
                  dt_missense[,c("POS", "dtype")],
                  dt_SV[,c("POS", "dtype")])

dt_vline_comb <- rbind(dt_vline |> mutate(data = "XtX"),
      dt_vline |> mutate(data = "C2"),
      dt_vline |> mutate(data = "IHS"))

dt_vline_comb$data <- factor(dt_vline_comb$data, levels = c("XtX", "C2", "IHS", "DP"))

## Read RepeatMasker ----------------------------------------------------------
dt_RE <- read.fwf(file = Repeatfile, 
                  widths = c(6,7,5,5,25,15,10,12,2,65,17,8,7,8,7,2), header = FALSE)

dt_RE <- data.table(CHR = dt_RE$V5,
                    startpos = dt_RE$V6,
                    endpos = dt_RE$V7,
                    type1 = dt_RE$V10,
                    type2 = dt_RE$V11,
                    ID = dt_RE$V15)

dt_RE <- dt_RE |> 
  filter(startpos > start_pos & endpos < end_pos)
dt_RE$data <- "DP"
dt_RE$data <- factor(dt_RE$data, levels = c("XtX", "C2", "IHS", "DP"))
dt_RE$dtype <- "Repeat (RepeatMasker)"

## Make plot of XtX, C2, IHS, DP ----------------------------------------------
p <- ggplot(data = dt_comb) +
  geom_point(mapping = aes(x = POS/1000000, y = value), size = 0.5, colour = dt_comb$colour) +  
  scale_x_continuous(expand = c(0,0), position = "top") +
  # scale_x_continuous(expand = c(0,0), breaks = c(20.55, 20.6, 20.65, 20.7, 20.75), position = "top") +
  facet_wrap(~data, scales = "free_y", nrow = 4, 
             strip.position = "left", 
             # labeller = as_labeller(c(A = "Currents (A)", V = "Voltage (V)") ) 
             )  +
  theme_bw() +
  ylab(NULL) +
  xlab("Position (Mb)") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(face = "bold", colour = "black", size = 15),
        axis.title.x = element_text(face = "bold", colour = "black"),
        panel.grid = element_blank())

# Add gene tracks -------------------------------------------------------------
p <- p +
  geom_rect(dt_genes_comb, mappin = aes(xmin=startPOS/1000000, xmax=endPOS/1000000,
                                   ymin= -Inf, ymax= Inf, fill = dtype), alpha = 0.2, col = NA) +
  geom_vline(dt_vline_comb, mapping = aes(xintercept = POS/1000000, col = dtype, lty = dtype)) +
  theme() +
  geom_rect(dt_RE, mappin = aes(xmin=startpos/1000000, xmax=endpos/1000000,
                                 ymin= -Inf, ymax= Inf, fill = dtype), alpha = 0.4, col = NA) +
  # geom_text(data = dt_lab, aes(label = label, x = x/1000000, y = y), 
  #           color = "black", fontface = "bold", size = 15/.pt) +
  scale_fill_manual(breaks = c("Gene", "Repeat (RepeatMasker)"), values = c("#F8766D", "#00BFC4")) +
  theme(legend.box = "horizontal") +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        # legend.justification = c(0.99, 0.99),
        legend.title = element_blank(),
        # legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "lines"),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1))



## Add second x axis track ------------------------------------------------
## This is so that the key positions can be matched with the heatmap
p <- p + scale_x_continuous(expand = c(0,0), breaks = x_axis_brks/1000000, position = "top",
                       sec.axis =  dup_axis(name = NULL,
                                            breaks = x_axis_brks2/1000000, labels = NULL))


# library(gtable)
# library(grid)
# g <- ggplotGrob(p)
# strips <- g$layout[grep("strip-l", g$layout$name), ]
# titles <- lapply(paste0("(", letters[seq_len(nrow(strips))], ")"), 
#                  textGrob, x = 0, hjust = 0, vjust = 1)
# g <- gtable_add_grob(g, grobs = titles, 
#                      t = strips$t - 1, b = strips$b - 2, 
#                      l = strips$l, r = strips$r)
# grid.newpage()
# grid.draw(g)
# p2 <- ggplot() +
#   geom_rect(dt_genes, mappin = aes(xmin=start_pos, xmax=end_pos,
#                                  ymin= 10, ymax= 20), col = "blue", fill = NA) +
#   scale_x_continuous(limits = c(start_pos, end_pos), expand = c(0,0)) + 
#   scale_y_continuous(expand = c(0,0)) +
#   theme_classic() +
#   theme(axis.line = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())


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
        axis.text.x = element_blank()) + 
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# phm
phm <- phm + geom_vline(xintercept = val,
                 linetype = "dashed", colour = "black", size = 0.5)

# phm 
# allset2plots <- p + phm + plot_layout(ncol = 1, nrow = 2, heights = c(5, 3))

allset2plots <- cowplot::plot_grid(p + theme(legend.position = "top", legend.justification = "center"), 
                                   phm #+ geom_text(x = diff(range(dt_AF$Rank)) * labrat, y = 11*0.99, label = "E)", fontface = "bold", size = 15/.pt)
                                   , 
                                   align = "v", 
                                   ncol = 1, 
                                   # axis = "lr", 
                                   rel_heights = c(13, 9))

# ggarrange(p, phm, ncol = 1, align = "v", heights = c(5,4))

png(outpng, width = 8.3, height = 10, units = "in", res = 600)
allset2plots
dev.off()


# png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/combined_chr8_XtX_C2_IHS_DP_AFV2.png", width = 8.3, height = 11.7, units = "in", res = 600)
# cowplot::plot_grid(p + theme(legend.position = "top", legend.justification = "center"), 
#                    phm + scale_x_continuous(expand = c(0,0), limits = c(1220, 2100))
#                    , 
#                    align = "v", 
#                    ncol = 1, 
#                    # axis = "lr", 
#                    rel_heights = c(5, 4), labels = c("A)", "B)"))
# dev.off()

  




