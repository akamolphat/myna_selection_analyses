library(vcfR)
library(data.table)
library(tidyverse)
inVCF <- read.vcfR("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.snps.vcf.gz")
gt <- extract.gt(inVCF, element = "GT", )
gt_pos <- gsub('^.*Superscaffold_chr8_*|_.*$', '', row.names(gt))

# POSITION of interest
# XtX
# Superscaffold_chr8_20694395_A_G
# Superscaffold_chr8_20673781_A_G
# C2
# Superscaffold_chr8_20694395_A_G # Same as XtX
# Superscaffold_chr8_20677064_A_C # Same as missense mutation on copy 2 AMY2A

ls_POS <- c(20677064, 20694395, 20673781)
gt_mat <- t(gt[gt_pos %in% ls_POS,])
gt_mat[gt_mat == "1/1" & !is.na(gt_mat)] <- 2
gt_mat[gt_mat == "0/1" & !is.na(gt_mat)] <- 1
gt_mat[gt_mat == "0/0" & !is.na(gt_mat)] <- 0
dtx <- data.table(gt_mat)
dtx$ID <- row.names(gt_mat)

dt_meta <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/metadata/TableS1.2.csv")
dt_meta$SV_genotype[dt_meta$SV_genotype == "1-Jan" & !is.na(dt_meta$SV_genotype)] <- 2
dt_meta$SV_genotype[dt_meta$SV_genotype == "0/1" & !is.na(dt_meta$SV_genotype)] <- 1
dt_meta$SV_genotype[dt_meta$SV_genotype == "0/0" & !is.na(dt_meta$SV_genotype)] <- 0

dt_merged <- merge(dt_meta, dtx, by.x = "Sample_ID", by.y = "ID", all.x = T)
# write_csv(dt_merged, file = "data/metadata/TableS1.2_complete.csv")


dt_merged <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/metadata/TableS1.2_complete.csv")

dt_merged_sub <- dt_merged %>%
  filter(Removed == "No") %>%
  mutate(SV_deletion_no = 2 - SV_insertion_no) %>%
  select(Sample_ID, pop3_abb, Superscaffold_chr8_20673781_A_G,
         Superscaffold_chr8_20677064_A_C, SV_deletion_no, Superscaffold_chr8_20694395_A_G) %>%
  pivot_longer(
    cols = c("Superscaffold_chr8_20673781_A_G",
             "Superscaffold_chr8_20677064_A_C", "SV_deletion_no", "Superscaffold_chr8_20694395_A_G"),
    names_to = "loci",
    values_to = "value"
  ) 

dt_merged_sub$POS <- as.numeric(sapply(strsplit(x = dt_merged_sub$loci, split = "_"), "[[" , 3))
dt_merged_sub$POS[dt_merged_sub$loci == "SV_deletion_no"] <- 20680000

dt_merged_sub <- dt_merged_sub %>% 
  mutate(Rank = dense_rank(POS)) %>%
  mutate(label = case_when(Rank == 3 ~ "20678727-20687524 ",
                           TRUE ~ as.character(POS))) #as.character is required to match column class

dt_merged_sub$value <- as.factor(dt_merged_sub$value)
# dt_merged_sub$pop <- gsub(pattern = "_", replacement = " ", x = dt_merged_sub$pop3_abb, perl = T)
# dt_merged_sub$pop[dt_merged_sub$pop == "Maharashtra subpopulation A"] <- "Maharashtra\nsubpop. A"
# dt_merged_sub$pop[dt_merged_sub$pop == "Madhya Pradesh"] <- "Madhya\nPradesh"
# dt_merged_sub$pop[dt_merged_sub$pop == "Great Barrier Island"] <- "Great Barrier\nIsland"

poplab <- c("Madhya\nPradesh", "Tamil Nadu", "Maharashtra", "Maharashtra\nsubpop. A", "Fiji", "Melbourne", "Cairns", "Napier", "Leigh", "Great\nBarrier Isl.", "South\nAfrica")
brk <- c("Madhya_Pradesh", "Tamil_Nadu", "Maharashtra", "Maharashtra_subpopulation_A", "Fiji", "Melbourne", "Cairns", "Napier", "Leigh", "Great_Barrier_Island", "South_Africa")
dt_poplab <- data.table(pop3_abb = brk, pop = factor(poplab, levels = rev(poplab)))
dt_merged_sub <- merge(dt_merged_sub, dt_poplab, by = "pop3_abb")
dt_pos <- unique(dt_merged_sub[,c("Rank", "label")])
## Genotype heatmap -------------------------------------------------------
gthm <- ggplot(dt_merged_sub, aes(x = Rank, y = Sample_ID, fill = value)) +
  geom_tile() + 
  labs(
    x = "Position",
    y = NULL
  ) +
  scale_x_continuous(breaks = dt_pos$Rank, labels = dt_pos$label, expand = c(0,0), 
                     sec.axis = sec_axis(~., 
                                         breaks = c(1,2,3,4), 
                                         labels = c("XtX peak", 
                                                    "C2 peak &\nAMY2A copy two\nmissense mutation","SV",
                                                    "XtX & C2 peak"))) +
  facet_grid(rows = vars(pop), switch = "y", scales = "free_y", space = "free_y") +
  scale_fill_manual(name = "Genotype", 
                    breaks = c(0,1,2, NA),
                    values = c("#fc8d59", "#ffffbf", "#91bfdb", "lightgrey"),
                    # values = c("#fee8c8", "#fdbb84", "#e34a33", "lightgrey"),
                    labels = c("0/0 (Insertion)", "0/1", "1/1 (Deletion)", "Not genotyped")) + 
  scale_y_discrete(expand = c(0,0)) + 
  geom_vline(xintercept = c(1.5, 2.5, 3.5), 
             colour = "black", linewidth = 0.1)+
  # ggtitle("Allele frequency vs population (standardised to Madhya Pradesh)") + 
  theme(strip.placement = "outside",
        panel.spacing = unit(0.1, "lines"),
        legend.background = element_rect(colour = "black"),
        legend.position = "top",
        # legend.direction = "vertical",
        legend.margin = margin(5, 12, 5, 12),
        legend.title = element_text(face = "bold"),
        axis.text.x.top = element_text(colour = "black", face = "bold", size = 10),
        strip.text.y.left = element_text(colour = "black", face = "bold", size = 10, angle = 0),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(colour = "black", face = "bold", size = 7),
        axis.title.x = element_text(colour = "black", face = "bold"),
        plot.margin = unit(c(9, 7, 7, 7), "points")) 

gthm
png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/chr8_GT_key_outlier_XtX_C2_SV_SNPs_top.png", width = 8.3, height = 10.5, units = "in", res = 600)
gthm 
dev.off()

pdf("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/chr8_GT_key_outlier_XtX_C2_SV_SNPs_top.pdf", width = 8.3, height = 10.5)
gthm 
dev.off()

