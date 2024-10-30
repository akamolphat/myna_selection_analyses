library(rehh)
library(data.table)
library(tidyverse)

hapfile <- "data/processed/BCFtools/vcftools_filtered/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.SVsnpsremoved.phased_shapeit5.snps.vcf.gz"
hh <- data2haplohh(hap_file = hapfile,
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

mrk.nr.start = which(mrk.names(hh) == "Superscaffold_chr8_20600034_C_T")

mrk.nr.end = which(mrk.names(hh) == "Superscaffold_chr8_20749985_C_T")

# Peak 1: Superscaffold_chr8_20669616_C_T
# Last SNP before SV: "Superscaffold_chr8_20678534_T_G"
# First SNP after SV: Superscaffold_chr8_20688007_A_T
# Peak 2: Superscaffold_chr8_20691861_A_G 
hh_subset <- subset(hh, select.mrk = mrk.nr.start:mrk.nr.end)
haplo <- hh_subset@haplo
popmap <- "data/processed/BCFtools/vcftools_filtered/pop_map_WGS_ALL.txt"
hh2longdt <- function(haplo, popmap){
  indnames <- rownames(haplo)
  snpnames <- colnames(haplo) 
  dt1 <- as.data.table(haplo)
  dt1$haplo_name <- indnames
  # dt1$POS <- gsub(x = gsub(x = snpnames, pattern = "Superscaffold_chr8_", replacement = ""), pattern = "_[A-Z]_[A-Z]", replacement = "")
  dt1 <- dt1 %>% pivot_longer(!haplo_name, names_to = "ID", values_to = "presence_absence")
  dt1$POS <- gsub(x = gsub(x = dt1$ID, pattern = "Superscaffold_chr8_", replacement = ""), pattern = "_[A-Z]_[A-Z]", replacement = "")
  dt_popmap <- fread(popmap, header = F, col.names = c("individual", "pop"))
  dt_trans <- data.table(haplo_name = c(paste(dt_popmap$individual, 1, sep = "_"),
                                        paste(dt_popmap$individual, 2, sep = "_")),
                         individual = rep(dt_popmap$individual, times = 2))
  dt_popmap <- merge(dt_trans, dt_popmap, by = "individual", all = T)
  dt_popmap <- dt_popmap %>%
    filter(haplo_name %in% unique(dt1$haplo_name))
  dt1a <- merge(dt_popmap, dt1, by = "haplo_name", all = T)
  return(dt1a)
}

dt1a <- hh2longdt(haplo, popmap = popmap)
dt1a$presence_absence <- as.factor(dt1a$presence_absence)
dt1a$presence_absence[is.na(dt1a$presence_absence)] <- "missing"
# dt1a$presence <- dt1a$presence_absence
# dt1a$presence[dt1a$presence_absence != 1] <- NA
dt1a$POS <- as.numeric(dt1a$POS)

dt1a$pop <- gsub(pattern = "_", replacement = " ", x = dt1a$pop, perl = T)
dt1a$pop[dt1a$pop == "Maharashtra subpopulation A"] <- "Maharashtra\nsubpop. A"
poporder <- rev(c("Madhya Pradesh", "Tamil Nadu", "Maharashtra", "Maharashtra\nsubpop. A", "Fiji", "Melbourne", "Cairns", "Napier", "Leigh", "Great Barrier Island", "South Africa"))
dt1a$pop <- factor(x = dt1a$pop,
                   levels = poporder)

dt_RE <- data.table(startpos=rep(20678727, times = 11),
                    endpos=rep(20687524, times = 11),
                    unique(dt1a$pop))
haplo_hm <- ggplot() +
  geom_line(data = dt1a, mapping = aes(x = POS/1000000, y = haplo_name, colour = individual), linewidth = 1.5) +
  scale_colour_manual(values = rep(c("lightblue", "khaki1"), times = 40)) +
  ## Mark key outlier SNP
  # Peak outlier SNP before and afer SV
  # XtX
  # Superscaffold_chr8_20694395_A_G
  # Superscaffold_chr8_20673781_A_G
  # C2
  # Superscaffold_chr8_20694395_A_G # Same as XtX
  # Superscaffold_chr8_20677064_A_C # Same as missense mutation on copy 2 AMY2A
  geom_vline(xintercept = 20673781/1000000, lty = "dashed", linewidth = 0.2) +
  geom_vline(xintercept = 20677064/1000000, lty = "dashed", linewidth = 0.2) +
  geom_vline(xintercept = 20694395/1000000, lty = "dashed", linewidth = 0.2) +
  ##
  geom_rect(dt_RE, mapping = aes(xmin=startpos/1000000, xmax=endpos/1000000,
                                 ymin= -Inf, ymax= Inf),
            alpha = 0.2,
            fill = "lightgrey", colour = NA) +
  geom_point(data = dt1a %>% filter(presence_absence == 1), mapping = aes(x = POS/1000000, y = haplo_name), colour = "red", size = 0.2) +
  # facet_wrap(~pop, 
  #            scales = "free_y",
  #            ncol = 1,
  #            strip.position = "right") +
  facet_grid(pop~., switch = "y", scales="free_y", space = "free_y") +
  scale_x_continuous(expand = c(0,0)) +
  # scale_y_discrete(limits=rev) +
  xlab("Position (Mb)") +
  ylab("Individual haplotype") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 5),
        legend.position = "none",
        strip.placement = "outside",
        strip.text.y = element_text(size = 7)) 


png("results/EHHS/chr8_SV_hap_struc.png", width = 8.3, height = 11.7, units = "in", res = 600)
haplo_hm
dev.off()
