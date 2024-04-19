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
  geom_vline(xintercept = 20669616/1000000, lty = "dashed", linewidth = 0.2) +
  geom_vline(xintercept = 20691861/1000000, lty = "dashed", linewidth = 0.2) +
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

# Plot iHH over the region
infile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/chr8_SVsnpsremoved_scanhh.txt"
scan <- fread(infile) %>% filter(CHR == "Superscaffold_chr8" & POSITION >= 20600000 & POSITION <= 20750000)
infile2 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/chr8_scanhh.txt"
scan2 <- fread(infile2) %>% filter(CHR == "Superscaffold_chr8" & POSITION >= 20600000 & POSITION <= 20750000)

png("results/EHHS/chr8_SV_iHH.png", width = 8.3, height = 8.3, units = "in", res = 600)
par(mfrow = c(2,1))
# Plot iHH over the region when SV is not included
plot(x = scan$POSITION/1000000, y = scan$IHH_A, xlab = "Position (Mb)", ylab = "iHH", col = NA, ylim = c(0, max(c(scan$IHH_A, scan$IHH_D))))
rect(xleft = 20.678727, xright = 20.687524, ybottom = 0, ytop = max(c(scan$IHH_A, scan$IHH_D)), col = "lightgrey", border = NA)
points(x = scan$POSITION/1000000, y = scan$IHH_A, col = "red")
points(x = scan$POSITION/1000000, y = scan$IHH_D, col = "blue", pch = 2)
abline(v = 20669616/1000000, lty = "dashed")
abline(v = 20691861/1000000, lty = "dashed")
title("iHH over chr8 outlier peak (SNPs from SV excluded)")
legend(x = 20.6, y = 35000, col = c("red", "blue"), pch = c(1,2), legend = c(expression("iHH"[Major]), expression("iHH"[Minor])))
legend(x = 20.6, y = 48000, lty = c("dashed"), legend = c("iHS outlier peak"))
# Plot iHH over the region when SNPs in SV were included
plot(x = scan2$POSITION/1000000, y = scan2$IHH_A, xlab = "Position (Mb)", ylab = "iHH", col = NA, ylim = c(0, max(c(scan2$IHH_A, scan2$IHH_D))))
rect(xleft = 20.678727, xright = 20.687524, ybottom = 0, ytop = max(c(scan2$IHH_A, scan2$IHH_D)), col = "lightgrey", border = NA)
points(x = scan2$POSITION/1000000, y = scan2$IHH_A, col = "red")
points(x = scan2$POSITION/1000000, y = scan2$IHH_D, col = "blue", pch = 2)
abline(v = 20667691/1000000, lty = "dashed")
abline(v = 20692138/1000000, lty = "dashed")
title("iHH over chr8 outlier peak (SNPs from SV included)")

dev.off()

# Plot iHS SV removed -----------------------------------------------------
dtx <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/WGS_IHS_chr8_SVremoved.txt") %>%
  filter(CHR == "Superscaffold_chr8" & POSITION >= 20600000 & POSITION <= 20750000)

# Plot iHH over the region ------------------------------------------------
infile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/chr8_SVsnpsremoved_scanhh.txt"
scan <- fread(infile) %>% filter(CHR == "Superscaffold_chr8" & POSITION >= 20600000 & POSITION <= 20750000)
png("results/EHHS/chr8_SV_iHS_iHH_snpsremoved.png", width = 8.3, height = 8.3, units = "in", res = 600)
par(mfrow = c(2,1))
plot(x = dtx$POSITION/1000000, y = dtx$LOGPVALUE, xlab = "Position (Mb)", ylab = "iHS")
rect(xleft = 20.678727, xright = 20.687524, ybottom = 0, ytop = max(c(scan$IHH_A, scan$IHH_D)), col = "lightgrey", border = NA)
abline(v = 20669616/1000000, lty = "dashed")
abline(v = 20691861/1000000, lty = "dashed")
title("iHS over chr8 outlier peak (SNPs from SV excluded)")

# Plot iHH over the region when SV is not included
plot(x = scan$POSITION/1000000, y = scan$IHH_A, xlab = "Position (Mb)", ylab = "iHH", col = NA, ylim = c(0, max(c(scan$IHH_A, scan$IHH_D))))
rect(xleft = 20.678727, xright = 20.687524, ybottom = 0, ytop = max(c(scan$IHH_A, scan$IHH_D)), col = "lightgrey", border = NA)
points(x = scan$POSITION/1000000, y = scan$IHH_A, col = "red")
points(x = scan$POSITION/1000000, y = scan$IHH_D, col = "blue", pch = 2)
abline(v = 20669616/1000000, lty = "dashed")
abline(v = 20691861/1000000, lty = "dashed")
title("iHH over chr8 outlier peak (SNPs from SV excluded)")
legend(x = 20.6, y = 35000, col = c("red", "blue"), pch = c(1,2), legend = c(expression("iHH"[Major]), expression("iHH"[Minor])))
legend(x = 20.6, y = 48000, lty = c("dashed"), legend = c("iHS outlier peak"))

dev.off()

# Plot iHS SV removed pos shifted -----------------------------------------
dtx <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/WGS_IHS_chr8_SVremoved_shifted.txt") %>%
  filter(CHR == "Superscaffold_chr8" & POSITION >= 20600000 & POSITION <= 20750000)

# Plot iHH over the region ------------------------------------------------
infile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/chr8_SVsnpsremoved_pos_shifted_scanhh.txt"
scan <- fread(infile) %>% filter(CHR == "Superscaffold_chr8" & POSITION >= 20600000 & POSITION <= 20750000)
png("results/EHHS/chr8_SV_iHS_iHH_pos_shifted.png", width = 8.3, height = 8.3, units = "in", res = 600)
par(mfrow = c(2,1))
plot(x = dtx$POSITION/1000000, y = dtx$LOGPVALUE, xlab = "Position (Mb)", ylab = "iHS")
abline(v = 20678727/1000000, lty = "dashed")
title("iHS over chr8 outlier peak (SNPs from SV excluded and position shifted)")

# Plot iHH over the region when SV is not included
plot(x = scan$POSITION/1000000, y = scan$IHH_A, xlab = "Position (Mb)", ylab = "iHH", col = NA, ylim = c(0, max(c(scan$IHH_A, scan$IHH_D))))
# rect(xleft = 20.678727, xright = 20.687524, ybottom = 0, ytop = max(c(scan$IHH_A, scan$IHH_D)), col = "lightgrey", border = NA)
points(x = scan$POSITION/1000000, y = scan$IHH_A, col = "red")
points(x = scan$POSITION/1000000, y = scan$IHH_D, col = "blue", pch = 2)
abline(v = 20678727/1000000, lty = "dashed")
title("iHH over chr8 outlier peak (SNPs from SV excluded and position shifted)")
legend(x = 20.6, y = 100000, col = c("red", "blue"), pch = c(1,2), legend = c(expression("iHH"[Major]), expression("iHH"[Minor])))
legend(x = 20.63, y = 100000, lty = c("dashed"), legend = c("SV position"))

dev.off()