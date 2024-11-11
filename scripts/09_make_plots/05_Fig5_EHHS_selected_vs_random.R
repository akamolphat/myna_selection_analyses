# Calculate EHHS for SNP of interest --------------------------------------
## Library ------
library(rehh)
library(data.table)
## Inputs -------
popmapfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/plink_QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10/pop_map.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.txt"
ls_mrk_id <- c("Superscaffold_chr8_20673781_A_G", "Superscaffold_chr8_20677064_A_C", "Superscaffold_chr8_20694395_A_G", # Selected SNPs
               "Superscaffold_chr8_21800137_T_C", "Superscaffold_chr8_19599941_G_A") # Random SNPs

# POSITION of interest
# XtX
# Superscaffold_chr8_20694395_A_G
# Superscaffold_chr8_20673781_A_G
# C2
# Superscaffold_chr8_20694395_A_G # Same as XtX
# Superscaffold_chr8_20677064_A_C # Same as missense mutation on copy 2 AMY2A

EHHS_from_hapfile <- function(hapfile, popmapfile, ls_mrk_id){
  # This function calculates EHHS for the listed markers (ls_mrk_id)
  # in the differeen populations as listed in the popmapfile
  # using the haplotype data from the hapfile
  require(rehh)
  require(data.table)
  hh <- data2haplohh(hap_file = hapfile,
                     polarize_vcf = FALSE,
                     vcf_reader = "data.table")
  
  # Read in pop_ind to subset hh
  dt_pop_map <- fread(popmapfile,
                      header = F,
                      col.names = c("ID", "pop"))
  
  extract_format_id <- function(dt_pop_map, hh, pop){
    ID_MP <- dt_pop_map$ID[dt_pop_map$pop %in% pop]
    ID_MP1 <- paste(ID_MP, "1", sep = "_")
    ID_MP2 <- paste(ID_MP, "2", sep = "_")
    ID_MP_subset <- rownames(hh@haplo)[rownames(hh@haplo) %in% c(ID_MP1, ID_MP2)]
    return(ID_MP_subset)
  }
  
  EHHS_pop <- function(dt_pop_map, hh, pop, MRK_ID, poplab){
    ID_MPR_subset <- extract_format_id(dt_pop_map = dt_pop_map, hh = hh, pop = pop)
    hh_MPR <- subset(hh, select.hap = ID_MPR_subset)
    res_MPR <- calc_ehhs(hh_MPR, mrk = MRK_ID, include_nhaplo = TRUE)
    dt1 <- res_MPR$ehhs
    dt1$pop <- poplab
    return(dt1)
  }
  
  # ls_mrk_id <- c("Superscaffold_chr8_20672140_A_G", "Superscaffold_chr8_20677064_A_C", # Selected SNPs
  #                "Superscaffold_chr8_19599941_G_A", "Superscaffold_chr8_21800137_T_C") # Random SNPs
  # 
  # ls_mrk_id <- c("Superscaffold_chr8_20692138_A_C", "Superscaffold_chr8_20672140_A_G", # Selected SNPs
  #                "Superscaffold_chr8_19599941_G_A", "Superscaffold_chr8_21800137_T_C") # Random SNPs
  
  ctr <- 1
  for (mrkid in ls_mrk_id){
    for (i in unique(dt_pop_map$pop)){
      print(i)
      if (ctr == 1){
        dt_comb <- EHHS_pop(dt_pop_map = dt_pop_map, hh = hh, pop = i, MRK_ID = mrkid, poplab = gsub(pattern = "_", replacement = " ", x = i, fixed = T))
        dt_comb$MRK_ID = mrkid
      } else {
        dt_x <- EHHS_pop(dt_pop_map = dt_pop_map, hh = hh, pop = i, MRK_ID = mrkid, poplab = gsub(pattern = "_", replacement = " ", x = i, fixed = T))
        dt_x$MRK_ID = mrkid
        dt_comb <- rbind(dt_comb, dt_x)
      }
      ctr <- ctr + 1
    }
    
  }
  return(dt_comb)
}

cleanpopnames <- function(dt_comb){
  # This function clean population names and assign native and invasive populations
  # and is very specific to this particular table
  dt_comb$pop[dt_comb$pop == "Maharashtra subpopulation A"] <- "Maharashtra subpop. A"
  dt_comb$range <- "Invasive"
  dt_comb$range[dt_comb$pop %in% c("Madhya Pradesh", "Maharashtra", "Maharashtra subpop. A", "Tamil Nadu")] <- "Native"
  dt_comb$outlier <- "Random"
  dt_comb$outlier[dt_comb$MRK_ID %in% c("Superscaffold_chr8_20673781_A_G", 
                                        "Superscaffold_chr8_20694395_A_G", "Superscaffold_chr8_20677064_A_C")] <- "Outlier"
  dt_comb$POS_ID <- paste("Position:", gsub(pattern = ".*Superscaffold_chr8_(.+)_.*_.*", "\\1", x = dt_comb$MRK_ID))
  dt_comb$MRK_POS <- as.numeric(gsub(pattern = "Position: ", replacement = "", dt_comb$POS_ID))
  dt_comb$POSITION2 <- dt_comb$POSITION - dt_comb$MRK_POS
  return(dt_comb)
}


## Original dataset --------------------------------------------------------
hapfile1 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.phased_shapeit5.snps.vcf.gz"
outEHHS1 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/EHHS_selected_vs_random_withSV.txt"
dt_comb <- EHHS_from_hapfile(hapfile = hapfile1,
                             popmapfile = popmapfile,
                             ls_mrk_id = ls_mrk_id)
### Clean up populations names  -----
dt_comb <- cleanpopnames(dt_comb)

### Output to txt -----
fwrite(dt_comb, file = outEHHS1, sep = "\t", quote = F)

## SVremoved dataset -------------------------------------------------------
hapfile2 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.SVsnpsremoved.phased_shapeit5.snps.vcf.gz"
outEHHS2 <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/EHHS_selected_vs_random_withoutSV.txt"
dt_comb2 <- EHHS_from_hapfile(hapfile = hapfile2,
                              popmapfile = popmapfile,
                              ls_mrk_id = ls_mrk_id)
### Clean up populations names  -----
dt_comb2 <- cleanpopnames(dt_comb2)

### Output to txt -----
fwrite(dt_comb2, file = outEHHS2, sep = "\t", quote = F)

# Making EHHS plots -------------------------------------------------------
# dt_comb <- fread(file = outEHHS1)
# dt_comb2 <- fread(file = outEHHS2)
ls_mrk_id <- c("Superscaffold_chr8_19599941_G_A", "Superscaffold_chr8_20677064_A_C", "Superscaffold_chr8_20694395_A_G", "Superscaffold_chr8_20673781_A_G",  # Selected SNPs
               "Superscaffold_chr8_21800137_T_C") #) # Random SNPs
# 4 panel plot of the following markers:
ls_mrk_id2 <- c("Superscaffold_chr8_20677064_A_C", "Superscaffold_chr8_20694395_A_G", "Superscaffold_chr8_20673781_A_G",  # Selected SNPs
                "Superscaffold_chr8_19599941_G_A") # "Superscaffold_chr8_21800137_T_C", Random SNPs

# 2 panel plot of the following markers:
ls_mrk_id3 <- c("Superscaffold_chr8_20694395_A_G", "Superscaffold_chr8_19599941_G_A")

CLSTorder <- c("Madhya Pradesh", "Tamil Nadu", "Maharashtra", "Maharashtra subpop. A", "Fiji", "Melbourne", "Cairns", "Napier", "Leigh", "Great Barrier Island", "South Africa")

## Original dataset --------------------------------------------------------
dt_comb <- fread(outEHHS1) %>%
  filter(MRK_ID %in% ls_mrk_id)
### Clean up -----
dt_comb$pop <- factor(dt_comb$pop, levels = CLSTorder)
override.linetype <- c(1,1,1,1,3,3,3,3,3,3,3)
dt_comb_sum <- dt_comb %>% 
  group_by(POS_ID) %>%
  summarise(diff_centre = max(c(abs(POSITION - MRK_POS))),
            MRK_POS = unique(MRK_POS))
dt_dummy <- unique(dt_comb[,c("outlier", "POS_ID")])
dt_fake_datapoints <- merge(dt_comb_sum, dt_dummy)
dt_fake_datapoints$EHHS <- 0
dt_fake_datapoints2 <- dt_fake_datapoints
dt_fake_datapoints$POSITION <- dt_fake_datapoints$MRK_POS + max(dt_comb_sum$diff_centre)
dt_fake_datapoints2$POSITION <- dt_fake_datapoints2$MRK_POS - max(dt_comb_sum$diff_centre)
dt_fake <- rbind(dt_fake_datapoints2, dt_fake_datapoints)
plot_panel <- function(dt_comb, dt_fake, override.linetype){
  p <- ggplot(data = dt_comb) +
    geom_line(mapping = aes(x = POSITION/1000000, y = EHHS, col = pop, lty = range)) + 
    geom_point(data = dt_fake, mapping = aes(x = POSITION/1000000, y = EHHS), col = NA) +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(data = unique(dt_comb[,c("outlier", "POS_ID", "MRK_POS")]),
               aes(xintercept = MRK_POS/1000000), lty = "dashed") +
    scale_linetype_manual(values = c(1,3), breaks = c("Native", "Invasive")) +
    guides(col = guide_legend(title = "Population", 
                              title.position="top", title.hjust = 0.5,
                              override.aes = list(linetype = override.linetype)), 
           lty = guide_legend(title = "Range", 
                              title.position="top", title.hjust = 0.5)) +
    xlab("Position (Mb)") +
    facet_wrap(outlier~POS_ID, scales = "free_x", ncol = 2) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.background = element_rect(colour = "black")) 
  return(p)
}
### Plot 5 panel plot ----
p0a <- plot_panel(dt_comb, dt_fake, override.linetype) +  
  facet_wrap(POS_ID~outlier, scales = "free_x", ncol = 2) 

png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/EHHS_chr8_selected_vs_random_5panel_withSV_revision1.png", 
    width = 8.3, height = 8.3, units = "in", res = 300)
p0a
dev.off()

### Zoomed in 5 panel plot ----
dt_fake_datapoints$POSITION <- dt_fake_datapoints$MRK_POS + 50000
dt_fake_datapoints2$POSITION <- dt_fake_datapoints2$MRK_POS - 50000
dt_fake2 <- rbind(dt_fake_datapoints2, dt_fake_datapoints)

dt_comb_50kb <- dt_comb %>% filter(POSITION2 <= 50000 & POSITION2 >= -50000)

p0b <- plot_panel(dt_comb_50kb, dt_fake2, override.linetype) +  
  facet_wrap(POS_ID~outlier, scales = "free_x", ncol = 2) 
png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/EHHS_chr8_selected_vs_random_5panel_zoomed_withSV_revision1.png", 
    width = 8.3, height = 8.3, units = "in", res = 300)
p0b
dev.off()

### Plot 4 panel plot ----
dt_comb_4panel <- dt_comb %>% 
  filter(MRK_ID %in% ls_mrk_id2)
dt_comb_4panel$POS_ID <- factor(dt_comb_4panel$POS_ID, levels = unique(dt_comb_4panel$POS_ID))
dt_fake_4panel <- dt_fake2 %>%
  filter(POS_ID %in% unique(dt_comb_4panel$POS_ID))
p1a <- plot_panel(dt_comb_4panel, dt_fake_4panel, override.linetype)
png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/EHHS_chr8_selected_vs_random_4panel_withSV_revision1.png", 
    width = 8.3, height = 8.3, units = "in", res = 300)
p1a
dev.off()
### Zoomed in 4 panel plot -----
dt_comb_50kb_4panel <- dt_comb_50kb %>% 
  filter(MRK_ID %in% ls_mrk_id2)
dt_fake2_sub_4panel <- dt_fake2 %>%
  filter(POS_ID %in% unique(dt_comb_50kb_4panel$POS_ID))

p1b <- plot_panel(dt_comb_50kb_4panel, dt_fake2_sub_4panel, override.linetype)
png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/EHHS_chr8_selected_vs_random_4panel_zoomed_withSV_revision1.png", 
    width = 8.3, height = 8.3, units = "in", res = 300)
p1b
dev.off()
### 2 panel plot -----
dt_comb_sub <- dt_comb %>% 
  filter(MRK_ID %in% ls_mrk_id3)
dt_fake_sub <- dt_fake2 %>%
  filter(POS_ID %in% unique(dt_comb_sub$POS_ID))

p2a <- plot_panel(dt_comb_sub, dt_fake_sub, override.linetype)
png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/EHHS_chr8_selected_vs_random_2panel_withSV_revision1.png", 
    width = 8.3, height = 5, units = "in", res = 300)
p2a
dev.off()

### Zoomed in 2 panel plot ------
dt_comb_50kb_sub <- dt_comb_50kb %>% 
  filter(MRK_ID %in% ls_mrk_id3)
dt_fake2_sub <- dt_fake2 %>%
  filter(POS_ID %in% unique(dt_comb_50kb_sub$POS_ID))

p2b <- plot_panel(dt_comb_50kb_sub, dt_fake2_sub, override.linetype)
png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/EHHS_chr8_selected_vs_random_2panel_zoomed_withSV_revision1.png", 
    width = 8.3, height = 5, units = "in", res = 300)
p2b
dev.off()

pdf("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/EHHS_chr8_selected_vs_random_2panel_zoomed_withSV_revision1.pdf", 
    width = 8.3, height = 5)
p2b
dev.off()
