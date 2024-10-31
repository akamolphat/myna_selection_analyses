library(data.table)
library(tidyverse)
dtcomb <- fread(file = "data/processed/baypass/mac10/combined/outliers/myna_baypass_mac10_combined_xtx_outliers_0.99999_cluster1_summary.txt")
dtcomb$outlier_statistics <- "XtX"
bp_region_dist <- 50000 # Cluster within 50kbp of each other
# con_id <- 1:7
con_id <- c(1,4)  # Only CON001 and CON004 were needed and used in this study. 
# CON001 = Native vs invasive
# CON004 = Native vs invasive populations founded by one-step introductions (i.e., Melbourne, Fiji, and South Africa) 
for (i in con_id){
  infile <- paste("data/processed/baypass/mac10/combined/outliers/myna_baypass_mac10_combined_IS_C2_GEA_summary_contrast_CON_00", 
                  i, 
                  "_outliers_0.99999_cluster_summary.txt", sep = "")
  dtcon <- fread(file = infile)
  dtcon$outlier_statistics <- paste("CON00", i, sep = "")
  dtcomb <- rbind(dtcomb, dtcon)
}

dtcomb5 <- dtcomb %>%
  filter(n > 5)

dtcomb10 <- dtcomb %>%
  filter(n > 10)

group_outliers_cluster <- function(dt){
  dt <- dt[with(dt, order(chr, start_pos, end_pos)), ]
  dt$overlap <- NA
  ctr <- 1
  maxposbool <- F
  for(i in 1:(nrow(dt) - 1) ){
    if (dt$chr[i] == dt$chr[i + 1]){
      if (maxposbool){
        maxpos
      } else {
        maxpos <- dt$end_pos[i]
      }
      startpos <- dt$start_pos[i + 1]
      if (maxpos >= startpos){
        if (maxpos > dt$end_pos[i + 1]) {
          maxposbool <- T
        } else {
          maxposbool <- F
        }
        dt$overlap[i] <- ctr
        dt$overlap[i + 1] <- ctr
      } else {
        # dt$overlap[i] <- ctr
        ctr <- ctr + 1
      }
    } else {
      # dt$overlap[i] <- ctr
      ctr <- ctr + 1
      maxposbool <- F
    }
  }
  
  dt <- dt[with(dt, order(outlier_statistics, chr, start_pos, end_pos)), ]
  return(dt)
}

dtcomb10 <- group_outliers_cluster(dt = dtcomb10)
dtcomb5 <- group_outliers_cluster(dt = dtcomb5)

dir.create(path = "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/combined_outliers/")
write_csv(x = dtcomb10, file = "results/combined_outliers/outliers_n10_summary.csv")
write_csv(x = dtcomb5, file = "results/combined_outliers/outliers_n5_summary.csv")

dtbed10 <- dtcomb10 |>
  dplyr::mutate(start_pos = start_pos - bp_region_dist,
                end_pos = end_pos + bp_region_dist) 
dtbed10$start_pos[dtbed10$start_pos < 0] <- 0

fwrite(x = dtbed10[,c("chr", "start_pos", "end_pos", "outlier_statistics", "cluster1")], ## One extra column for ease of use
       file = "results/combined_outliers/outliers_n10_summary.bed", col.names = F, sep = "\t") 

dtbed5 <- dtcomb5 |>
  dplyr::mutate(start_pos = start_pos - bp_region_dist,
                end_pos = end_pos + bp_region_dist) 
dtbed5$start_pos[dtbed5$start_pos < 0] <- 0

fwrite(x = dtbed5[,c("chr", "start_pos", "end_pos", "outlier_statistics", "cluster1")], ## One extra column for ease of use
       file = "results/combined_outliers/outliers_n5_summary.bed", col.names = F, sep = "\t") 


