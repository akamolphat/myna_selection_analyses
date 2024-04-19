library(data.table)
library(tidyverse)
dt <- fread(file = "results/combined_outliers/outliers_n10_summary_genelist.txt",
            header = F,
            col.names = c("chr", "start_pos", "end_pos", "outlier_statistics", "cluster1", "gene_start_pos", "gene_end_pos", "gene_name"))

dtsum <- dt %>% 
  group_by(chr, start_pos, end_pos, cluster1) %>%
  summarise(genes = paste(gene_name, collapse = ", ")) 

dtmeta <- read_csv(file = "results/combined_outliers/outliers_n10_summary.csv")

dtmeta <- merge(dtmeta, dtsum[,c("chr", "cluster1", "genes")], by = c("chr", "cluster1"), all.x = T)

write_csv(dtmeta, file = "results/combined_outliers/outliers_n10_summary_genenames.csv")


dt <- fread(file = "results/combined_outliers/outliers_n5_summary_genelist.txt",
            header = F,
            col.names = c("chr", "start_pos", "end_pos", "outlier_statistics", "cluster1", "gene_start_pos", "gene_end_pos", "gene_name"))

dtsum <- dt %>% 
  group_by(chr, start_pos, end_pos, cluster1) %>%
  summarise(genes = paste(gene_name, collapse = ", ")) 

dtmeta <- read_csv(file = "results/combined_outliers/outliers_n5_summary.csv")

dtmeta <- merge(dtmeta, dtsum[,c("chr", "cluster1", "genes")], by = c("chr", "cluster1"), all.x = T)

write_csv(dtmeta, file = "results/combined_outliers/outliers_n5_summary_genenames.csv")
