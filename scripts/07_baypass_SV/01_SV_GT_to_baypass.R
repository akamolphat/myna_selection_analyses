library(tidyverse)
library(data.table)
dt <- read_csv("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/Potential_repeats/repeat_ID_482137/Genotyping_SV_manually.csv")
dtsum <- dt %>%
  group_by(pop) %>%
  summarise(insertion = sum(SV_insertion_no, na.rm = T),
            deletion = sum((2-SV_insertion_no), na.rm = T)) 

matsum <- as.matrix(dtsum[,c("insertion", "deletion")])

mat_d <- matrix(t(matsum), nrow = 1, byrow = T)

fwrite(mat_d, file = "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/myna_baypass_SV.txt", sep = " ", row.names = F, col.names = F)
