dtx <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/WGS_IHS.txt")
png("results/EHHS/WGS_IHS_pvalue_histogram.png", width = 7.48031, height = 3.5, units = "in", res = 600)
hist((10^(-dtx$LOGPVALUE)), xlab = expression(~italic(p)*"-value"), main = expression('iHS'~italic(p)*'-value histogram'))
dev.off()

min(dtx$FDR_qval[dtx$LOGPVALUE >= 6], na.rm = T)

