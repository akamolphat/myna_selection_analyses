library(data.table)
library(rehh)
chr_order <- c("chr1", "chr1A", "chr2", "chr3", "chr4", "chr4A", "chr5a", "chr5b", "chr5c", "chr6", "chr7", "chr8", "chr9", 
               "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
               "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr29", 
               "chr30", "chr31", "chr32")
ctr <- 1
for(i in chr_order) {
  # perform scan on a single chromosome (calculate iHH values)
  infile <- paste("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/", i, "_scanhh.txt", sep = "")
  scan <- fread(infile)
  rownames(scan) <- scan$V1
  scan$V1 <- NULL
  if (ctr == 1) {
    wgscan <- scan
  } else {
    wgscan <- rbind(wgscan, scan)
  }
  ctr <- ctr + 1
}
# calculate genome-wide iHS values
wgscan.ihs <- ihh2ihs(wgscan)

dtx <- wgscan.ihs$ihs
dtx$ID <- row.names(dtx)
dtx$FDR_qval <- p.adjust(p = (10^(-dtx$LOGPVALUE)), method = "fdr")

IHSout <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/WGS_IHS.txt"
fwrite(dtx, file = IHSout, sep = "\t", quote = F)


range(dtx$FDR_qval[dtx$LOGPVALUE > 6], na.rm = T) 
range(dtx$FDR_qval[dtx$LOGPVALUE <= 6], na.rm = T)


