library(gdsfmt)
library(SNPRelate)
library(data.table)
setwd("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/")

# Thinned WGS data -------------------------------------------------------------
#biallelic by default
snpgdsVCF2GDS("variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.5kbthinned.snps.vcf.gz", 
              "variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.5kbthinned.gds")
snpgdsSummary("variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.5kbthinned.gds")
genofile = snpgdsOpen("variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.5kbthinned.gds")

sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Read popmap for colouring 
popmap <- fread("pop_map_WGS_ALL.txt", col.names = c("ID", "pop"), header = F, skip = 0)
# popmap <- popmap[match(sample_id, popmap[["ID"]])]
# add.gdsn(genofile, "popmap", popmap)

# #LD based SNP pruning
# set.seed(1000)
# snpset = snpgdsLDpruning(genofile,ld.threshold = 0.5)
# snp.id=unlist(snpset)

popmap <- popmap[match(sample_id, popmap[["ID"]])]
# Make FST heatmap
# v <- snpgdsFst(genofile, population=as.factor(popmap$pop), method="W&C84", autosome.only = F)

# distance matrix - use IBS
dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL, 
                         snp.id = NULL,
                         autosome.only = F,
                         remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)

# dissMatrix$ibs[pop.idx, pop.idx]

snpgdsClose(genofile)

snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)
popmap <- popmap[match(snpHCluster$sample.id, popmap[["ID"]])]

#samp.group=as.factor(pop_code)


cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group= as.factor(popmap$pop), 
                        col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=1:11,label.H=FALSE, label.Z=TRUE, 
                        verbose=TRUE)
png("WGS_tree.png", width = 11.7, height = 8.3, units = "in", res = 600)
# par(xpd=TRUE)
# par(mar = c(6, 5, 4, 6.5))
snpgdsDrawTree(cutTree, main = "WGS",
               type = "dendrogram",
               edgePar=list(col=rgb(0.5,0.5,0.5,0.75),
                            t.col="black"),
               y.label.kinship=F,
               leaflab="none")

legend("bottomright", legend=levels(as.factor(popmap$pop)), col=1:nlevels(as.factor(popmap$pop)), 
       pch=1:11, ncol=2)

dev.off()



# key outlier region -----------------------------------------------------------
snpgdsVCF2GDS("variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.20625-20699kb.snps.vcf.gz", 
              "variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.20625-20699kb.gds")
snpgdsSummary("variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.20625-20699kb.gds")
genofile = snpgdsOpen("variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.20625-20699kb.gds")

sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Read popmap for colouring 
popmap <- fread("pop_map_WGS_ALL.txt", col.names = c("ID", "pop"), header = F, skip = 0)
# popmap <- popmap[match(sample_id, popmap[["ID"]])]
# add.gdsn(genofile, "popmap", popmap)

# #LD based SNP pruning
# set.seed(1000)
# snpset = snpgdsLDpruning(genofile,ld.threshold = 0.5)
# snp.id=unlist(snpset)

# distance matrix - use IBS
dissMatrix2  =  snpgdsIBS(genofile , sample.id=NULL, 
                         snp.id = NULL,
                         autosome.only = F,
                         remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)
snpgdsClose(genofile)

snpHCluster2 =  snpgdsHCluster(dissMatrix2, sample.id=NULL, need.mat=TRUE, hang=0.01)
popmap <- popmap[match(snpHCluster2$sample.id, popmap[["ID"]])]

# samp.group=as.factor(pop_code)

cutTree2 = snpgdsCutTree(snpHCluster2, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group= as.factor(popmap$pop), 
                        col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=1:11,label.H=FALSE, label.Z=TRUE, 
                        verbose=TRUE)
png("keyoutlier_tree.png", width = 11.7, height = 8.3, units = "in", res = 600)

snpgdsDrawTree(cutTree2, main = "Key outlier region",
               type = "dendrogram",
               edgePar=list(col=rgb(0.5,0.5,0.5,0.75),
                            t.col="black"),
               y.label.kinship=F,
               leaflab="none")

# legend("bottom", legend=levels(as.factor(popmap$pop)), col=1:nlevels(as.factor(popmap$pop)), 
#        pch=1:11, ncol=3)


dev.off()

labels <- gsub("_", " ", levels(as.factor(popmap$pop)))
labels[labels == "Maharashtra subpopulation A"] <- "Maharashtra subpop. A"

png("WGS_vs_keyoutlier_tree.png", height = 10, width = 8.3, units = "in", res = 600)
par(mfrow = c(2,1), mai = c(0.1, 0.1, 0.1, 0.1))

snpgdsDrawTree(cutTree, main = "WGS",
               type = "dendrogram",
               edgePar=list(col=rgb(0.5,0.5,0.5,0.75),
                            t.col="black"),
               y.label.kinship=F,
               leaflab="none")

legend("bottomright", legend=labels, col=1:nlevels(as.factor(popmap$pop)), 
       pch=1:11, ncol=2)

# keyoutlier region
snpgdsDrawTree(cutTree2, main = "Key outlier region",
               type = "dendrogram",
               edgePar=list(col=rgb(0.5,0.5,0.5,0.75),
                            t.col="black"),
               y.label.kinship=F,
               leaflab="none")
dev.off()


