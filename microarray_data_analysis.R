# This R script can be used for microarray data analysis.
#
# Date: 2018/11/13
# Upgrade: 2021/1/31
# Author: C.W.Weng(Jeff)

### install packages
install.packages("gplots")

### use the bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

### install additional bioconductor libraries
BiocManager::install("GEOquery")
BiocManager::install("affy")
BiocManager::install("simpleaffy")
BiocManager::install("affyPLM")
BiocManager::install("affydb")
BiocManager::install("limma")
BiocManager::install("gcrma")
BiocManager::install("annotate")
BiocManager::install("hgu133a.db")
BiocManager::install("hgu133b.db")
BiocManager::install("hgu133plus2.db")
BiocManager::install("hugene10stv1cdf")
BiocManager::install("hugene10stv1probe")
BiocManager::install("hugene10stprobeset.db")
BiocManager::install("hugene10sttranscriptcluster.db")

### load the necessary libraries
library(GEOquery)
library(affy)
library(simpleaffy)
library(affyPLM)
library(affydb,character.only=T)
library(limma)
library(gcrma)
library(annotate)
library(hgu133a.db)
library(hgu133b.db)
library(hgu133plus2.db)
library(hugene10stv1cdf)
library(hugene10stv1probe)
library(hugene10stprobeset.db)
library(hugene10sttranscriptcluster.db)
library(gplots)

### set working directory
setwd("/Users/cwweng/Documents/Coding/R/microarray_data_analysis")

### download the CEL file package by GEO Series ID (GSE - Prefix of ID)
getGEOSuppFiles(GEO = "GSE21363", baseDir = getwd())

### unpack the CEL files
untar("GSE21363/GSE21363_RAW.tar", exdir = "data")
cels <- list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)
celpath <- paste(getwd(), "data", sep="/")
setwd(celpath)

### process the CEL files
celfiles <- c("GSM533844.CEL", "GSM533845.CEL", "GSM533846.CEL", "GSM533847.CEL", "GSM533848.CEL", "GSM533849.CEL")
raw.data <- ReadAffy(filenames = celfiles)
pData(raw.data)$Treatment <- rep(c("Day0","Day8"),each=3)

### quality control
# boxplot for the base 2 logarithm of expression values
n.cel  <- length(celfiles)
cols <- rainbow(n.cel*1.2)
boxplot(raw.data, col = cols, xlab = "Array Index", ylab="log2(density)")
# compare the original expressions and the logarithm scale
hist(raw.data, lty=1:3, col=cols)
legend("topright", legend=sampleNames(raw.data), lty=1:3, col=cols, box.col="transparent", xpd=T)
# check the control markers
data.qc <- qc(raw.data)
plot(data.qc)
# observe the RNA degradation plot
data.deg <- AffyRNAdeg(raw.data)
plotAffyRNAdeg(data.deg, col = cols)
legend("topleft", rownames(pData(raw.data)), col=cols, lwd=1, inset = 0.05)

### preprocessing and normalization
# normalize by GCRMA algorithm
data.gcrma <- gcrma(raw.data)
eset <- exprs(data.gcrma)
boxplot(eset, col = cols, main = "gcRMA")
# normalize by RMA algorithm
data.rma <- rma(raw.data)
eset <- exprs(data.rma)
boxplot(eset, col = cols, main = "RMA")
# quantile normalization
eset <- normalizeQuantiles(eset)
boxplot(eset, col = cols, main = "Quantile")

### identify the differential expressed genes
conditions_ <- factor(raw.data$Treatment)
design <- model.matrix(~-1+conditions_)
contrast.matrix <- makeContrasts(contrasts = "conditions_Day8 - conditions_Day0", levels = design)
fit <- lmFit(eset, design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
dif <- topTable(fit2, coef="conditions_Day8 - conditions_Day0", n = nrow(fit2), lfc = log2(1.5))
dif1 <- dif[dif[,"adj.P.Val"]<0.01,]

### annotate the differential expressed genes
affydb <- annPkgName(raw.data@annotation, type="db")
dif1$symbols <- getSYMBOL(rownames(dif1), affydb)
dif1$EntrezID <- getEG(rownames(dif1),affydb)
dif2 <- dif1[(!is.na(dif1$symbols)),]

### hierarchical clustering of the differential expressed genes
eset1 <- eset[(rownames(dif2)),]
row.names(eset1) <- dif2$symbols
heatmap.2(as.matrix(eset1), col = redgreen(75), cexRow = 0.2, cexCol = 0.5, scale ="row", trace = "none", key = T, keysize = 1.2, density.info = "none")

# output the expression profile of differential expressed genes
write.table(eset1, file="DEGs_Table.txt", quote=F, sep="\t")

# output the annotation profile of differential expressed genes
write.table(dif2, file="DEGs_annotation_Table.txt", quote=F, sep="\t")
