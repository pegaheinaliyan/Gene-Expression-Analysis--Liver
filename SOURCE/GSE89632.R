setwd("C:/Users/roozane/Desktop/liver/GSE89632")
library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)

series <- "GSE89632"


#download data
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")

#choose platform
if (length(gset) > 1) idx <- grep(platfrom, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] 
#gset[[1]]


#grouping data
gsms <- paste0("232222223323322333232333222232333332300000000000000002000200000")

groups <- c()
for (i in 1:nchar(gsms)) { 
  if (substr(gsms,i,i) == "0") {
    groups[i] <- "NOH"
  }else if(substr(gsms,i,i) == "1"){
    groups[i] <- "OH"
  }else if(substr(gsms,i,i) == "2"){
    groups[i] <- "ST"
  }else if(substr(gsms,i,i) == "3"){
    groups[i] <- "NASH"
  }
}


#expression matrix 
expMat <- exprs(gset)

#log2 if requierd
#ex <- log2(ex+1) 
#exprs(gset) <- ex
#pdf("Results/boxPlot.pdf", width =  170)
png("Results/boxPlot.png", width =  3000, height = 500)
boxplot(expMat)
dev.off()

#----HeatMap
#pdf("Results/CorrolationHeatMap.pdf", width =15, height = 15)
png("Results/CorrolationHeatMap.png", width =2000, height = 2000)
pheatmap(cor(expMat), labels_row = groups, labels_col = groups, fontsize = 20)
dev.off()

#----PCA
pca <- prcomp(expMat)
#pdf("Results/PCA.pdf")
png("Results/PCA_Variances.png", 1000,1000)
plot(pca)
dev.off()

png("Results/PCA.png",1000,1000)
plot(pca$x[,1:2])
dev.off()

#----PCA-Scaled
expMat_scaled <- t(scale(t(expMat) , scale = FALSE))
pca <- prcomp(expMat_scaled)
#pdf("Results/PCA_Scaled.pdf")
png("Results/PCA_Scaled_Variances.png", 1000,1000)
plot(pca)
dev.off()
png("Results/PCA_Scaled.png",1000,1000)
plot(pca$x[,1:2])
dev.off()

pcr <- data.frame(pca$rotation[,1:3], Group = groups)
#pdf("Results/PCA_Samples.pdf")
png("Results/PCA_Samples.png", 2000, 2000)
ggplot(pcr, aes(PC1, PC2, color = Group )) + geom_point(size = 10) + scale_color_brewer(palette = "Paired") + theme_bw() + theme(legend.text=element_text(size=rel(1.5)),text = element_text(size=30))
dev.off()

###### Differential Expression Analysis

#setting new Groups
#NOH 0 - OH 1 - ST 2 - NASH 3
analyzeGroup1 <- "0"
analyzeGroup2 <- "3"
txtname <- paste("Results/",series,"-",analyzeGroup1,"-",analyzeGroup2,".txt", sep = "")


sml <- c()
for (i in 1:nchar(gsms)) { 
  if (substr(gsms,i,i) == analyzeGroup1 || substr(gsms,i,i) == analyzeGroup2){
    sml[i] <- substr(gsms,i,i)  
  }else{
    sml[i] <- "X"
  }
}
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]


groups <- c()

for (i in 1:length(sml)) { 
  if (sml[i] == "0") {
    groups[i] <- "NOH"
  }else if(sml[i] == "1"){
    groups[i] <- "OH"
  }else if(sml[i] == "2"){
    groups[i] <- "ST"
  }else if(sml[i] == "3"){
    groups[i] <- "NASH"
  }
}


groups <- factor(groups)
gset$description <- groups
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(groups)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(NOH - NASH, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number= Inf)
tT <- subset(tT, select=c("ID", "Symbol","adj.P.Val","logFC", "P.Value"))
tT <- na.omit(tT)
write.table(tT, txtname, row.names=F, sep="\t" , quote = F)

###################BIO 12
#PathWay Analysis
aml.up <- subset(tT, logFC >1 & adj.P.Val <0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Symbol,"///")))
upFileName <- paste("Results/",series,"-",analyzeGroup1,"-",analyzeGroup2,"-UP.txt", sep = "")
write.table(aml.up.genes, file = upFileName, quote = F, row.names = F, col.names = F)

aml.down <- subset(tT, logFC < -1 & adj.P.Val <0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Symbol,"///")))
downFileName <- paste("Results/",series,"-",analyzeGroup1,"-",analyzeGroup2,"-DOWN.txt", sep = "")
write.table(aml.down.genes, file = downFileName, quote = F, row.names = F, col.names = F)

