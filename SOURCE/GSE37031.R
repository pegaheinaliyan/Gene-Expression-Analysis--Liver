
setwd("C:/Users/roozane/Desktop/liver/GSE37031")
library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)




series <- "GSE37031"
platform <- "GPL14788"

#download data
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")

####choose platform
if (length(gset) > 1) idx <- grep(platfrom, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] 
#gset[[1]]


#grouping data
gr <- c(rep("NASH",8),rep("Control",7))

#expression matrix 
ex <- exprs(gset)

#delete rows that include NA
ex<-na.omit(ex)



#log2 if requierd
#ex <- log2(ex+1) 
#exprs(gset) <- ex

#boxplot
pdf("Results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off()

#normaize if requird
#ex <- normalizeQuantiles(ex)
#exprs(gset) <- ex

#correlation HeatMap
pdf("Results/CoreHeatmap.pdf", width = 15, height = 15)
pheatmap(cor(ex), labels_row = gr, labels_col = gr, color = greenred(256), border_color = NA)
dev.off()

#principle component analysis(PCA)
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

#subtract mean from each Gen too see the variation
#func scale works on column so transport the ex
ex.scale <- t(scale(t(ex),scale = F))
pc <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc$x[,1:2])
dev.off()

# PCR for Samples
pcr <- data.frame(pc$r[,1:3],Group=gr) ##the samples are in rotation(pc$r)
pdf("Results/PCA_samples.pdf")
ggplot(pcr,aes(PC1, PC2, color = Group)) + geom_point(size=3) +theme_bw()
dev.off()

#Diferrential Expression Analysis
#we needs category instead of string
gr <- factor(gr)
gset$description <- gr

design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gr)

fit <- lmFit(gset, design)

cont.matrix <- makeContrasts(Control-NASH,levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

write.table(tT,"Results/Control-NASH(all).txt", row.names=F, sep="\t",quote = F)



