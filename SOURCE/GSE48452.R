setwd("C:/Users/roozane/Desktop/liver/GSE48452")
library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)

series <- "GSE48452"
platform <- "GPL11532"

#download data
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")

#choose platform
if (length(gset) > 1) idx <- grep(platfrom, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] 
#gset[[1]]


#grouping data
gr <- c(rep("Control",5),"NASH","HealthyObse",rep("Control",3),"NASH",rep("HealthyObse",5),"NAFLD","NASH",rep("NAFLD",2),rep("HealthyObse",3),"NAFLD",
        "HealthyObse",rep("NASH",2),"HealthyObse","Control","NAFLD","HealthyObse",rep("NASH",4 ),"HealthyObse","NASH","HealthyObse","NASH",
        rep("Control",2),"NAFLD",rep("HealthyObse",3),"NASH","HealthyObse","NASH","Control",rep("HealthyObse",2),"NAFLD","HealthyObse","NAFLD",
        "Control","NAFLD","NASH","NAFLD","HealthyObse","NAFLD","HealthyObse","Control","HealthyObse","NASH","NAFLD",rep("NASH",2),"NAFLD","HealthyObse",
        "NASH","NAFLD",rep("HealthyObse",2))

#expression matrix 
ex <- exprs(gset)

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

write.table(tT,"Results/Control-Nash.txt", row.names=F, sep="\t",quote = F)
