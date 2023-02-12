### Ariana Tumor Recurrence Project
### Aime for mrna and cna data, confounders = ER status, n_components = 150
### setwd("/Users/Raghav/Pine Biotech/Research Fellows/Ariana/ISEF 2023/AIME-related-main/BRCA_AR_main/")


library(AIME)
library(MLmetrics)
library(pROC)
library(keras)
library(ggplot2)
library(tidyverse)
library(e1071)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("CancerSubtypes")
library(CancerSubtypes)

## Set the working directory to the data folder
#setwd("C:/Users/miz_r/Documents/Ariana/ScienceFair2022-2023/BRCA_AR/BRCA_AR")
setwd("/Users/Raghav/Pine Biotech/Research Fellows/Ariana/ISEF 2023/AIME-related-main/BRCA_AR_main/")
temp.path<-"/Users/Raghav/Pine Biotech/Research Fellows/Ariana/ISEF 2023/AIME-related-main/BRCA_AR_main/temp"
source('fdrgamma.r')

### Load gene expression table
gene <- read.table("data_mrna_select_brca.txt", header = TRUE, sep = "\t", row.names = 1)
gene_sel = FSbyVar(gene, cut.type = "topk", 1000)
gene2<-t(gene_sel)

### Load copy number variation table
cna <- read.table("data_cna_select_brca.txt", header = TRUE, sep = "\t", row.names = 1)
cna_sel = FSbyVar(cna, cut.type = "topk", 1000)
cna2<-t(cna_sel)

dim(gene2)
dim(cna2)

## z-score normalization
for(i in 1:ncol(gene2)) gene2[,i]<-(gene2[,i]-mean(gene2[,i]))/sd(gene2[,i])
for(i in 1:ncol(cna2)) cna2[,i]<-(cna2[,i]-mean(cna2[,i]))/sd(cna2[,i])

#####  
clin <- read.table("clinical_data_select.txt", header = TRUE, sep = "\t", row.names = 1)

ER<-rep(NA, nrow(clin))
ER[clin$ER.Status == "Positive"]<-1
ER[clin$ER.Status == "Negative"]<-0

conf<-cbind(ER)

sel<-which(apply(is.na(conf),1,sum) == 0)

all.colors<-c("red","green","blue","cyan", "yellow","orange","grey50", "red","white","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4")

g<-aime.select(data.in=gene2[sel,], data.out=cna2[sel,], confounder=conf[sel,], all.in.layers=2:5, all.out.layers=2:5, all.dropouts = c(0.1,0.2,0.3), ncomp = 150, repeats=3, col=all.colors[as.numeric(as.factor(ER))], cor.cut=0.5, kurtosis.cut=0.5, skew.cut=0.5)


rec<-new("list")
for(n in 1:10)
{
  b<-aime(data.in =gene2, data.out=cna2, in.layers=5, out.layers=4,  max.epochs=100, max.dropout=0.2, importance.permutations=2, ncomp=150, pairwise.importance = TRUE)
  rec[[n]]<-b
  save(rec, file="mrna to cna aime layers 4 5 dropout 0.2 no confounder.bin")
}

rec<-new("list")
for(n in 1:10)
{
  b<-aime(data.out=cna2, data.in=gene2, confounder=conf, in.layers=5, out.layers=4,  max.epochs=100, max.dropout=0.2, importance.permutations=2, ncomp=150, pairwise.importance = TRUE)
  rec[[n]]<-b
  save(rec, file="mrna to cna aime layers 5 4 dropout 0.2 with confounder.bin")
}


### Load AIME output
load("mrna to cna aime layers 5 4 dropout 0.2 with confounder.bin")

### Verify PCA output for confounder factor adjustment
embed1.pca <- prcomp(rec[[1]][["embeded"]])
p1 <- autoplot(embed1.pca, main="PCA for mrna-cna embedded data", 
               data=clin, 
               shape ="Disease.Free.Status", colour = "ER.Status.By.IHC", size = 3.5) + theme_light(base_size = 14)
png("aime_output_pca.png", width = 8, height = 6, units = "in", res = 300)
p <- plot_grid(p1, ncol=1)
plot_grid(p, ncol=1) # rel_heights values control title margins
ttemp <- dev.off()

