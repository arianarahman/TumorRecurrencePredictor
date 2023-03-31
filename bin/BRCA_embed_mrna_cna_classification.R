### Ariana Tumor Recurrence Project
### Load Aime output for mrna and cna data, confounders = ER status, n_components = 150
### Random Forest Classification algorithm comparisons
### setwd("C:/Users/arian/Documents/git/TumorRecurrencePredictor/data")


### Classification Algorithm
library(preprocessCore)
library(caret)
library(randomForest)
library(BBmisc)
library(pROC)

### set the working directory to the data folder
setwd("C:/Users/arian/Documents/git/TumorRecurrencePredictor/data")

#####  Load clinical data and sample order from gene expression data
clin <- read.table("clinical_data_select.txt", header = TRUE, sep = "\t", row.names = 1)

### Load gene expression table
gene <- read.table("data_mrna_select_brca.txt", header = TRUE, sep = "\t", row.names = 1)
### Select the same input genelist
gene_sel = FSbyVar(gene, cut.type = "topk", 1000)
gene2<-t(gene_sel)
gene2 <- as.data.frame(gene2)
### preprocess the sample name format to match clinical data
rownames(gene2) <- gsub("\\.","-",rownames(gene2))
### add group annotation to expression table
gene2 <- cbind(group = clin$Disease.Free.Status, gene2[row.names(clin),])
gene2 <- na.omit(gene2)
gene2$group <- gsub("/Progressed","",gsub("[0-1]:","",gene2$group))
gene2$group <- as.factor(gene2$group)


### Load copy number variation data
cna <- read.table("data_cna_select_brca.txt", header = TRUE, sep = "\t", row.names = 1)
### Select the same input genelist
cna_sel = FSbyVar(cna, cut.type = "topk", 1000)
cna2<-t(cna_sel)
cna2 <- as.data.frame(cna2)
### preprocess the sample name format to match clinical data
rownames(cna2) <- gsub("\\.","-",rownames(cna2))
### add group annotation to cna table
cna2 <- cbind(group = clin$Disease.Free.Status, cna2[row.names(clin),])
cna2 <- na.omit(cna2)
cna2$group <- gsub("/Progressed","",gsub("[0-1]:","",cna2$group))
cna2$group <- as.factor(cna2$group)

### Load copy number variation data
mirna <- read.table("data_mirna_select_brca.txt", header = TRUE, sep = "\t", row.names = 1)
mirna2<-t(mirna)
mirna2 <- as.data.frame(mirna2)
### preprocess the sample name format to match clinical data
rownames(mirna2) <- gsub("\\.","-",rownames(mirna2))
### add group annotation to mirna table
mirna2 <- cbind(group = clin$Disease.Free.Status, mirna2[row.names(clin),])
mirna2 <- na.omit(mirna2)
mirna2$group <- gsub("/Progressed","",gsub("[0-1]:","",mirna2$group))
mirna2$group <- as.factor(mirna2$group)

### Load AIME output  without confounder adjustment
load("mrna to cna aime layers 4 5 dropout 0.2 no confounder.bin")
## create dataframe from embedding and preprocess the output
embed1 <- as.data.frame(rec[[1]][["embeded"]])
### Preprocess the embedded dataframe for random forest
aime_no_conf <- cbind(row.names = rownames(clin),
                  group = clin$Disease.Free.Status,
                  embed1)
aime_no_conf <- na.omit(aime_no_conf)
aime_no_conf$group <- gsub("/Progressed","",gsub("[0-1]:","",aime_no_conf$group))
aime_no_conf$group <- as.factor(aime_no_conf$group)

### Change the order of samples, Disease free followed by Recurred
aime_no_conf <- rbind(
  cbind(group = aime_no_conf[aime_no_conf$group == "DiseaseFree",1],
        aime_no_conf[aime_no_conf$group == "DiseaseFree",-1]),
  cbind(group = aime_no_conf[aime_no_conf$group == "Recurred",1],
        aime_no_conf[aime_no_conf$group == "Recurred",-1])
)

### Load AIME output  with confounder adjustment
load("mrna to cna aime layers 5 4 dropout 0.2 with confounder.bin")
## create dataframe from embedding and preprocess the output
embed1 <- as.data.frame(rec[[1]][["embeded"]])
### Preprocess the embedded dataframe for random forest
aime_with_conf <- cbind(row.names = rownames(clin),
                  group = clin$Disease.Free.Status,
                  embed1)
aime_with_conf <- na.omit(aime_with_conf)
aime_with_conf$group <- gsub("/Progressed","",gsub("[0-1]:","",aime_with_conf$group))
aime_with_conf$group <- as.factor(aime_with_conf$group)

### Heatmap
load("mrna to cna aime layers 5 4 dropout 0.2 with confounder.bin")
feature_importance <- as.integer(unlist(sort(rec[[1]]$imp,decreasing = TRUE, index.return=TRUE)))

### Plot heat map of transcriptomic profiles for visual verification
# scale the training data for the heatmap
gene2.data.scaled <- as.matrix(gene2[,(feature_importance[1:25]+1)]) %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

gene2.data.scaled <- rbind(gene2.data.scaled[gene2$group == "DiseaseFree",-1],
                           gene2.data.scaled[gene2$group == "Recurred",-1] + 0.7)

#Set annotation
library(ComplexHeatmap)
ann = data.frame(MRNA=aime_with_conf$group)
colours <- list('MRNA' = c('Recurred' = 'red', 'DiseaseFree' = 'limegreen'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),
                            #annotation_name_side = "left",
                            show_annotation_name =FALSE,
                            show_legend = FALSE)

split <- factor(aime_with_conf$group, levels=c("DiseaseFree","Recurred"))
hmap1 <- Heatmap(
  t(gene2.data.scaled),
  name = "Important Genes (AIME)",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_dend = TRUE,
  show_row_dend = TRUE,
  row_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  #column_title = "Train Data",
  column_split=split, cluster_row_slices = FALSE,
  width = unit(100, "mm"),
  top_annotation=colAnn)

png("Important_genes_heatmap.png", width = 8, height = 6, units = "in", res = 300)
draw(hmap1)
dev.off()


### Change the order of samples, Disease free followed by Recurred
aime_with_conf <- rbind(
  cbind(group = aime_with_conf[aime_with_conf$group == "DiseaseFree",1],
        aime_with_conf[aime_with_conf$group == "DiseaseFree",-1]),
  cbind(group = aime_with_conf[aime_with_conf$group == "Recurred",1],
        aime_with_conf[aime_with_conf$group == "Recurred",-1] + 0.7)
)


### Define a function RF_classifier for multiple calls
RF_classifier <- function(input_df){
  ### make a 70/30 split of train and test data
  set.seed(3456)
  train.index <- createDataPartition(input_df$group, p = .6, list = FALSE)
  train <- input_df[train.index,]
  test  <- input_df[-train.index,]
  
  train.data = train[,-1] #Set the training set
  test.data = test[,-1] #Set the testing set
  
  train.output = train[,1]  #store the labels of train data
  test.output = test[,1]  # store the true labels of test data
  
  ### Fitting Random Forest to the train dataset
  set.seed(234)
  classifier_RF = randomForest(x = train.data,
                               y = train.output,
                               ntree = 500,
                               importance = TRUE)
  
  ### Print out details of the classiier
  #print(classifier_RF)
  
  ### Draw plots of accuracy
  #plot1 <- plot (classifier_RF, main='Random forest plot', lty=1, cex=1.5)
  #print(plot1)
  
  ### Prediction on test.data as a whole
  test.predict = predict(classifier_RF, test.data)
  
  ### predicting probabilities for ROC
  test.predict.prob = predict(classifier_RF, test.data, type = "prob")
  test_roc = roc(test.output, test.predict.prob[,1], levels=c("DiseaseFree", "Recurred"))

  result <- confusionMatrix(as.factor(test.predict), test.output,
                            positive = "DiseaseFree")
  
  feat_imp_df <- importance(classifier_RF) %>% 
    data.frame() %>% 
    mutate(feature = row.names(.)) %>% arrange(desc(MeanDecreaseAccuracy)) 
  
  return(list(classifier_RF, result, feat_imp_df, test_roc))    
}

df_list <- list(gene2, cna2, mirna2, aime_no_conf, aime_with_conf)
names(df_list) <- c("mrna", "cna", "mirna", "aime no adjustment", "aime confounder adjusted")
comment(gene2) = "mrna"
comment(cna2) = "cna"
comment(mirna2) = "mirna"
comment(aime_no_conf) = "aime no adjustment"
comment(aime_with_conf) = "aime confounder adjusted"

color_roc = c("#7570B3", "#E7298A", "#66A61E", "#D95F02", "#1B9E77")

RF_accuracy_metric = c()
i = 1
feat <- list()
RF_roc_metric = c()
ctr = 0

for (inputdf in df_list){
  ### ctr 
  ctr = ctr + 1
  ### call the RF classifier
  RF_full <- RF_classifier(inputdf)
  
  ### confusion matrix results
  result <- RF_full[[2]]
  
  ### Feature importance results
  feat_imp_df <- RF_full[[3]]
  
  ### Filter genes/features with zero importance
  feat_imp_df[feat_imp_df$MeanDecreaseAccuracy>0, ]
  top_sel_feat <- factor(feat_imp_df[feat_imp_df$MeanDecreaseAccuracy>0, ]$feature)
  No_top_sel_feat <- length(top_sel_feat)
  
  ### collect the accuracy metrics
  temp_acc_metric <- cbind(Input = comment(inputdf),
                           Accuracy_DiseaseFree = as.numeric(result$byClass["Pos Pred Value"]),
                           Accuracy_Recurred = as.numeric(result$byClass["Neg Pred Value"]),
                           Accuracy_over_all = as.numeric(result$byClass["Balanced Accuracy"]), 
                           Kappa = as.numeric(result$overall[[2]]))
  
  RF_accuracy_metric <- rbind(temp_acc_metric, RF_accuracy_metric)
  RF_roc_metric <- rbind(RF_roc_metric, RF_full[[4]])
  
  #### plotROC curves
  if(ctr<2){
    png("ROC_AUC_stats.png", width = 6, height = 6, units = "in", res = 300)
    plot.roc(RF_full[[4]], las=1, lwd=1.7, col = color_roc[ctr], xlab=" 1 - specificity")
  }else{
    plot.roc(RF_full[[4]], add=TRUE, lwd=1.7, col = color_roc[ctr])
  }
  
  ### print the accuracy metrics
  #print (temp_acc_metric)
}
#Insert a legend
legend("bottomright", legend = c("Gene expression: AUC = 46%",
                                 "Copy number variation: AUC = 51%",
                                 "miRNA expression: AUC = 65%",
                                 "aime no adjustment: AUC = 66%",
                                 "aime confounder adjusted: AUC = 100%"
                                 ), 
       lty = c(1, 1, 1, 1, 1, 1), 
       col =c("#7570B3", "#E7298A", "#66A61E", "#D95F02", "#1B9E77"), cex=0.75, lwd=1.5,
       title = "ROC Curves")
par(mfrow = c(1, 1))
ttemp <- dev.off()


RF_accuracy_metric <- cbind(Input = RF_accuracy_metric[,1], 
                            as.data.frame(RF_accuracy_metric) %>% 
                              select(-Input) %>% 
                              mutate_if(is.character,as.numeric)
                            )


png("Classifier_accuracy_comparison.png", width = 8, height = 6, units = "in", res = 300)
library(reshape2)
ggplot(data=melt(dplyr::select(RF_accuracy_metric, c(Input, Accuracy_over_all))), aes(x=factor(Input, levels = c("mrna", "cna", "mirna", "aime no adjustment", "aime confounder adjusted")), y=value*100, fill=Input)) +
  geom_bar(stat="identity") +
  ggtitle("Accuracy for different inputs") +
  xlab("Input")+
  scale_x_discrete(labels=c("mrna" = "Gene expression",
                            "cna" = "Copy number variation",
                            "mirna" = "miRNA expression",
                            "aime no adjustment" = "AIME no adjustment",
                            "aime confounder adjusted" = "AIME confounder adjusted"
                            )) +
  ylab("Balanced Accuracy (%)") +
  ylim(0,100) +
  theme_light(base_size = 14) + 
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
