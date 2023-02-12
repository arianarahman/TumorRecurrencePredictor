### Ariana Tumor Recurrence Project
### BRCA datasets from CBioportal exploratory analysis
### setwd("/Users/Raghav/Pine Biotech/Research Fellows/Ariana/ISEF 2023/Data Exploration")

library(tidyverse)
library(ggfortify)
library(VennDiagram)
library(cowplot)

### Go the CBioportal and download clinical data for patient ID common to RNASEQ, Mutation, Copy Number Variation, Methylation data
### https://bit.ly/3JJ6IpK

### set working directory to data folder
setwd("/Users/Raghav/Pine Biotech/Research Fellows/Ariana/ISEF 2023/Data Exploration/PCA_unified/")

### Load the common clinical data
clin.data.tcga.legacy.comm.all.data <- read_tsv("brca_tcga_clinical_data.tsv")

# rename colnames
colnames(clin.data.tcga.legacy.comm.all.data) <- make.names(names(clin.data.tcga.legacy.comm.all.data))

# Dimension
dim(clin.data.tcga.legacy.comm.all.data)

# patient distribution for clinical characteristics disease free/recurred
### Remove A,B,C sub lineation of stages and space between Stage and I, II, III
clin.data.tcga.legacy.comm.all.data$Derived.Tumor.Stage = gsub(" |A|B|C", "",clin.data.tcga.legacy.comm.all.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)
clin_char1 <- table(clin.data.tcga.legacy.comm.all.data$Disease.Free.Status, clin.data.tcga.legacy.comm.all.data$Derived.Tumor.Stage)

png("patient_distribution_tum_stage.png", width = 8, height = 6, units = "in", res = 300)
par(mar=c(5.1,5.1,4.1,2.1))
barplot(clin_char1,
        main = "Tumor Recurrence by Disease Stage",
        xlab = "No. Patients",
        xlim = c(0,350),
        col = c("green","red"),
        density=25,
        las = 2,
        horiz = TRUE
)
legend("topright",
       c("Disease Free","Recurred"),
       fill = c("green","red"),
       density = 25
)
tmp <- dev.off()


### subset only necessary columns
clin.data.tcga.legacy.comm.all.data.select <- clin.data.tcga.legacy.comm.all.data %>% dplyr::select(Sample.ID,
                                                                                                    Patient.ID,
                                                                                                    Disease.Free..Months.,
                                                                                                    Disease.Free.Status,
                                                                                                    Year.Cancer.Initial.Diagnosis,
                                                                                                    Derived.Tumor.Stage,
                                                                                                    ER.Status.By.IHC)

write.table(clin.data.tcga.legacy.comm.all.data.select, file = "clinical_data_select.txt", quote = FALSE, row.names = FALSE, sep = "\t")

### convert the patient id format in clinical data to meet mrna data format
clin.data.tcga.legacy.comm.all.data.select$Patient.ID <- paste0(clin.data.tcga.legacy.comm.all.data.select$Patient.ID, "-01")
clin.data.tcga.legacy.comm.all.data$Patient.ID <- paste0(clin.data.tcga.legacy.comm.all.data$Patient.ID, "-01")


## Load all omics data downloded from Cbioportal
mrna.data <- read_tsv("data_mrna_seq_v2_rsem_brca.txt")
mrna.data <- mrna.data %>% dplyr::select(-Entrez_Gene_Id)
mirna.data <- read_tsv("data_mirnaseq_firebrowse_brca.txt")
cna.data <- read_tsv("data_linear_cna_brca.txt")
cna.data <- cna.data %>% dplyr::select(-Entrez_Gene_Id)
mut.data <- read_tsv("data_mutation_brca.txt")

### check the dimension of omics data
dim(mrna.data)
dim(mirna.data)
dim(cna.data)
dim(mut.data)

# rename colnames
colnames(mrna.data) <- make.names(names(mrna.data))
colnames(mirna.data) <- make.names(names(mirna.data))
colnames(cna.data) <- make.names(names(cna.data))
colnames(mut.data) <- make.names(names(mut.data))

### check column names
names(mrna.data)[1:10]
names(mirna.data)[1:10]
names(cna.data)[1:10]
names(mut.data)[1:10]

### reformat column names
colnames <- gsub("\\.","-",names(mrna.data)[-1])
names(mrna.data) <- c("Hugo_Symbol", colnames)
   
colnames <- gsub("\\.","-",names(mirna.data)[-1])
names(mirna.data) <- c("Hugo_Symbol", colnames)

colnames <- gsub("\\.","-",names(cna.data)[-1])
names(cna.data) <- c("Hugo_Symbol", colnames)

colnames <- gsub("\\.","-",names(mut.data)[-1])
names(mut.data) <- c("Hugo_Symbol", colnames)


### Patient sharing between omics data
venn.diagram(
  x = list(names(mrna.data)[-1], names(mirna.data)[-1], names(cna.data)[-1], names(mut.data)[-1]),
  category.names = c("mrna (1100)", "mirna (452)", "cna (1081)", "mut (982)"),
  filename = "Patient_distribution_omics_data.png",
  output = FALSE,
  imagetype="png" ,
  height = 800 ,
  width = 800 ,
  resolution = 800,
  compression = "lzw",
  lwd=0.2,
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
  alpha = 0.50,
  fontfamily = "sans",
  cex = 0.15,
  cat.cex = 0.12,
  cat.dist = c(0.15, 0.15, 0.1, 0.1),
  cat.pos = c(-60, 60, -30, 30),
  #cat.default.pos = "outer",
  disable.logging = TRUE,
  margin = 0.12
)

### Common patients tumor regression status
common_patients_omics_data <- Reduce(intersect, list(names(mrna.data)[-1], names(mirna.data)[-1], names(cna.data)[-1], names(mut.data)[-1]))

### clinical characteristics of common patients
plot_regression_status <- function(select_clin_char1, plot_file_name) {
  clin_char1 <- table(select_clin_char1$Disease.Free.Status, select_clin_char1$Derived.Tumor.Stage)
  print("Processing...")
  print(plot_file_name)
  print(clin_char1)
  par(mar=c(5.1,5.1,4.1,2.1))
  png(plot_file_name, width = 8, height = 6, units = "in", res = 300)
  barplot(clin_char1,
          main = "Tumor Recurrence by Disease Stage",
          xlab = "No. Patients",
          xlim = c(0,plyr::round_any(max(clin_char1), 50, f = ceiling)),
          col = c("green","red"),
          density=25,
          las = 2,
          horiz = TRUE
  )
  legend("topright",
         c("Disease Free","Recurred"),
         fill = c("green","red"),
         density = 25
  )
  tmp <- dev.off()
}

### subset omics data for selected patients
clin_char1.select <- clin.data.tcga.legacy.comm.all.data.select %>% dplyr::filter(Patient.ID %in% names(mrna.data)[-1])
plot_file_name = "mrna_patient_distribution_tum_stage.png"
plot_regression_status (clin_char1.select, plot_file_name)
### subset rnaseq data for the selected patients
mrna.data.select <- mrna.data %>% 
  dplyr::select(c("Hugo_Symbol", clin_char1.select$Patient.ID)) %>%
  na.exclude(.)
### the total number of duplicated genes is 17 hence we can exclude them.
mrna.data.select <- as.data.frame(mrna.data.select[!duplicated(mrna.data.select$Hugo_Symbol),])
### remove imputed expression values in the processed data
count <- apply(mrna.data.select[,-1], 1, function(x)length(unique(x)))
mrna.data.select <- mrna.data.select[count>1,]
write.table(mrna.data.select, file = "data_mrna_select_brca.txt", quote = FALSE, row.names = FALSE, sep = "\t")


clin_char1.select <- clin.data.tcga.legacy.comm.all.data.select %>% dplyr::filter(Patient.ID %in% names(mirna.data)[-1])
plot_file_name = "mirna_patient_distribution_tum_stage.png"
plot_regression_status (clin_char1.select, plot_file_name)
### subset mirna data for the selected patients
mirna.data.select <- mirna.data %>% 
  dplyr::select(c("Hugo_Symbol", clin_char1.select$Patient.ID)) %>%
  na.exclude(.)
### the total number of duplicated genes is 17 hence we can exclude them.
mirna.data.select <- as.data.frame(mirna.data.select[!duplicated(mirna.data.select$Hugo_Symbol),])
### remove imputed expression values in the processed data
count <- apply(mirna.data.select[,-1], 1, function(x)length(unique(x)))
mirna.data.select <- mirna.data.select[count>1,]
write.table(mirna.data.select, file = "data_mirna_select_brca.txt", quote = FALSE, row.names = FALSE, sep = "\t")


clin_char1.select <- clin.data.tcga.legacy.comm.all.data.select %>% dplyr::filter(Patient.ID %in% names(cna.data)[-1])
plot_file_name = "cna_patient_distribution_tum_stage.png"
plot_regression_status (clin_char1.select, plot_file_name)
### subset cna data for the selected patients
cna.data.select <- cna.data %>% 
  dplyr::select(c("Hugo_Symbol", clin_char1.select$Patient.ID)) %>%
  na.exclude(.)
### the total number of duplicated genes is 17 hence we can exclude them.
cna.data.select <- as.data.frame(cna.data.select[!duplicated(cna.data.select$Hugo_Symbol),])
### remove imputed expression values in the processed data
count <- apply(cna.data.select[,-1], 1, function(x)length(unique(x)))
cna.data.select <- cna.data.select[count>1,]
write.table(cna.data.select, file = "data_cna_select_brca.txt", quote = FALSE, row.names = FALSE, sep = "\t")


clin_char1.select <- clin.data.tcga.legacy.comm.all.data.select %>% dplyr::filter(Patient.ID %in% names(mut.data)[-1])
plot_file_name = "mut_patient_distribution_tum_stage.png"
plot_regression_status (clin_char1.select, plot_file_name)
### subset mutation data for the selected patients
mut.data.select <- mut.data %>% 
  dplyr::select(c("Hugo_Symbol", clin_char1.select$Patient.ID)) %>%
  na.exclude(.)
### the total number of duplicated genes is 17 hence we can exclude them.
mut.data.select <- as.data.frame(mut.data.select[!duplicated(mut.data.select$Hugo_Symbol),])
### remove imputed expression values in the processed data
count <- apply(mut.data.select[,-1], 1, function(x)length(unique(x)))
mut.data.select <- mut.data.select[count>1,]
write.table(mut.data.select, file = "data_mut_select_brca.txt", quote = FALSE, row.names = FALSE, sep = "\t")


clin_char1.select <- clin.data.tcga.legacy.comm.all.data.select %>% dplyr::filter(Patient.ID %in% common_patients_omics_data)
plot_file_name = "common_patient_distribution_tum_stage.png"
plot_regression_status (clin_char1.select, plot_file_name)

### check the number of genes we are loosing when excluding na values
### Not much difference in this data
dim(mrna.data) - dim(na.exclude(mrna.data))
dim(mirna.data) - dim(na.exclude(mirna.data))
dim(cna.data) - dim(na.exclude(cna.data))
dim(mut.data) - dim(na.exclude(mut.data))

### check the number of subsetted data
dim(mrna.data.select)
dim(mirna.data.select)
dim(cna.data.select)
dim(mut.data.select)

### make the gene id as row names and exclude from the dataframe
row.names(mrna.data.select) <- mrna.data.select$Hugo_Symbol
row.names(mirna.data.select) <- mirna.data.select$Hugo_Symbol
row.names(cna.data.select) <- cna.data.select$Hugo_Symbol
row.names(mut.data.select) <- mut.data.select$Hugo_Symbol

## Exclude the redundant hugo_symbol (gene symbol) from the dataframe 
mrna.data.select <- mrna.data.select[,-1]
mirna.data.select <- mirna.data.select[,-1]
cna.data.select <- cna.data.select[,-1]
mut.data.select <- mut.data.select[,-1]

### transpose the dataframe for exploratory analysis
mrna.data.select.t <- t(mrna.data.select)
mirna.data.select.t <- t(mirna.data.select)
cna.data.select.t <- t(cna.data.select)
mut.data.select.t <- t(mut.data.select)

########### Principal Component Analysis
### run pca for the data
mrna.pca <- prcomp(mrna.data.select.t, scale. = TRUE, center = TRUE)
mirna.pca <- prcomp(mirna.data.select.t, scale. = TRUE, center = TRUE)
cna.pca <- prcomp(cna.data.select.t, scale. = TRUE, center = TRUE)
mut.pca <- prcomp(mut.data.select.t, scale. = TRUE, center = TRUE)

### plot the results, color for ER status and shape for recurrent status
p1 <- autoplot(mrna.pca, main="PCA for mrna data", 
               data=clin.data.tcga.legacy.comm.all.data, 
               shape ="Disease.Free.Status", colour = "ER.Status.By.IHC", size = 3.5) + theme_light(base_size = 14)

p2 <- autoplot(mirna.pca, main="PCA for mirna data", 
               data=clin.data.tcga.legacy.comm.all.data %>% dplyr::filter(Patient.ID %in% names(mirna.data)[-1]), 
               shape ="Disease.Free.Status", colour = "ER.Status.By.IHC", size = 3.5) + theme_light(base_size = 14)

p3 <- autoplot(cna.pca, main="PCA for cna data", 
               data=clin.data.tcga.legacy.comm.all.data, 
               shape ="Disease.Free.Status", colour = "ER.Status.By.IHC", size = 3.5) + theme_light(base_size = 14)

p4 <- autoplot(mut.pca, main="PCA for mutation data", 
               data=clin.data.tcga.legacy.comm.all.data, 
               shape ="Disease.Free.Status", colour = "ER.Status.By.IHC", size = 3.5) + theme_light(base_size = 14)

png("omics_data_pca.png", width = 8, height = 24, units = "in", res = 300)
p <- plot_grid(p1, p2, p3, p4, ncol=1)
plot_grid(p, ncol=1) # rel_heights values control title margins
dev.off()

