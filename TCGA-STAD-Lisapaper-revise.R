library(TCGAWorkflowData)
# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
#library("genefilter")
library("DESeq2")
library(pheatmap)
library(colorRamp2)
library(ggrepel)
library(IOBR)
library(maftools)
##QUERY DATA FIRST


query_TCGA = GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq")
coad_res<-getResults(query_TCGA)
head(coad_res$sample_type)


summary(factor(coad_res$sample_type)) # summary of distinct tissues types present in this study
##QUERY DATA SELECTING SAMPLE TYPE
query_TCGA = GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  data.type = 'Gene Expression Quantification',
  workflow.type = "STAR - Counts",
  sample.type = c('Primary Tumor','Solid Tissue Normal'))


##DOWNLOAD TCGA DATA
GDCdownload(query = query_TCGA)

tcga_data = GDCprepare(query_TCGA,summarizedExperiment = T)

tcga_data@colData$deceased <- ifelse(tcga_data@colData$vital_status == "Alive", 0, 1)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
tcga_data@colData$overall_survival <- ifelse(tcga_data@colData$vital_status == "Alive",
                                             tcga_data@colData$days_to_last_follow_up,
                                             tcga_data@colData$days_to_death)

#check the size of the object
dim(tcga_data)
# Define the conditions to filter by
colnames(colData(tcga_data))

# Replace NAs with FALSE in the logical condition
# subset_condition <- with(tcga_data@colData, 
#                          ifelse(is.na(shortLetterCode), FALSE, 
#                                 shortLetterCode == 'TP' & age_at_index < 50) | 
#                            shortLetterCode == 'NT')
stad <- subset(tcga_data, select = colData(tcga_data)$shortLetterCode == "TP"  | colData(tcga_data)$shortLetterCode == "NT")

rownames(stad)<-(stad@rowRanges@elementMetadata$gene_name)

matrix <- assay(stad,'tpm_unstrand')
matrix <- matrix[rowSums(matrix) != 0,]
matrix <- log2(matrix + 1)




sample<-data.frame(ID=row.names(stad@colData),
                   age= stad@colData$age_at_index,subtype=stad@colData$paper_Molecular.Subtype,
                   KRAS =stad@colData$paper_KRAS.mutation,
                   MSI =stad@colData$paper_MSI.status,
                   CIMP =stad@colData$paper_CIMP.Category,
                   TP53 =stad@colData$paper_TP53.mutation,
                   submittedID=stad@colData$submitter_id,
                   group=stad@colData$sample_type,
                   sex=stad@colData$gender,
                   stage = stad@colData$ajcc_pathologic_stage,
                   T_stage=stad@colData$ajcc_pathologic_t,
                   M_stage = stad@colData$ajcc_pathologic_m,
                   N_stage=stad@colData$ajcc_pathologic_n, row.names = row.names(stad@colData))
sample$stage <- gsub("[ABC]", "", sample$stage)
sample$T_stage <- gsub("[abc]", "", sample$T_stage)
sample$M_stage <- gsub("[abc]", "", sample$M_stage)
sample$N_stage <- gsub("[abc]", "", sample$N_stage)


##perform DEG analysis
ddsSE <- DESeqDataSet(stad, design = ~ sample_type)


keep <- rowSums(counts(ddsSE)) >= 10
ddsSE <- ddsSE[keep,]
ddsSE$sample_type <-relevel(ddsSE$sample_type,ref = 'Solid Tissue Normal')

ddsSE <- DESeq(ddsSE,test = 'Wald')

resultsNames(ddsSE)


res <- results(ddsSE)
res1<-results(ddsSE,name = 'sample_type_Primary.Tumor_vs_Solid.Tissue.Normal')
df_res1<-as.data.frame(res1)
df_res1$Gene <-rownames(df_res1)
dim(res)
df_res1<-na.omit(df_res1)


#### take a list òf DEG
df_res_final<-df_res1[abs(df_res1$log2FoldChange)>1 & df_res1$padj<0.05,]
#writexl::write_xlsx(df_res_final,path = 'TCGA_stad_DEG.xlsx')

df_res_final1<-readxl::read_xlsx('TCGA_stad_DEG.xlsx')
###read DEG result


### cox-analysis of target gene
patients <- sample[sample$group=='Primary Tumor',]
stad_tumor <- subset(stad, select = colData(stad)$shortLetterCode == "TP" )

#selected_genes<- Reduce(intersect,list(common_genes,fibroblast_genes,metabolism))
status <-stad_tumor$deceased


time <-stad_tumor$overall_survival

#time <- ifelse(status == "1", time, max(time, na.rm = TRUE))

mat_tumor <- assay(stad_tumor)
mat_tumor <- mat_tumor[rowSums(mat_tumor) != 0,]
mat_tumor <- log2(mat_tumor + 1)

mat_tumor_saved_OS<-data.frame(time=time, status=status,t(mat_tumor))
sample_tumor<-sample[colnames(mat_tumor),]

saveRDS(sample_tumor,file = 'tcga_tumor_info.rds')  ## this is for 'clinical analysis'
saveRDS(mat_tumor,file = 'tcga_stad_mat_tumor.rds')
saveRDS(mat_tumor_saved_OS,file = 'TCGA_stad_OS_and_matrix.rds')
##
##### for saving time, we read the exp matrix that we saved before


mat_tumor<- readRDS('tcga_stad_mat_tumor.rds')
mat_tumor_saved_OS<-readRDS('TCGA_stad_OS_and_matrix.rds')

list_overlap<-readRDS('list_overalap_meta_caf.rds')
list36<-c("GPX3","ADH1B","ACSM5","MSRB3","PPP1R3B",
          "GXYLT2","NPR1","CPT1C","GPX8","PDE9A","PLOD2","GUCY1A2",
          "P4HA3",'DPYS','AOX1','RBP4','CH25H','CYP1B1','NT5E','GFPT2',
          'GAMT','ALDH1A2','PDE1B','LRAT','GPX7','PDE1A','ASPA', 'ACSS3','TUSC3','GALNT16',
          'GGT5','ST3GAL6','PRDM6','ENTPD2','VIM','PDK4')
list_overlap<-unique(append(list_overlap,list36))
stad_matrix <- as.data.frame(t(mat_tumor[list_overlap,]))## list overlap taken from the intersection between metabolism and DEG of our cell type

status <-mat_tumor_saved_OS$status


time <-mat_tumor_saved_OS$time##
stad_OS<-data.frame(stad_matrix,time=time, status=status)



library(survival)
library(glmnet)
# perform univariate cox regreGBMion for each gene
# Fit univariate Cox models for each column in 'stad_matrix' and extract relevant statistics
cox_uni_stad <- lapply(1:ncol(stad_matrix), function(i) {
  fit <- coxph(Surv(time, status) ~ stad_matrix[, i])  # Cox model for each variable
  
  # Extract coefficients and statistics
  coef <- summary(fit)$coefficients[1, "coef"]
  HR <- summary(fit)$coefficients[1, "exp(coef)"]
  se <- summary(fit)$coefficients[1, "se(coef)"]
  z <- summary(fit)$coefficients[1, "z"]
  pval <- summary(fit)$coefficients[1, "Pr(>|z|)"]
  
  # Calculate 95% CI bounds for the HR
  CI_low <- exp(coef - 1.96 * se)
  CI_up <- exp(coef + 1.96 * se)
  
  # Return the results as a named vector
  c(coef = coef, HR = HR, pval = pval, CI_low = CI_low, CI_up = CI_up)
})

# Combine results into a data frame
cox_uni_stad <- do.call(rbind, cox_uni_stad)  # Bind rows into a data frame
colnames(cox_uni_stad) <- c("coef", "HR", "pvalue", "95%CI_low", "95%CI_up")  # Set column names
rownames(cox_uni_stad) <- colnames(stad_matrix)  # Set row names based on matrix column names

# Convert to a data frame
cox_uni_stad <- as.data.frame(cox_uni_stad)


# select genes with p-value < 0.05
selected_genes_stad_tcga <- row.names(cox_uni_stad)[cox_uni_stad$pval < 0.01]
cox_uni_stad<-cox_uni_stad[selected_genes_stad_tcga,]
cox_uni_stad<-rownames_to_column(cox_uni_stad,var = 'ID')
final_stad_matrix<-stad_OS[, selected_genes_stad_tcga]
write.csv(cox_uni_stad,file = 'stad_uni_cox_full_36.csv')


####forest plot for these significant genes
p1<- sig_forest(cox_uni_stad, signature = "ID", 
                pvalue = 'pvalue',
                HR='HR',
                CI_low_0.95 = "95%CI_low",
                CI_up_0.95 = "95%CI_up",
                color_option =2,
                n = 70)



