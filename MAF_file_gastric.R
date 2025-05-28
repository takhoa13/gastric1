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
library(GEOquery)
library(dplyr)
library(maftools)
library(TCGAmutations)

##########
### read input data, in this case 'risk _scores' object in 'forestplot_nomogram_clinical.R'
input_clinical<-risk_scores$tcga
high_risk_ids <- input_clinical$submittedID[input_clinical$Risk == 'high risk']
low_risk_ids <- input_clinical$submittedID[input_clinical$Risk == 'low risk']

  
maf<-TCGAmutations::tcga_load(study = "STAD",source = "Firehose")
maf_lowrisk<-subsetMaf(maf = maf, tsb = low_risk_ids, mafObj = T)
maf_highrisk<-subsetMaf(maf = maf, tsb = high_risk_ids, mafObj = T)
###plot oncoplot

oncoplot(maf = maf_lowrisk, top = 10,fontSize = 1,colors = vc_cols)


#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt <- mafCompare(m1 = maf_lowrisk, m2 = maf_highrisk, m1Name = 'Low risk', m2Name = 'High risk', minMut = 5)
print(pt.vs.rt)
pt.vs.rt$results<-pt.vs.rt$results[pt.vs.rt$results$or!='Inf',]
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.0001)
?forestPlot

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Spectral')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

coOncoplot(m1 = maf_lowrisk, m2 = maf_highrisk, m1Name = 'Low risk', m2Name = 'High risk', removeNonMutated = TRUE,colors = vc_cols)
coBarplot(m1 = maf_lowrisk, m2 = maf_highrisk, m1Name = 'Low risk', m2Name = 'High risk',legendTxtSize= 0.7,pctSize=0.7,geneSize = 0.7,colors = vc_cols,borderCol = "black")
somaticInteractions(maf = maf_highrisk, top = 10, pvalue = c(0.05, 0.01),showSum = F,showSigSymbols = T,leftMar = 8)
plotmafSummary(maf = maf_highrisk, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


#mafSurvival(maf = maf_lowrisk, genes = 'NEB', time = 'days_to_last_followup', Status = 'vital_status', isTCGA = TRUE)


##### make mutation matrix
mut_list <- make_mut_matrix(maf = maf, isTCGA   = T, category = "multi")


drugInteractions(maf = maf_highrisk, fontSize = 0.75,plotType = 'bar')


dev.off()
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
maf_lowrisk.tnm = trinucleotideMatrix(maf = maf_lowrisk, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
plotApobecDiff(tnm = maf_lowrisk.tnm, maf = maf_lowrisk, pVal = 0.2)

