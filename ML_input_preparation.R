library(limma)
library(GEOquery)
library(org.Hs.eg.db)
library(dplyr)
library(ggrepel)
library(IOBR)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(dplyr)
library(tidyr)
library(enrichplot)


#####################

tcga<-readRDS('TCGA_stad_OS.rds')
gse59<-readRDS('gse15459_matrix_OS.rds')
gse53<-readRDS('gse26253_matrix_OS.rds')

### Get the common genes
common_genes <- Reduce(intersect, list(colnames(tcga),
                                       colnames(gse59), 
                                       colnames(gse53)))
saveRDS(common_genes,file = 'common_gene_3dataset.rds')
### subset the matrices

tcga<-tcga[,common_genes]
gse59<-gse59[,common_genes]
gse53<-gse53[,common_genes]

#### make a list of input dataframe
input<-list('tcga'=tcga,
            'gse53'=gse53,
            'gse59'=gse59)
           

metabolism<-unlist(signature_metabolism)


#saveRDS(input,file = 'input_ML_4dataframes.RDS')


### take a list of DEG of fibroblast (or other cell type) from single cell dataset
markers.list<-readxl::read_xlsx('caf_marker.xlsx')
markers.list<-as.data.frame(markers.list)
fibroblast_genes <- markers.list[markers.list$cluster == 'CAF_6', "gene"]

#list_overlap<-Reduce(intersect,list(common_genes,df_res_final$Gene,metabolism)) #### DEG with log2fc>1 TCGA dataset

list_overlap<-Reduce(intersect,list(common_genes,fibroblast_genes,metabolism)) #### DEG with log2fc>1 TCGA dataset
saveRDS(list_overlap,file = 'list_overalap_meta_caf.rds')
### after applying Cox-univariate (p<0.01)

list_overlap_uni<-c('time','status',selected_genes_stad_tcga) # refer to TCGA-STAD-new.R

#input<-readRDS('input_ML_5dataframes.RDS')

##### subset the input

subset_input <- lapply(input, function(df) {
  valid_columns <- intersect(names(df), list_overlap_uni)
  df[, valid_columns, drop = FALSE]  # Subset by column names
})
subset_input$tcga$time<-round(subset_input$tcga$time/30)

saveRDS(subset_input,file = 'subset_input_ML_dataframes_after_univariate.RDS')


#This is the input for 'Machine_gastric.R



