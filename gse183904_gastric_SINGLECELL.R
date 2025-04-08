library(Matrix)
library(Seurat)
#library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(gridExtra)
#library(SCENIC)
library(celldex)
library(scRNAseq)
library(ggpubr)
library(hdf5r)
library(RColorBrewer)
library(devtools)
library(AnnoProbe)
#library(rjags)
#library(infercnv)
library(presto)
library(SCP)
library(reticulate)
library(forcats)
library(copykat)
library(harmony)
library(R.utils)



# Initialize the first Seurat object with the first sample
setwd('/home/lenovo/Documents/NCKH_Duy Tan/Gastric/')
first_sample_folder <- list.dirs(recursive = FALSE)[1]

  # Get the file paths for the count matrix and barcode files
# Path to the folder containing the CSV files
#sample_folder <- "path/to/your/folder"  # Update this with your folder path

# List all CSV files in the folder
count_matrix_file <- list.files(sample_folder, pattern = "\\.csv$", full.names = TRUE)

# Initialize the main Seurat object as NULL (to hold merged data later)
seurat_obj <- NULL

# Loop over each CSV file, create a Seurat object, and merge them
for (file in count_matrix_file) {
  # Read the count matrix file
  count_matrix <- read.csv(file, header = TRUE, row.names = 1)
  
  # Create a new Seurat object for the current sample
  sample_seurat <- CreateSeuratObject(counts = count_matrix)
  
  # Extract the base name of the file without extension and directory for unique cell IDs
  cell_id <- tools::file_path_sans_ext(basename(file))
  
  # Rename cells in the Seurat object with the file name (so each cell ID is unique)
  sample_seurat <- RenameCells(object = sample_seurat, add.cell.id = cell_id)
  
  # If this is the first file, initialize seurat_obj with sample_seurat
  if (is.null(seurat_obj)) {
    seurat_obj <- sample_seurat
  } else {
    # If seurat_obj is already initialized, merge the current sample into it
    seurat_obj <- merge(seurat_obj, y = sample_seurat)
  }
  
  # Optionally print the progress
  cat("Processed and merged:", cell_id, "\n")
}

# Check the merged Seurat object
print(seurat_obj)



seurat_obj$sample <- rownames(seurat_obj@meta.data)
# split sample column
seurat_obj@meta.data <- separate(seurat_obj@meta.data, col = 'sample', into = c('GSE', 'Sample','library','barcode'), 
                                sep = '_')

# calculate mitochondrial percentage
seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj, pattern='^MT-')


saveRDS(seurat_obj,file = 'gse183904_gastric_raw.rds')

###################3

#######
# read RDS seurat object again

seurat_obj<-readRDS('./gse183904_gastric_raw.rds')


# filtering
merged_seurat_filtered <- subset(seurat_obj, subset = nCount_RNA > 500 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 20)

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered,nfeatures=2000)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered,resolution=0.1)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- RunTSNE(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- JoinLayers(merged_seurat_filtered)

### perform Harmony
merged_seurat_filtered <- RunHarmony(merged_seurat_filtered, "patient")

# Create a new column named 'group' in the metadata
sample_group_mapping <- c(
  "GSM5573466" = "Normal",   # sample1
  "GSM5573467" = "Tumor",    # sample2
  "GSM5573468" = "Tumor",    # sample3
  "GSM5573469" = "Normal",   # sample4
  "GSM5573470" = "Tumor",    # sample5
  "GSM5573471" = "Normal",   # sample6
  "GSM5573472" = "Tumor",    # sample7
  "GSM5573473" = "Tumor",    # sample8
  "GSM5573474" = "Normal",   # sample9
  "GSM5573475" = "Tumor",    # sample10
  "GSM5573476" = "Normal",   # sample11
  "GSM5573477" = "Tumor",    # sample12
  "GSM5573478" = "Tumor",    # sample13
  "GSM5573479" = "Tumor",    # sample14
  "GSM5573480" = "Tumor",    # sample15
  "GSM5573481" = "Tumor",    # sample16
  "GSM5573482" = "Tumor",    # sample17
  "GSM5573483" = "Tumor",    # sample18
  "GSM5573484" = "Tumor-Peritonium",    # sample19 (Peritonium tissue)
  "GSM5573485" = "Tumor-Peritonium",    # sample20 (Peritonium tissue)
  "GSM5573486" = "Normal",   # sample21
  "GSM5573487" = "Tumor",    # sample22
  "GSM5573488" = "Normal",   # sample23
  "GSM5573489" = "Tumor",    # sample24
  "GSM5573490" = "Normal",   # sample25
  "GSM5573491" = "Tumor",    # sample26
  "GSM5573492" = "Tumor",    # sample27
  "GSM5573493" = "Tumor",    # sample28
  "GSM5573494" = "Tumor",    # sample29
  "GSM5573495" = "Tumor",    # sample30
  "GSM5573496" = "Normal",   # sample31
  "GSM5573497" = "Tumor",    # sample32
  "GSM5573498" = "Tumor",    # sample33
  "GSM5573499" = "Tumor",    # sample34
  "GSM5573500" = "Normal",   # sample35
  "GSM5573501" = "Tumor",    # sample36
  "GSM5573502" = "Normal-Peritonium",   # sample37 (Peritonium tissue)
  "GSM5573503" = "Tumor-Peritonium",    # sample38 (Peritonium tissue)
  "GSM5573504" = "Tumor",    # sample39
  "GSM5573505" = "Tumor"     # sample40
)

patient_mapping <- c(
  "sample1" = "NGCII518",   
  "sample2" = "NGCII518",    
  "sample3" = "NGCII519",    
  "sample4" = "NGCII520",   
  "sample5" = "NGCII520",    
  "sample6" = "NGCII521",   
  "sample7" = "NGCII521",    
  "sample8" = "NGCII524",   
  "sample9" = "NGCII525",   
  "sample10" = "NGCII525",    
  "sample11" = "NGCII527",   
  "sample12" = "NGCII527",    
  "sample13" = "NGCII529",    
  "sample14" = "NGCII502",    
  "sample15" = "NGCII511",   
  "sample16" = "NGCII499",   
  "sample17" = "NGCII509",    
  "sample18" = "NGCII498",    
  "sample19" = "NGCII011",    
  "sample20" = "NGCII015",    
  "sample21" = "NGCII513",  
  "sample22" = "NGCII513",    
  "sample23" = "NGCII514",   
  "sample24" = "NGCII514",   
  "sample25" = "NGCII512",  
  "sample26" = "NGCII512",    
  "sample27" = "NGCII510",    
  "sample28" = "NGCII531",    
  "sample29" = "NGCII522",   
  "sample30" = "NGCII533",    
  "sample31" = "NGCII536",   
  "sample32" = "NGCII539",   
  "sample33" = "NGCII540",    
  "sample34" = "NGCII541",    
  "sample35" = "NGCII538",   
  "sample36" = "NGCII538",    
  "sample37" = "NGCII538",   
  "sample38" = "NGCII538",    
  "sample39" = "NGCII543",    
  "sample40" = "NGCII545"     
)

gender_mapping <- c(
  "NGCII518"='M',   
  "NGCII519"='M',    
  "NGCII520"= 'M',    
  "NGCII521"='M',   
  "NGCII524"= 'M',   
  "NGCII525"='M',    
  "NGCII527"= 'M',    
  "NGCII529"='F',    
  "NGCII502"='M',    
  "NGCII511"='M',   
  "NGCII499"='M',   
  "NGCII509"='M',    
  "NGCII498"='F',    
  "NGCII011"='F',    
  "NGCII015"='M',    
  "NGCII513"='M',  
  "NGCII514"='F',   
  "NGCII512"='F',    
  "NGCII510"='F',    
  "NGCII531"='F',    
  "NGCII522"='M',   
  "NGCII533"='M',    
  "NGCII536"='M',   
  "NGCII539"='F',   
  "NGCII540"='M',    
  "NGCII541"='M',    
  "NGCII538"='F',   
  "NGCII543"='F',    
  "NGCII545"='M'     
)

location_mapping <- c(
  "NGCII518"='Antrum',   
  "NGCII519"='NA',    
  "NGCII520"= 'Distal',    
  "NGCII521"='Antrum',   
  "NGCII524"= 'NA',   
  "NGCII525"='Antrum',    
  "NGCII527"= 'NA',    
  "NGCII529"='Antrum',    
  "NGCII502"='Cardia',    
  "NGCII511"='Antrum',   
  "NGCII499"='Pylorus',   
  "NGCII509"='Cardia',    
  "NGCII498"='Antrum',    
  "NGCII011"='Antrum',    
  "NGCII015"='Peritoneum',    
  "NGCII513"='Antrum',  
  "NGCII514"='Body-antrum',   
  "NGCII512"='Antrum',    
  "NGCII510"='Antrum',    
  "NGCII531"='Proximal',    
  "NGCII522"='Body-antrum',   
  "NGCII533"='Antrum',    
  "NGCII536"='Antrum',   
  "NGCII539"='Antrum',   
  "NGCII540"='Antrum',    
  "NGCII541"='Antrum',    
  "NGCII538"='Antrum',   
  "NGCII543"='Antrum',    
  "NGCII545"='Antrum'     
)
# Assign the group based on the mapping
merged_seurat_filtered@meta.data$group <- sample_group_mapping[merged_seurat_filtered@meta.data$GSE]
merged_seurat_filtered@meta.data$patient <- patient_mapping[merged_seurat_filtered@meta.data$Sample]
merged_seurat_filtered@meta.data$gender <- gender_mapping[merged_seurat_filtered@meta.data$patient]
merged_seurat_filtered@meta.data$location <- location_mapping[merged_seurat_filtered@meta.data$patient]


saveRDS(merged_seurat_filtered,file = 'gse183904_gastric_filtered.rds')


#Read RDS object
merged_seurat_filtered<-readRDS('gse183904_gastric_filtered.rds')


###subset data to pairs according to clinical factors

metadata<-merged_seurat_filtered@meta.data
merged_seurat_filtered<-subset(merged_seurat_filtered, subset = (group%in%c('Tumor','Normal')))

# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'patient')
p1
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'group')

DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'group',label = F,
        label.size = 4, label.color = "black")

FeaturePlot(merged_seurat_filtered, reduction = 'umap', features = c('REG4','ACKR1'))


#### DEG
#run DEG for each cluster
Idents(merged_seurat_filtered) <- "seurat_clusters"
DefaultAssay(merged_seurat_filtered) <- "RNA"
markers.list = FindAllMarkers(merged_seurat_filtered,logfc.threshold = 0.5,test.use = "wilcox",only.pos = T)
writexl::write_xlsx(markers.list,path = 'DEG183904_new.xlsx')
cluster_annotations <- list(
  '0' = 'Lymphocyte',
  '1' = 'Plasma cell',
  '2' = 'Epithelial cell',
  '3' = 'Macrophage',
  '4' = 'Epithelial cell',
  '5' = 'Endothelial cell',
  '6' = 'Fibroblast', # ASC / Neuron / ODC markers all present,
  '7' = 'Fibroblast',
  '8' = 'B cell',
  '9' = 'Mast cell',
  '10' = 'Pericyte',
  '11' = 'Lymphocyte',
  '12' = 'Epithelial cell'
)

merged_seurat_filtered@meta.data$celltype <- unlist(cluster_annotations[merged_seurat_filtered$seurat_clusters])
merged_seurat_filtered$celltype<-factor(merged_seurat_filtered$celltype,levels = c('Lymphocyte','Plasma cell',
                                                                                   'Epithelial cell','Macrophage',
                                                                                   'Endothelial cell','Fibroblast',
                                                                                   'B cell','Mast cell','Pericyte'))


#############
#############

DimPlot(merged_seurat_filtered, reduction = "umap", label = T, pt.size = 0.5,label.color = "black",
        group.by ='celltype',cols = 'Set3',label.size =3 ,label.box = T,repel = T) 


##### visualize

library(dplyr)
merged_seurat_filtered <- RunDEtest(merged_seurat_filtered, group_by = "celltype", BPPARAM = BiocParallel::SerialParam())


de_filter <- filter(merged_seurat_filtered@tools$DEtest_celltype$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)

de_top <- de_filter %>%
  group_by(gene) %>%
  top_n(1, avg_log2FC) %>%
  group_by(group1) %>%
  top_n(3, avg_log2FC)

Idents(merged_seurat_filtered)<-merged_seurat_filtered$celltype


#### downgrade the Assay of seurat object
merged_seurat_filtered[["RNA"]] <- as(merged_seurat_filtered[["RNA"]], "Assay")

merged_seurat_filtered <- AnnotateFeatures(merged_seurat_filtered, species = "Homo_sapiens", db = c("TF", "CSPA"))

library(magick)

ht5 <- GroupHeatmap(merged_seurat_filtered,
                    features = de_top$gene,  group.by = "celltype",group_palette = "Set3",
                    heatmap_palette = "YlOrRd",cell_split_palette = "Set3",
                    cell_annotation = c("gender", "location"), cell_annotation_palette = c("Set3", "Paired"),
                    cell_annotation_params = list(width = unit(10, "mm")),
                   # feature_annotation = c("TF", "CSPA"),
                   # feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
                    add_dot = TRUE, add_bg = TRUE, add_reticle = TRUE,
                    flip = T, column_title_rot = 45, nlabel = 0, show_row_names = T,show_column_names = T
)
ht5

ht6 <- GroupHeatmap(merged_seurat_filtered,
                    features = de_top$gene, feature_split = de_top$gene, group.by = "celltype",group_palette = "Set3",
                    heatmap_palette = "YlOrRd",cell_split_palette = "Set3",
                    cell_annotation = c("gender", "location"), cell_annotation_palette = c("Set3", "Paired"),
                    cell_annotation_params = list(width = unit(20, "mm")),
                    # feature_annotation = c("TF", "CSPA"),
                    # feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
                    add_dot = TRUE, add_bg = TRUE, add_reticle = TRUE,
                    flip = F, column_title_rot = 45, nlabel = 0, show_row_names = F,show_column_names = T
)
ht6

merged_seurat_filtered <- RunEnrichment(
  srt = merged_seurat_filtered, group_by = "celltype", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05",BPPARAM = BiocParallel::SerialParam())

EnrichmentPlot(
  srt = merged_seurat_filtered, group_by = "celltype", group_use = "Fibroblast",
  plot_type = "network", network_layout	= 'kk',network_adjiter = 50,network_adjscale = 100,
  theme_use = "theme_blank",legend.position = "top",
  legend.direction = "horizontal"
)

CellStatPlot(merged_seurat_filtered, stat.by = "group", group.by = "celltype", 
             label = TRUE,palette = 'Set2',
             position = "dodge",theme_use = 'theme_scp') %>% panel_fix(height = 4, width = 10)

EnrichmentPlot(srt = merged_seurat_filtered, group_by = "celltype", 
               compare_only_sig = TRUE,plot_type = "comparison",split_by = c("Groups"))


merged_seurat_filtered <- RunGSEA(
  srt = merged_seurat_filtered, group_by = "celltype", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05",BPPARAM = BiocParallel::SerialParam()
)
#GSEAPlot(srt = merged_seurat_filtered, group_by = "CellType", group_use = "Endocrine", id_use = "GO:0007186")


library(BiocParallel)
library(org.Hs.eg.db)

register(BiocParallel::SerialParam())
DEGs <- merged_seurat_filtered@tools$DEtest_celltype$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]

ht <- FeatureHeatmap(
  srt = merged_seurat_filtered, group.by = "celltype", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Homo_sapiens", db = c("GO_BP"), anno_terms = TRUE, 
  #feature_annotation = c("TF", "CSPA"), 
  #feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 7, width = 4, group_palette = "Set3",slot = 'data'
)
print(ht$plot)



###### Addmodule score (using the positive genes list from machine_gastric.R)
list_input<-readRDS('list_overalap_meta_caf.rds')

merged_seurat_filtered <- AddModuleScore(object = merged_seurat_filtered, features =list_input, name = "meta_list")
### Mask all un-needed cell type
VlnPlot(merged_seurat_filtered, features = "meta_list1", pt.size = 0,
        group.by = "celltype", cols = RColorBrewer::brewer.pal(n = length(unique(merged_seurat_filtered$celltype)), "Set3")) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Add boxplots without outliers
  stat_kruskal_test()

FeaturePlot(object = merged_seurat_filtered, features ='meta_list1')#+
scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu")))



############
#### subset fibroblast
############
caf <-subset(merged_seurat_filtered, subset = (celltype =='Fibroblast'))

#caf <-subset(merged_seurat_filtered, subset = (celltype =='Fibroblast'&group=='Tumor'))

saveRDS(caf, file = 'caf_cluster_tumor.rds')

caf<-readRDS('caf_cluster.rds')

comparison = list(c('Tumor','Normal'))
caf <- AddModuleScore(object = caf, features =list_input, name = "meta_list")
### Mask all un-needed cell type
VlnPlot(caf, features = "meta_list1", pt.size = 0,
        group.by = "group",split.by = 'group', cols = RColorBrewer::brewer.pal(n = length(unique(merged_seurat_filtered$celltype)), "Set3")) + 
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Add boxplots without outliers
  stat_anova_test(comparisons = comparison)+ ylim(-1,3)+
  labs(title = "Fibroblast")

FeaturePlot(object = merged_seurat_filtered, features ='meta_list1')#+
scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu")))

caf<-subset(caf, subset = (group == 'Tumor'))
#caf<-readRDS('caf_cluster_tumor.rds')


### recluster again
caf <- NormalizeData(object = caf)
caf <- FindVariableFeatures(object = caf,nfeatures=2000)
caf <- ScaleData(object = caf)
caf <- RunPCA(object = caf)
caf <- FindNeighbors(object = caf, dims = 1:20)
caf <- FindClusters(object = caf,resolution=0.05)
caf <- RunUMAP(object = caf, dims = 1:20)
caf <- RunTSNE(object = caf, dims = 1:20)

### perform Harmony
caf <- RunHarmony(caf, "patient")


DimPlot(caf, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE,
        label.size = 4, label.color = "black")

cluster_mapping <- c(
  "0" = "CAF_1",   # sample1
  "1" = "CAF_2",    # sample2
  "2" = "CAF_3",    # sample3
  "3" = "CAF_4",   # sample4
  "4" = "CAF_5",    # sample5
  "5" = "CAF_6")
caf@meta.data$cluster <- cluster_mapping[caf@meta.data$seurat_clusters]
Idents(caf)<-caf$cluster

DimPlot(caf, reduction = 'umap', group.by = 'cluster',label = TRUE, cols = 'Set2',
        label.size = 4, label.color = "black")


##### check survival of each cluster in CAF
## first, we need to calculate DEG of each cluster
caf_marker<-FindAllMarkers(caf,logfc.threshold = 0.5,test.use = "wilcox",only.pos = F)

## take a list of 50 genes in each cluster
caf_marker.surv=caf_marker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

writexl::write_xlsx(caf_marker,path = 'caf_marker.xlsx')

#remove genes from signature gene lists that do not exist in expression matrix ##otherwise the code won't work
library(scrabble)
mat_tumor<-readRDS('TCGA_stad_OS_and_matrix.rds')
sample_tumor<-readRDS('tcga_tumor_info.rds')

Counts=as.data.frame(t(mat_tumor[,-c(1,2)])) 
geneset.all <- caf_marker.surv[ caf_marker.surv$gene %in% rownames(Counts),]
geneset <- geneset.all[!duplicated(geneset.all$gene),]
genesetlist=geneset.all%>% split(x =.$gene, f = .$cluster,drop = F)
Scounts=Counts[geneset$gene,]
colCenter = function(m, by = 'mean') {
  m = as.matrix(m)
  if (by == 'mean')  by = T
  else if (by == 'median') by = matrixStats::colMedians(m)
  else stopifnot(is.numeric(by) & length(by) == ncol(m))
  scale(m, center = by, scale = F)
}
sc=score(mat=Scounts,
         groups=genesetlist,
         binmat = NULL,
         bins = NULL,
         controls = NULL,
         bin.control = F,
         center = T,
         nbin = 30,
         n = 100,
         replace = T)

library(survival)
library(survminer)

ObjName= "HumanStad-AllSamples-"
sc<-as.data.frame(sc)

TCGA.sc=cbind(mat_tumor[,c(1,2)],sc)


#for loop to assign "positive" or "negative" to cell scores then draw survival plots
data.subset=colnames(sc)
data.subset=sort(data.subset, decreasing = FALSE)
TCGA.sc2 =TCGA.sc
surv_object <- Surv(time = TCGA.sc$time, event = TCGA.sc$status)

Final <- list()


for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  TCGA.sc2 <- TCGA.sc2%>% mutate(Expression.Level = ifelse(TCGA.sc2[YY]>=0, "Positive", "Negative"))
  TCGA.sc2$Expression.Level <-factor(TCGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = TCGA.sc2)

  XX <- ggsurvplot(fit1, data = TCGA.sc2, pval = TRUE,
                   pval.cExpression.Levelrd = c(50, 1),
                   pval.size=6,legend="top",
                   risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, 
                   palette = c("green4","Red"),
                   font.x = c(14, "bold"),
                   font.y = c(14, "bold"),
                   legend.labs = c("Not-Enriched", "Enriched"),
                   font.tickslab = c(12, "plain"), 
                   xlab = "Time")
  Final[[i]] = XX
}
pdf(file ="survival plots TCGA allFigure.pdf", height = 8, width =10)

arrange_ggsurvplots(Final,print = TRUE, ncol = 3,nrow = 2,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()

caf <- RunSlingshot(srt = caf, group.by = "cluster", reduction = "umap")

CellDimPlot(caf, group.by = "cluster", reduction = "umap", lineages = paste0("Lineage", 1), lineages_span = 0.6)

caf <- RunDEtest(caf, group_by = "cluster", BPPARAM = BiocParallel::SerialParam())


de_filter <- filter(caf@tools$DEtest_cluster$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)

de_top <- de_filter %>%
  group_by(gene) %>%
  top_n(1, avg_log2FC) %>%
  group_by(group1) %>%
  top_n(3, avg_log2FC)

Idents(caf)<-caf$cluster

caf <- RunEnrichment(
  srt = caf, group_by = "cluster",  db = "GO_BP", #db = "MSigDB_H", 
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05")

EnrichmentPlot(srt = caf, group_by = "cluster", db = "GO_BP",#db = "MSigDB_H", 
               plot_type = "comparison")%>% panel_fix(height = 8, width = 2.5)
