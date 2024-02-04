### 0. OVERVIEW ####
#AIMS:
# INTEGRATION FOR 2nd ROUND OF SUBCLUSTERING for non-CD4+ and CD4+T cells
# this code will filter the gene and perform RPCA integration

# SECTIONS:
# 1. SET UP: load the library,integration function, and Seurat object(s), and set the directory, parameters and variables
# 2. INTEGRATION: filter the gene (and re-normalize the dataset) and perform integration (RPCA)

# NOTES: 
# before running this code please:
# 1. check and change parameters and variables accordingly (in section 1)
# 2. write the input directory for the Seurat object (dir) and the output directory (for the results)

### 1. SET UP ####
##### load the libraries #####
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)

##### set the directory #####
dir = "/mnt/sod2/csb6/AIDA/AIDA_DataFreeze_v2"
plotdir = "/mnt/sod2/csb6/eliora/AIDA_data_freeze_v2/20230827_T_NK_subclustering/CD4_subclustering" #output directory
dir.create(plotdir)
setwd(plotdir)

##### set up parameter and variable ##### 
batch_varibale = "Library"     #the batch will be based on this column name
batch_check = "Country"        #plot UMAP split and group based on this column name - to check the integration result
name_cell_type = ""            #please insert: "T_NK_CD4" or "T_NK_nonCD4"

##### load integration function #####
#RPCA is used as the method of integration

integration.RPCA <- function(file.to.integrate, batch.labels, dir.to.save, cell.name, K.parameter=NULL, nPC.parameter=NULL){
  #set the parameter and the input object
  if(is.null(K.parameter)){K.param = 100} else {K.param = K.parameter}
  print(paste0("K.parameter is ",K.param))
  finalSample <- file.to.integrate
  finalSample@active.assay <- "RNA" #setting the default assay as RNA
  
  #normalize and find the variable features per batch 
  Seurat.List <- lapply(X = SplitObject(finalSample, split.by = batch.labels), FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  features <- SelectIntegrationFeatures(object.list = Seurat.List)
  Seurat.List <- lapply(X = Seurat.List, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE,npcs = 50)
  })
  x <- table(finalSample[[batch.labels]])
  reference_dataset <- which(names(Seurat.List) == names(which(x == max(x))))
  print("Done with feature selection")
  
  #finding anchors 
  if(is.null(nPC.parameter)){nPC.parameter = 30} else {nPC.parameter = nPC.parameter}
  anchors <- FindIntegrationAnchors(object.list =Seurat.List, anchor.features = features, 
                                    reduction = "rpca", reference = reference_dataset, dims = 1:nPC.parameter)
  print("Done with anchor selection")
  
  #integration
  Combined <- IntegrateData(anchorset = anchors, k.weight = K.param, dims = 1:nPC.parameter)
  print("Done with Integration, Now saving RDS")
  saveRDS(Combined, file.path(dir.to.save,paste0(cell.name,"_integrated_RPCA.rds")))
  print("RPCA done")
  
  return(Combined)
} 

integration.RPCA.UMAP <- function(file.to.integrate, dir.to.save, batch.labels, cell.name, nPC.parameter=NULL){
  #set the parameter and the input object
  if(is.null(nPC.parameter)){nPC.parameter = 30} else {nPC.parameter = nPC.parameter}
  finalSample <- file.to.integrate
  finalSample@active.assay <- 'integrated'
  finalSample <- finalSample %>% ScaleData() %>% RunPCA() %>% RunUMAP(reduction = "pca", dims = 1:nPC.parameter)
  
  #plot UMAP
  print("plotting UMAP...")
  number_unique = nrow(unique(finalSample[[batch.labels]]))
  png(file.path(dir.to.save, paste0(cell.name,"_integrated_RPCA_group_by_",batch.labels,".png")), type = "cairo", height=500, width= 600)
  print(DimPlot(finalSample, group.by = batch.labels, reduction = 'umap')+ labs(title = paste0(cell.name,"_integrated_RPCA_group_by_",batch.labels)))
  dev.off()
  png(file.path(dir.to.save, paste0(cell.name,"_integrated_RPCA_split_by_",batch.labels,".png")), type = "cairo", height=500, width= 500*number_unique)
  print(DimPlot(finalSample, group.by = batch.labels, split.by = batch.labels, reduction = 'umap')+ labs(title = paste0(cell.name,"_integrated_RPCA_split_by_",batch.labels)))
  dev.off()
  
  return(finalSample)
}

##### load the rds #####
#we use the original seurat object to maintain the original list of genes 
print("loading the seurat object ...")
cell_list = readRDS("") #the directory of the list of cells 
T_seurat = readRDS(file.path(dir,"AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T.RDS"))
NK_seurat = readRDS(file.path(dir,"AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK.RDS"))

##subset seurat objects and merge them together 
T_seurat = subset(T_seurat, cells = intersect(Cells(T_seurat),cell_list))
NK_seurat = subset(NK_seurat, cells = intersect(Cells(NK_seurat),cell_list))
seurat_subset <- merge(T_seurat,NK_seurat)
gene_cell_number =dim(seurat_subset)

### 2. INTEGRATION ####
#### gene filter ####
#filter the gene: subset the seurat object to include genes expressed in >= 0.1% of cells 
#gene filter should be performed before integrating the Seurat object 
vec_cells_per_gene <- rowSums(seurat_subset@assays$RNA@counts > 0)
vec_genes_to_keep <- names(vec_cells_per_gene)[which(vec_cells_per_gene >= dim(seurat_subset)[2]*0.001)]
seurat_subset <- seurat_subset[vec_genes_to_keep, ]

#re-normalize
seurat_subset <- NormalizeData(seurat_subset, normalization.method = "LogNormalize", scale.factor = 10000)

#compile the number of gene and cells
gene_cell_number = rbind(gene_cell_number,dim(seurat_subset))
rownames(gene_cell_number) = c("ori","after_gene_filter")
colnames(gene_cell_number) = c("number_gene","number_cells")
write.csv(gene_cell_number, paste0(name_cell_type,"_number_of_cells_gene.csv"))

#### integrate the data - RPCA #### 
#run the integration function
rds.integrated <- integration.RPCA(seurat_subset, batch_varibale, plotdir, name_cell_type)
#plot UMAP: checking batch effect
rds.integrated = integration.RPCA.UMAP(rds.integrated, plotdir, batch_check, name_cell_type)
print("seurat object is filtered and integrated...")




