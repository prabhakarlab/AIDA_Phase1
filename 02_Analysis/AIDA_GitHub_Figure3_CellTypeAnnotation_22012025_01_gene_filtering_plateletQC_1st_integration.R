### 0. OVERVIEW ####
## AIMS: 
# this code will filter the gene, remove the flagged_platelet_doublets, 
# and perform integration for the 1st round of sub-clustering (please refer another script for 2nd round of sub-clustering)

## SECTIONS: 
# 1. SET UP: load the library, seurat object(s), and integration function, and set the directory, parameters and variables
# 2. FILTER GENE: include genes expressed in >= 0.1% of cells 
# 3. PLATELET QC: flagged and removed the suspected platelet doublets based on the genes' expression
# 4. INTEGRATION: filter the gene (and re-normalize the dataset) and perform integration (RPCA)

## NOTES: 
# before running this code please:
# 1. check and change the parameter and variable accordingly (in section 1)
# 2. write the directory of the input of the seurat object (dir) and the output 

### 1. SET UP ####
##### load the libraries #####
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)

##### set the directory #####
dir = "/mnt/sod2/csb6/AIDA/AIDA_DataFreeze_v2" #input directory
plotdir = "" #output directory
dir.create(plotdir)
setwd(plotdir)

##### set up parameter and variable ##### 
#for platelet QC
marker_genes = c("ITGA2B","PF4","TUBB1","PPBP")
threshold = 0.3
check.labels = 'Country'       #the metadata column for country information 

#for integration
batch_varibale = "Library"     #the batch will be based on this column name
batch_check = "Country"        #plot UMAP split and group based on this column name - to check the integration result
name_cell_type = ""            #please insert: "B" , "T_NK" , or "pDC_Myeloid"

##### load the rds #####
print("loading the seurat object ...")
seurat_all = readRDS(file.path(dir,""))
colnames(seurat_all[[]])      #check the metadata

  #for lineages with 2 seurat objects
  rds_1 = readRDS(file.path(dir,""))
  rds_2 = readRDS(file.path(dir,""))
  seurat_all <- merge(rds_1,rds_2)

##### load integration function #####
#RPCA is used as the method if integration

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
  
### 2. FITLER GENE ####
print("seurat object is loaded, fitlering the gene now...")
##### filter the gene #####
#make a list of genes expressed in >= 0.1% of cells
vec_cells_per_gene <- rowSums(seurat_all@assays$RNA@counts > 0)
vec_genes_to_keep <- names(vec_cells_per_gene)[which(vec_cells_per_gene >= dim(seurat_all)[2]*0.001)]

  #for lineages with 2 seurat objects
  vec_cells_per_gene_1 <- rowSums(rds_1@assays$RNA@counts > 0)
  vec_genes_to_keep_1 <- names(vec_cells_per_gene_1)[which(vec_cells_per_gene_1 >= dim(rds_1)[2]*0.001)]
  vec_cells_per_gene_2 <- rowSums(rds_2@assays$RNA@counts > 0)
  vec_genes_to_keep_2 <- names(vec_cells_per_gene_2)[which(vec_cells_per_gene_2 >= dim(rds_2)[2]*0.001)]
  vec_genes_to_keep = union(vec_genes_to_keep_1,vec_genes_to_keep_2)
  
#subset the seurat object to contain only filtered genes
seurat_subset <- seurat_all[vec_genes_to_keep, ]
gene_cell_number = rbind(dim(seurat_all),dim(seurat_subset))

##### re-normalize #####
DefaultAssay(seurat_subset)<-'RNA'
seurat_subset <- NormalizeData(seurat_subset, normalization.method = "LogNormalize", scale.factor = 10000)

### 3. PLATELET QC ####
#flagging cells with platelet gene expression
#the criteria: top 30% based on the sum of genes' expression (rank only the non-0 total UMI)
print("seurat object filtered and re-normalized, platelet QC is starting...")

#### run UMAP - just for visualization #####
PBMC_r <-seurat_subset
PBMC_r <-FindVariableFeatures(PBMC_r, selection.method = "vst", nfeatures = 2000) 
PBMC_r <-ScaleData(PBMC_r) %>% RunPCA() %>% RunUMAP(reduction = "pca", dims = 1:30) #scaling is done for the top variable features

#### checking the expression of platelet's marker genes #####
#create feature plots 
DefaultAssay(PBMC_r) <- 'RNA'
png(file = file.path(plotdir,paste0(name_cell_type,"_UMAP_of_seurat_marker_genes_QC.png")),width=600, height=600, type = "cairo")
FeaturePlot(PBMC_r, slot = 'data', features = unlist(marker_genes), ncol = 2) 
dev.off()

#check the number of expressing cells and average expression
marker.df <- data.frame()
for (i in (1:length(marker_genes))){
  gene_ij <- marker_genes[i]
  if (gene_ij %in% rownames(PBMC_r[["RNA"]]@data)) {
    marker.ij <- GetAssayData(object = PBMC_r, assay = "RNA", slot = "data")[marker_genes[i],]
    number.ij <- sum(marker.ij>0)
    exp_per.ij <- round(sum(marker.ij>0)*100/length(marker.ij),4)
    exp_mean.ij.all<- round(mean(marker.ij),4)
    exp_mean.ij<- round(sum(marker.ij)/sum(marker.ij>0),4)
    cut.off.val <- min(head(sort(marker.ij,decreasing=TRUE),n=round(threshold*number.ij,0)))
  }else {
    exp_per.ij  = 0
    exp_mean.ij = 0
    exp_mean.ij.all = 0
    cut.off.val=0}
  marker.df.tmp <- data.frame(gene=gene_ij, 
                              num_expressed_cells = number.ij, 
                              percentage_expressed_cells=exp_per.ij, 
                              mean_expression_expressed_only = exp_mean.ij, 
                              mean_expression_all = exp_mean.ij.all,
                              cut_off_expression = cut.off.val)
  marker.df <- rbind(marker.df,marker.df.tmp)
}
marker.df$percentage_expressed_cells <- paste0(round(marker.df$percentage_expressed_cells,2),"%")
write.csv(marker.df,file.path(plotdir,paste0(name_cell_type,"_Expression_percentage_mean_QC.csv")))
rm(gene_ij,exp_per.ij, exp_mean.ij.all, exp_mean.ij, marker.df.tmp, marker.ij)

#### flagging the platelet #####
exp.df = GetAssayData(object = PBMC_r, assay = "RNA", slot = "data")[marker_genes,]
sum.per.cell = apply(exp.df,2, FUN=sum)
number.of.exp.cells = sum(sum.per.cell>0)
number.cells = round(threshold*number.of.exp.cells,0)
cut.off = min(head(sort(sum.per.cell,decreasing=TRUE),n=number.cells))
barcodes.removed.sum <- names(head(sort(sum.per.cell,decreasing=TRUE),n=number.cells))
saveRDS(barcodes.removed.sum,file.path(plotdir,paste0(name_cell_type,"_Barcode_list_of_platelet_flagged_cells.rds")))
save(number.cells,cut.off,file= paste0(name_cell_type,"_cut_off_expresion_and_number_of_cells.rda"))

#### check the resource of the flagged cells ####
subset.flagged = subset(seurat_subset, cells = barcodes.removed.sum)
write.csv(data.frame(table(subset.flagged[[check.labels]])), paste0(name_cell_type,"_",check.labels,"_removed_flagged_cells.csv"))

#### create the UMAP highlighting flagged cells ####
png(paste0(name_cell_type,"_UMAP_flagged_platelet.png"), type = "cairo",height=500, width= 500)
DimPlot(PBMC_r,group.by = NULL,pt.size = 0.5,cells.highlight = barcodes.removed.sum,
        cols.highlight = "#DE2D26") + NoLegend() + labs(title = paste0(name_cell_type,"_platelet_flagged_cells"))
dev.off()
print("UMAP is created, list of flagged cells is saved, removing falgged cells now...")

#### removing the flagged cells #####
seurat_subset <- subset(seurat_subset, cells = barcodes.removed.sum, invert = TRUE)
saveRDS(seurat_subset,paste0(name_cell_type,"_seurat_object_filtered_normalized_Platelet_QC.rds"))
rm(PBMC_r)
print("flagged cells are removed, proceed to integatring object...")

### 3. INTEGRATION ####
#### gene filter ####
#filter the gene: subset these datasets to include genes expressed in >= 0.1% of cells 
#gene filter should  be performed before integrating the Seurat object 
vec_cells_per_gene <- rowSums(seurat_subset@assays$RNA@counts > 0)
vec_genes_to_keep <- names(vec_cells_per_gene)[which(vec_cells_per_gene >= dim(seurat_subset)[2]*0.001)]
seurat_subset <- seurat_subset[vec_genes_to_keep, ]

#re-normalize
seurat_subset <- NormalizeData(seurat_subset, normalization.method = "LogNormalize", scale.factor = 10000)

#compile the number of gene and cells
gene_cell_number = rbind(gene_cell_number,dim(seurat_subset))
rownames(gene_cell_number) = c("ori","after_gene_filter_1","after_plateletQC_gene_filter2")
colnames(gene_cell_number) = c("number_gene","number_cells")
write.csv(gene_cell_number, paste0(name_cell_type,"_number_of_cells_gene.csv"))

#### integrate the data - RPCA #### 
#run the integration function
rds.integrated <- integration.RPCA(seurat_subset, batch_varibale, plotdir, name_cell_type)
#plot UMAP: checking batch effect
rds.integrated = integration.RPCA.UMAP(rds.integrated, plotdir, batch_check, name_cell_type)
print("seurat object is filtered and integrated...")


