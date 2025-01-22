### AIDA Phase 1 Data Freeze v2 Figure 1 Cell state enrichment analysis: version 22 January 2025

### Pipeline for integration of all cells from all AIDA libraries

##### Preparing environment for dataset processing ##### 

library(gridExtra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(Seurat)

##### USER INPUT REQUIRED: FIELD 1 libraries of interest #####
### Enter working directory, and locations of RCAv2-related files ###
### Enter a a list of vectors of library names ###
### For each list member, loop through each library for their singlets ###

setwd("/mnt/sdd/AIDA_Pipeline_v2/Step02_Production/")

df_global_proj_immune <- read.table("/mnt/sdd/AIDA_Pipeline_v2/files_global/RCAv2_rownames_of_global_projection_immune_cells.txt", 
                                    header = FALSE)

##### INPUT REQUIRED: READ IN METADATA FOR APPROVED, CONSENTED DONORS AS df_consented_metadata #####
##### INPUT REQUIRED: READ IN FILE OF APPROVED, CONSENTED DONORS AS df_consented_genotypes #####

df_consented_metadata <- read.table(file = "/mnt/sdd/AIDA_Pipeline_v2/files_global/AIDA_Phase1_DonorMetadata_ApprovedGenotypes.txt", sep = "\t", header = TRUE)
rownames(df_consented_metadata) <- df_consented_metadata$DCP_ID
df_consented_genotypes <- read.table(file = "/mnt/sdd/AIDA_Pipeline_v2/files_global/AIDA_Phase1_ApprovedDonorGenotypes.txt", sep = "\t", header = TRUE)

# Wrangle metadata and genotypes files into single data frame

df_consented_genotypes <- df_consented_genotypes[, c("Single_cell_ID", "Genotyping_ID")]
df_consented_genotypes <- unique(df_consented_genotypes)
for (i in colnames(df_consented_metadata)[2:dim(df_consented_metadata)[2]]){
  df_consented_genotypes[, i] <- df_consented_metadata[df_consented_genotypes$Single_cell_ID, i]
}
rownames(df_consented_genotypes) <- df_consented_genotypes$Genotyping_ID

### Sanity checked through sample of df_consented_genotypes

df_cells_filter_QC <- read.table(file = "/mnt/sdd/AIDA_Pipeline_v2/Step02_Production/AIDA_Phase2_DataFrame_RCAv2_QC_Filter_Status.txt", sep = "\t", header = TRUE)
rownames(df_cells_filter_QC) <- df_cells_filter_QC$barcode

list_AIDA_libraries <- list()

file_libraries <- read.table("/mnt/sdd/AIDA_Pipeline_v2/files_global/AIDA_Phase1_CodePipeline_Assignment.txt", 
                             header = TRUE, sep = "\t")
vec_library_names <- file_libraries$Library

vec_gene_names <- c()

for (str_name in vec_library_names){
  
  ##### USER INPUT REQUIRED: FIELD 2 input directories: Genes and QC metrics #####
  ### Enter details of directories containing input files from AIDA_Pipeline_v2_Step01 ###
  ### Read in genes that are present at >=0.1% of cells per library ###
  ### Read in data frame containing metadata and QC metrics per library (pMito, NODG, nUMI) ###
  ### Find union of genes across libraries for RCAv2 projection ###
  
  df_gene_names <- read.table(file = paste0("/mnt/sdd/AIDA_Pipeline_v2/data/", str_name, "_5GEX/", str_name, "_DataFrame_Genes_passing_filter.txt"),
                              header = TRUE, sep = "\t")
  vec_gene_names <- c(vec_gene_names, df_gene_names$Genes)
  df_metadata_genotype <- read.table(file = paste0("/mnt/sdd/AIDA_Pipeline_v2/data/", str_name, "_5GEX/", str_name, "_DataFrame_Singlets_SingletDoublet_QC_Annotations.txt"),
                                     header = TRUE, sep = "\t")
  df_metadata_genotype <- df_metadata_genotype[, c("barcode_name", "orig.ident", "Genotyping_ID")]
  
  if (match(str_name, vec_library_names) == 1){
    df_all_metadata_genotype  <- df_metadata_genotype
  }
  if (match(str_name, vec_library_names) > 1){
    df_all_metadata_genotype <- rbind(df_all_metadata_genotype, df_metadata_genotype)
  }
}

vec_gene_names <- unique(vec_gene_names)

rownames(df_all_metadata_genotype) <- df_all_metadata_genotype$barcode_name
df_all_metadata_genotype$Genotyping_ID_noprefix <- sapply(df_all_metadata_genotype$Genotyping_ID, 
                                                          FUN = function(x){paste0(strsplit(x, split = "_")[[1]][-1], collapse = "_")})
df_cells_filter_QC$Country <- NA
df_cells_filter_QC$Library <- NA
df_cells_filter_QC$Genotyping_ID <- NA
df_cells_filter_QC$Library <- df_all_metadata_genotype[rownames(df_cells_filter_QC), "orig.ident"]
df_cells_filter_QC$Genotyping_ID <- df_all_metadata_genotype[rownames(df_cells_filter_QC), "Genotyping_ID_noprefix"]
df_cells_filter_QC$Country <- sapply(df_cells_filter_QC$Library, 
                                     FUN = function(x){strsplit(x, split = "_")[[1]][1]})
sum(is.na(df_cells_filter_QC$Country))
sum(is.na(df_cells_filter_QC$Library))
sum(is.na(df_cells_filter_QC$Genotyping_ID))

df_all_metadata_genotype$Country <- NA
df_all_metadata_genotype$Country <- sapply(df_all_metadata_genotype$orig.ident, 
                                           FUN = function(x){strsplit(x, split = "_")[[1]][1]})

##### For overview of integrated object #####

vec_barcodes_to_keep <- df_cells_filter_QC$barcode[which(df_cells_filter_QC$QC_Filter == "Retain")]

for (str_name in vec_library_names[1:length(vec_library_names)]){
  seurat_library <- readRDS(paste0("/mnt/sdd/AIDA_Pipeline_v2/data/", str_name, "_5GEX/", str_name, "_Seurat_Consented_Singlets_AllGenes.RDS"))
  vec_library_barcodes <- colnames(seurat_library)[colnames(seurat_library) %in% vec_barcodes_to_keep]
  # For keeping genes expressed at at least 0.1% of cells in at least one library 
  seurat_library_subset <- seurat_library[vec_gene_names, vec_library_barcodes]
  seurat_library_subset$RCAv2_QC <- df_cells_filter_QC[colnames(seurat_library_subset), "RCAv2_QC"]
  seurat_library_subset$RCA.proj.cell <- df_cells_filter_QC[colnames(seurat_library_subset), "RCA.proj.cell"]
  seurat_library_subset$pMito <- df_cells_filter_QC[colnames(seurat_library_subset), "pMito"]
  seurat_library_subset$NODG <- df_cells_filter_QC[colnames(seurat_library_subset), "NODG"]
  seurat_library_subset$nUMI <- df_cells_filter_QC[colnames(seurat_library_subset), "nUMI"]
  seurat_library_subset$Country <- df_cells_filter_QC[colnames(seurat_library_subset), "Country"]
  seurat_library_subset$Library <- df_cells_filter_QC[colnames(seurat_library_subset), "Library"]
  seurat_library_subset$Genotyping_ID <- df_cells_filter_QC[colnames(seurat_library_subset), "Genotyping_ID"]
  seurat_library_subset$DCP_ID <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "Single_cell_ID"]
  seurat_library_subset$BMI <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "BMI"]
  seurat_library_subset$Age <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "Age"]
  seurat_library_subset$Sex <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "Sex"]
  seurat_library_subset$ethnicity <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "ethnicity"]
  # For saving Seurat libraries as list members for integration 
  list_AIDA_libraries[[str_name]] <- seurat_library_subset
  print(str_name)
}

# For saving list of Seurat libraries for integration
saveRDS(list_AIDA_libraries, "Step03_list_AIDA_libraries.RDS")
# list_AIDA_libraries <- readRDS("Step03_list_AIDA_libraries.RDS")

# Identify library with highest number of cells passing cell type-specific QC

vec_names_list_AIDA_libraries <- names(list_AIDA_libraries)
vec_numbers_list_AIDA_libraries <- c()

for (i in 1:length(list_AIDA_libraries)){
  print(paste0(names(list_AIDA_libraries)[i], ": ", 
               ncol(list_AIDA_libraries[[i]]), " cells"))
  vec_numbers_list_AIDA_libraries <- c(vec_numbers_list_AIDA_libraries, 
                                       ncol(list_AIDA_libraries[[i]]))
}

df_cells_per_libraries <- data.frame(Library = vec_names_list_AIDA_libraries, 
                                     Cells = vec_numbers_list_AIDA_libraries)

str_reference_dataset <- df_cells_per_libraries$Library[which(df_cells_per_libraries$Cells == max(df_cells_per_libraries$Cells))]

# Normalise and identify variable features

list_AIDA_libraries <- lapply(X = list_AIDA_libraries, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
}
)

# Select features for integration, and run PCA on each Seurat object in the list

vec_features <- SelectIntegrationFeatures(object.list = list_AIDA_libraries)
saveRDS(vec_features, file = "Step03_IntegrationVariableFeatures_allbatches.RDS")

# Reduce each individual library to the integration features
# This is required due to the overall size of the dataset - if we retain all
# the ~19k to ~20k features per library that are expressed at >0.1% in at least
# one library, the dataset becomes too large and hits an error, such as 
### Cholmod error 'problem too large' at file ../Core/cholmod_sparse.c, line 89

list_AIDA_libraries <- list()

for (str_name in vec_library_names){
  seurat_library <- readRDS(paste0("/mnt/sdd/AIDA_Pipeline_v2/data/", str_name, "_5GEX/", str_name, "_Seurat_Consented_Singlets_AllGenes.RDS"))
  vec_library_barcodes <- colnames(seurat_library)[colnames(seurat_library) %in% vec_barcodes_to_keep]
  seurat_library_subset <- seurat_library[vec_features, vec_library_barcodes]
  seurat_library_subset$RCAv2_QC <- df_cells_filter_QC[colnames(seurat_library_subset), "RCAv2_QC"]
  seurat_library_subset$RCA.proj.cell <- df_cells_filter_QC[colnames(seurat_library_subset), "RCA.proj.cell"]
  seurat_library_subset$pMito <- df_cells_filter_QC[colnames(seurat_library_subset), "pMito"]
  seurat_library_subset$NODG <- df_cells_filter_QC[colnames(seurat_library_subset), "NODG"]
  seurat_library_subset$nUMI <- df_cells_filter_QC[colnames(seurat_library_subset), "nUMI"]
  seurat_library_subset$Country <- df_cells_filter_QC[colnames(seurat_library_subset), "Country"]
  seurat_library_subset$Library <- df_cells_filter_QC[colnames(seurat_library_subset), "Library"]
  seurat_library_subset$Genotyping_ID <- df_cells_filter_QC[colnames(seurat_library_subset), "Genotyping_ID"]
  seurat_library_subset$DCP_ID <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "Single_cell_ID"]
  seurat_library_subset$BMI <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "BMI"]
  seurat_library_subset$Age <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "Age"]
  seurat_library_subset$Sex <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "Sex"]
  seurat_library_subset$ethnicity <- df_consented_genotypes[seurat_library_subset$Genotyping_ID, "ethnicity"]
  seurat_library_subset@assays$RNA@var.features <- vec_features
  list_AIDA_libraries[[str_name]] <- seurat_library_subset
  print(str_name)
}

list_AIDA_libraries <- lapply(X = list_AIDA_libraries, FUN = function(x){
  x <- ScaleData(x, features = vec_features)
  x <- RunPCA(x, features = vec_features)
}
)

# Use library with highest number of cells as reference dataset for integration

reference_dataset <- which(names(list_AIDA_libraries) == str_reference_dataset)
anchors <- FindIntegrationAnchors(object.list = list_AIDA_libraries, 
                                  reference = reference_dataset, 
                                  reduction = "rpca", 
                                  dims = 1:30)
saveRDS(anchors, "Step03_Integration_Anchors_RPCA.rds")

# Integrate data

integrated_PBMC <- IntegrateData(anchorset = anchors, dims = 1:30)

### Subset out India and European cells for this analysis

df_AIDA_metadata <- integrated_PBMC[[]]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$Country == "IN"), ]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$ethnicity == "European"), ]
seurat_AIDA_subset <- integrated_PBMC[ , rownames(df_AIDA_metadata)]

seurat_AIDA_subset <- FindNeighbors(object = seurat_AIDA_subset, 
                                    reduction = "pca", dims = 1:30, 
                                    k.param = 500, verbose = TRUE, 
                                    compute.SNN = FALSE)
### Above step takes ~1-2 hours...

graph_nnMat_final <- Matrix::Matrix(seurat_AIDA_subset@meta.data[, "integrated_nn"], sparse = TRUE)
seurat_AIDA_subset$UMAP_1 <- seurat_AIDA_subset@reductions$umap@cell.embeddings[, 1]
seurat_AIDA_subset$UMAP_2 <- seurat_AIDA_subset@reductions$umap@cell.embeddings[, 2]

### Comparisons of ethnicities

vec_Asian_ethnicity <- c("Japanese", "Korean", "SG_Chinese", "SG_Malay", "SG_Indian", "Thai")

for (ref in vec_Asian_ethnicity){
  
  ref_num <- length(which(seurat_AIDA_subset$ethnicity==ref))
  
  # calculate number of reference cells in each cell's 500 neighbours
  graph_nnMat_final.ref <- graph_nnMat_final[, names(which(seurat_AIDA_subset$ethnicity==ref))]
  graph_nnMat_final.ref.sum <- Matrix::rowSums(graph_nnMat_final.ref)
  
  # calculate number of other cells in each cell's 500 neighbours
  graph_nnMat_final.comp <- graph_nnMat_final[, names(which(seurat_AIDA_subset$ethnicity!=ref))]
  graph_nnMat_final.comp.sum <- Matrix::rowSums(graph_nnMat_final.comp)
  
  comp.ref.df <- data.frame(comp = graph_nnMat_final.comp.sum,
                            ref = graph_nnMat_final.ref.sum)
  
  comp_FC <- length(which(seurat_AIDA_subset$ethnicity!=ref)) / ref_num
  comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
  comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)
  
  comp.ref.df[(comp.ref.df$comp == 0 & comp.ref.df$ref == 0), "log2Ratio"] <- 0
  
  ### Flipping the sign to look for enrichment of cells from reference
  seurat_AIDA_subset$log2Ratio <- -comp.ref.df$log2Ratio
  seurat_AIDA_subset$log2Ratio[which(seurat_AIDA_subset$log2Ratio > 2)] <- 2
  seurat_AIDA_subset$log2Ratio[which(seurat_AIDA_subset$log2Ratio < -2)] <- -2
  
  ggplot()+
    geom_point(seurat_AIDA_subset[[]], 
               mapping = aes(x=UMAP_1, y=UMAP_2, 
                             color=log2Ratio), 
               size=0.2, stroke=0)+
    scale_color_gradientn(colors = c("darkorange","#e8ded3","#E8E8E8","#d5e0e8","darkblue" ),
                          values=c(1, .54, .5, .46, 0),
                          limits=c(-2,2))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) + 
    ggtitle(paste0("Cell state enrichment in ", ref, " donors"))
  ggsave(paste0("Figure1_AIDA_CellState_500kNN_", ref, ".png"), device = "png",
         width = 7, height = 6)
  
}

### Looking at cell state enrichment for Male versus Female sex

ref <- "Female"
ref_num <- length(which(seurat_AIDA_subset$Sex==ref))

# calculate number of reference cells in each cell's 500 neighbours
graph_nnMat_final.ref <- graph_nnMat_final[, names(which(seurat_AIDA_subset$Sex==ref))]
graph_nnMat_final.ref.sum <- Matrix::rowSums(graph_nnMat_final.ref)

# calculate number of other cells in each cell's 500 neighbours
graph_nnMat_final.comp <- graph_nnMat_final[, names(which(seurat_AIDA_subset$Sex!=ref))]
graph_nnMat_final.comp.sum <- Matrix::rowSums(graph_nnMat_final.comp)

comp.ref.df <- data.frame(comp = graph_nnMat_final.comp.sum,
                          ref = graph_nnMat_final.ref.sum)

comp_FC <- length(which(seurat_AIDA_subset$Sex!=ref)) / ref_num
comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)

comp.ref.df[(comp.ref.df$comp == 0 & comp.ref.df$ref == 0), "log2Ratio"] <- 0

### Flipping the sign to look for enrichment of cells from reference
seurat_AIDA_subset$log2Ratio <- -comp.ref.df$log2Ratio
seurat_AIDA_subset$log2Ratio[which(seurat_AIDA_subset$log2Ratio > 2)] <- 2
seurat_AIDA_subset$log2Ratio[which(seurat_AIDA_subset$log2Ratio < -2)] <- -2

ggplot()+
  geom_point(seurat_AIDA_subset[[]], 
             mapping = aes(x=UMAP_1, y=UMAP_2, 
                           color=log2Ratio), 
             size=0.2, stroke=0)+
  scale_color_gradientn(colors = c("darkorange", "#e8ded3", "#E8E8E8", "#d5e0e8", 
                                   "darkblue" ),
                        values=c(1, .54, .5, .46, 0),
                        limits=c(-2, 2))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Cell state enrichment in ", ref, " donors"))
ggsave(paste0("Figure1_AIDA_CellState_500kNN_", ref, ".png"), device = "png",
       width = 7, height = 6)

seurat_AIDA_subset$age_bin <- 1
seurat_AIDA_subset$Age <- as.numeric(seurat_AIDA_subset$Age)

vec_donor_ages <- c()

for (donor in unique(seurat_AIDA_subset$DCP_ID)){
  vec_donor_ages <- c(vec_donor_ages, unique(seurat_AIDA_subset$Age[which(seurat_AIDA_subset$DCP_ID == donor)]))
}

# unique(seurat_AIDA_subset$DCP_ID, seurat_AIDA_subset$Age)
summary(seurat_AIDA_subset$Age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 19.00   32.00   40.00   41.58   49.00   77.00 

summary(vec_donor_ages)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 19.00   33.00   41.00   42.04   50.00   77.00 

seurat_AIDA_subset$age_bin[which(seurat_AIDA_subset$Age >= 19 & seurat_AIDA_subset$Age < 33)] <- "Ages_19_to_32"
seurat_AIDA_subset$age_bin[which(seurat_AIDA_subset$Age >= 33 & seurat_AIDA_subset$Age < 41)] <- "Ages_33_to_40"
seurat_AIDA_subset$age_bin[which(seurat_AIDA_subset$Age >= 41 & seurat_AIDA_subset$Age < 50)] <- "Ages_41_to_49"
seurat_AIDA_subset$age_bin[which(seurat_AIDA_subset$Age >= 50 & seurat_AIDA_subset$Age <= 77)] <- "Ages_50_to_77"

vec_age <- c("Ages_19_to_32", "Ages_33_to_40", "Ages_41_to_49", "Ages_50_to_77")

for (ref in vec_age){
  
  ref_num <- length(which(seurat_AIDA_subset$age_bin==ref))
  
  # calculate number of reference cells in each cell's 500 neighbours
  graph_nnMat_final.ref <- graph_nnMat_final[, names(which(seurat_AIDA_subset$age_bin==ref))]
  graph_nnMat_final.ref.sum <- Matrix::rowSums(graph_nnMat_final.ref)
  
  # calculate number of other cells in each cell's 500 neighbours
  graph_nnMat_final.comp <- graph_nnMat_final[, names(which(seurat_AIDA_subset$age_bin!=ref))]
  graph_nnMat_final.comp.sum <- Matrix::rowSums(graph_nnMat_final.comp)
  
  comp.ref.df <- data.frame(comp = graph_nnMat_final.comp.sum,
                            ref = graph_nnMat_final.ref.sum)
  
  comp_FC <- length(which(seurat_AIDA_subset$age_bin!=ref)) / ref_num
  comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
  comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)
  
  comp.ref.df[(comp.ref.df$comp == 0 & comp.ref.df$ref == 0), "log2Ratio"] <- 0
  
  ### Flipping the sign to look for enrichment of cells from reference
  seurat_AIDA_subset$log2Ratio <- -comp.ref.df$log2Ratio
  seurat_AIDA_subset$log2Ratio[which(seurat_AIDA_subset$log2Ratio > 2)] <- 2
  seurat_AIDA_subset$log2Ratio[which(seurat_AIDA_subset$log2Ratio < -2)] <- -2
  
  ggplot()+
    geom_point(seurat_AIDA_subset[[]], 
               mapping = aes(x=UMAP_1, y=UMAP_2, 
                             color=log2Ratio), 
               size=0.2, stroke=0)+
    scale_color_gradientn(colors = c("darkorange", "#e8ded3", "#E8E8E8", "#d5e0e8", 
                                     "darkblue" ),
                          values=c(1, .54, .5, .46, 0),
                          limits=c(-2,2))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) + 
    ggtitle(paste0("Cell state enrichment in donors spanning ", ref))
  ggsave(paste0("Figure1_AIDA_CellState_500kNN_", ref, ".png"), device = "png",
         width = 7, height = 6)
  
}
