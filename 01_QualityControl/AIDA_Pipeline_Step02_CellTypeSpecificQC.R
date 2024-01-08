### Standard pipeline for cell type-specific QC of all AIDA libraries

##### Preparing environment for dataset processing ##### 

library(RCAv2)
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
df_RCAv2_to_RCAv2QC <- read.table(file = "/mnt/sdd/AIDA_Pipeline_v2/files_global/RCAv2_Names_to_RCAv2_QC_Annotations_Immune_Global_Novershtern_Monaco.txt", sep = "\t", header = TRUE)
rownames(df_RCAv2_to_RCAv2QC) <- df_RCAv2_to_RCAv2QC$RCAv2_Names

list_AIDA_libraries <- list()

file_libraries <- read.table("/mnt/sdd/AIDA_Pipeline_v2/files_global/AIDA_Phase1_CodePipeline_Assignment.txt", 
                             header = TRUE, sep = "\t")

for (str_code_part in unique(file_libraries$CellTypeQC_Part)){
  list_AIDA_libraries[[str_code_part]] <- file_libraries$Library[which(file_libraries$CellTypeQC_Part == str_code_part)]
}

for (str_name_in_list in names(list_AIDA_libraries)){
  
  vec_str_name <- list_AIDA_libraries[[str_name_in_list]]
  vec_gene_names <- c()
  
  ##### USER INPUT REQUIRED: FIELD 2 input directories: Genes and QC metrics #####
  ### Enter details of directories containing input files from AIDA_Pipeline_v2_Step01 ###
  ### Read in genes that are present at >=0.1% of cells per library ###
  ### Read in data frame containing QC metrics per library (pMito, NODG, nUMI) ###
  ### Find union of genes across libraries for RCAv2 projection ###
  
  for (str_name in vec_str_name){
    
    df_gene_names <- read.table(file = paste0("/mnt/sdd/AIDA_Pipeline_v2/data/", str_name, "_5GEX/", str_name, "_DataFrame_Genes_passing_filter.txt"),
                                header = TRUE, sep = "\t")
    vec_gene_names <- c(vec_gene_names, df_gene_names$Genes)
    df_QC <- read.table(file = paste0("/mnt/sdd/AIDA_Pipeline_v2/data/", str_name, "_5GEX/", str_name, "_DataFrame_Singlets_SingletDoublet_QC_Annotations.txt"),
                        header = TRUE, sep = "\t")
    df_QC <- df_QC[, c("orig.ident", "barcode_name", 
                       "pMito", "NODG", "nUMI")]
    
    if (match(str_name, vec_str_name) == 1){
      df_all_QC  <- df_QC
    }
    if (match(str_name, vec_str_name) > 1){
      df_all_QC <- rbind(df_all_QC, df_QC)
    }
    
  }
  
  vec_gene_names <- unique(vec_gene_names)
  
  # Remove genes "^MT-|^RPS|^RPL|^HSP" from RCAv2 projection step
  vec_gene_names_to_keep <- vec_gene_names[grep(pattern = "^MT-|^RPS|^RPL|^HSP", 
                                                x = vec_gene_names, invert = TRUE)]
  
  for (str_name in vec_str_name){
    
    ##### USER INPUT REQUIRED: FIELD 3 input directories: Singlets from each library #####
    ### Enter details of directories containing input files from AIDA_Pipeline_v2_Step01 ###
    ### Read in RCA object for each library ###
    ### Combine gene-cell matrices for each subset of libraries ###
    
    rca_PBMC <- readRDS(paste0("/mnt/sdd/AIDA_Pipeline_v2/data/", str_name, "_5GEX/", str_name, "_RCA_after_doublet_and_RBC_removal.rds"))
    
    matrix_rawdata <- rca_PBMC$raw.data[vec_gene_names_to_keep, ]
    
    if (match(str_name, vec_str_name) == 1){
      matrix_all_rawdata  <- matrix_rawdata
    }
    if (match(str_name, vec_str_name) > 1){
      matrix_all_rawdata <- Seurat::RowMergeSparseMatrices(matrix_all_rawdata, matrix_rawdata)
    }
    
  }
  
  RCA_all <- createRCAObject(rawData = matrix_all_rawdata)
  
  # Normalise data
  
  RCA_all <- dataLogNormalise(RCA_all)
  
  ##### PRINTING: DATA NORMALISATION #####
  print(str_name_in_list)
  print("Normalisation of data completed.")
  
  # Obtain RCAv2 projection results against immune component in RCAv2 global panel
  
  RCA_all <- dataProject(RCA_all, method = "GlobalPanel",
                         corMeth = "pearson", scale = TRUE)
  df_global_proj <- as.data.frame(RCA_all$projection.data)
  df_global_proj <- df_global_proj[df_global_proj_immune$V1, ]
  
  ##### PRINTING: RCAv2 GLOBAL PROJECTION #####
  print(str_name_in_list)
  print("RCAv2 global panel projection completed.")
  
  # Obtain RCAv2 projection results against Novershtern and Monaco panels
  
  RCA_all <- dataProjectMultiPanel(RCA_all, 
                                   method = list("NovershternPanel", "MonacoPanel"),
                                   scale = TRUE, corMeth = "pearson")
  df_Novershtern_Monaco_proj <- as.data.frame(RCA_all$projection.data)
  
  ##### PRINTING: RCAv2 PROJECTION #####
  print(str_name_in_list)
  print("RCAv2 Novershtern and Monaco panels projection completed.")
  
  # Combine immune component in RCAv2 global panel, Novershtern, and Monaco panels
  
  df_all_projection <- rbind(df_global_proj, df_Novershtern_Monaco_proj)
  df_all_projection <- as.matrix(df_all_projection)
  df_all_projection <- as(df_all_projection, "dgCMatrix")
  
  # Assign combined RCAv2 projection result to RCA object
  
  RCA_all$projection.data <- df_all_projection
  
  # Estimate the most probable cell type label for each cell using optimised code (apply)
  
  # RCA_all <- estimateCellTypeFromProjection(RCA_all, confidence = NULL)
  
  function_estimate_celltype <- function(x){
    return(names(x)[which(x == max(x))])
  }
  
  RCA_all$cell.Type.Estimate.per.cell <- as.list(apply(X = RCA_all$projection.data, 
                                                       MARGIN = 2, 
                                                       FUN = function_estimate_celltype))
  
  ##### PRINTING: RCAv2 CELL TYPE ESTIMATION #####
  print(str_name_in_list)
  print("Estimation of cell type completed.")
  
  ##### Perform PCA and clustering in RCAv2 projection space #####
  
  seurat_RCA_all <- Seurat::CreateSeuratObject(RCA_all$raw.data)
  
  # Code for PCA of RCAv2 projection space adapted from Seurat RunPCA function, 
  # which by default computes the PCA on the cell x gene matrix, 
  # hence the transposed projection data (cell x RCAv2 cell type) being used in irlba
  
  npcs <- 30
  npcs <- min(npcs, nrow(RCA_all$projection.data) - 1)
  pca.results <- irlba::irlba(A = t(RCA_all$projection.data), nv = npcs)
  feature.loadings <- pca.results$v
  sdev <- pca.results$d/sqrt(max(1, ncol(RCA_all$projection.data) - 1))
  projection <- pca.results$u %*% diag(pca.results$d)
  
  rownames(x = feature.loadings) <- rownames(x = RCA_all$projection.data)
  colnames(x = feature.loadings) <- paste0("PC_", 1:npcs)
  rownames(x = projection) <- colnames(x = RCA_all$projection.data)
  colnames(x = projection) <- colnames(x = feature.loadings)
  
  # Save the sum of (variances per cell across all RCAv2 cell type projections) under misc
  
  total.variance <- sum(apply(X = RCA_all$projection.data, MARGIN = 2, FUN=var))
  seurat_RCA_all@reductions[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = projection,
    loadings = feature.loadings,
    assay = "RNA",
    stdev = sdev,
    key = "PC_",
    misc = list(total.variance = total.variance))
  
  png(paste0(str_name_in_list, "_RCA_all_ElbowPlot.png"))
  print(ElbowPlot(object = seurat_RCA_all, ndims = 30))
  dev.off()
  
  ##### PRINTING: PCA USING irlba #####
  print(str_name_in_list)
  print("PCA completed.")
  
  # Perform clustering in RCAv2 projection space
  
  seurat_RCA_all <- Seurat::FindNeighbors(object = seurat_RCA_all, dims = 1:20)
  
  # Adjust clustering resolution based on number of cells in dataset
  
  if (length(colnames(seurat_RCA_all)) >= 200000){
    str_resolution <- 3
  } 
  if (length(colnames(seurat_RCA_all)) > 300000){
    str_resolution <- 4
  } 
  if (length(colnames(seurat_RCA_all)) < 200000){
    str_resolution <- 2
  }
  
  print(str_name_in_list)
  print("Number of cells: ")
  print(length(colnames(seurat_RCA_all)))
  print("Clustering resolution: ")
  print(str_resolution)
  
  seurat_RCA_all <- Seurat::FindClusters(seurat_RCA_all, resolution = str_resolution)
  
  # Annotate each cluster from RCAv2 projection space by majority vote of RCAv2 label 
  
  df_metadata <- seurat_RCA_all[[]]
  df_cell_type_estimate <- data.frame(row.names = colnames(RCA_all$raw.data),
                                      RCA.proj.cell = unlist(RCA_all$cell.Type.Estimate.per.cell))
  df_merged_cluster_celltypeestimate <- merge(df_metadata, df_cell_type_estimate, 
                                              by = "row.names")
  df_merged_cluster_celltypeestimate <- df_merged_cluster_celltypeestimate[, c("seurat_clusters", "RCA.proj.cell")]
  
  df_merged_cluster_celltypeestimate$count <- 1
  df_merged_cluster_celltypeestimate <- aggregate(count ~ ., df_merged_cluster_celltypeestimate, FUN = sum)
  df_merged_cluster_celltypeestimate <- df_merged_cluster_celltypeestimate %>% group_by(seurat_clusters) %>% top_n(n = 5, wt = count)
  df_merged_cluster_celltypeestimate <- df_merged_cluster_celltypeestimate[order(df_merged_cluster_celltypeestimate$seurat_clusters), ]
  
  df_RCAv2_cluster_to_RCAv2_QC <- FetchData(seurat_RCA_all, vars = c("seurat_clusters"))
  
  # Identify majority vote of RCAv2 annotation for each cluster (RCAv2 projection data clustering) 
  
  vec_RCAv2_cluster_annotations <- (df_merged_cluster_celltypeestimate %>% filter(count == max(count)))$RCA.proj.cell
  
  # Check through annotations to verify that per-cluster annotations look reasonable
  
  sink(file = paste(str_name_in_list, "RCAv2_Annotation_per_RCAv2_space_Cluster.txt", sep = "_"))
  print(str_name_in_list)
  for (i in 0:(length(vec_RCAv2_cluster_annotations)-1)){
    print(df_merged_cluster_celltypeestimate %>% filter(seurat_clusters == i))
  }
  sink()
  
  # Add RCAv2 majority vote label (per cluster based on RCAv2 projection data) to each cell
  
  df_RCAv2_cluster_to_RCAv2_QC$RCAv2_QC <- NA
  
  for (i in 0:(length(vec_RCAv2_cluster_annotations)-1)){
    df_RCAv2_cluster_to_RCAv2_QC[which(df_RCAv2_cluster_to_RCAv2_QC$seurat_clusters==i), "RCAv2_QC"] <- df_RCAv2_to_RCAv2QC[vec_RCAv2_cluster_annotations[i+1], "RCAv2_QC"]
  }
  
  ##### Save data frame of the RCAv2 majority-vote label and QC metrics (pMito, NODG, nUMI) per cell #####
  
  rownames(df_all_QC) <- df_all_QC$barcode_name
  df_RCAv2_cluster_to_RCAv2_QC$pMito <- df_all_QC[rownames(df_RCAv2_cluster_to_RCAv2_QC), "pMito"]
  df_RCAv2_cluster_to_RCAv2_QC$NODG <- df_all_QC[rownames(df_RCAv2_cluster_to_RCAv2_QC), "NODG"]
  df_RCAv2_cluster_to_RCAv2_QC$nUMI <- df_all_QC[rownames(df_RCAv2_cluster_to_RCAv2_QC), "nUMI"]
  df_RCAv2_cluster_to_RCAv2_QC$RCA.proj.cell <- df_cell_type_estimate[rownames(df_RCAv2_cluster_to_RCAv2_QC), "RCA.proj.cell"]
  df_RCAv2_cluster_to_RCAv2_QC$barcode <- rownames(df_RCAv2_cluster_to_RCAv2_QC)
  
  write.table(df_RCAv2_cluster_to_RCAv2_QC, paste(str_name_in_list, "DataFrame_RCAv2_QC_Annotation_pMito_NODG_nUMI.txt", sep = "_"), 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  ##### Plot marker gene DotPlots as well as UMAPs by cluster and by RCAv2 annotation label #####
  
  seurat_RCA_all$RCAv2_QC <- df_RCAv2_cluster_to_RCAv2_QC[colnames(seurat_RCA_all), "RCAv2_QC"]
  seurat_RCA_all <- NormalizeData(seurat_RCA_all, normalization.method = "LogNormalize", scale.factor = 10000)
  
  list_features <- list()
  list_features[["Monocytes, DCs"]] <- c("CD14", "FCGR3A", "CLEC9A", "CLEC10A", "SIGLEC6", "ITM2C")
  list_features[["Plasma B, B"]] <- c("MZB1", "MS4A1", "TCL1A", "TNFRSF17")
  list_features[["T"]] <- c("CD3D", "TRGV9", "GZMK", "GZMB", "CD4", "IL7R", 
                            "CCR7", "FOXP3", "CD8A", "TRAV1-2", "ITGB1", "KLRB1")
  list_features[["NK"]] <- c("NKG7", "GNLY", "NCAM1")
  list_features[["Others"]] <- c("PPBP", "HBB", "CD34")
  
  pdf(paste0(str_name_in_list, "_DotPlot_Marker_Genes_Clusters.pdf"), width = 25, height = 15)
  print(DotPlot(seurat_RCA_all, features = list_features, group.by = "seurat_clusters"))
  dev.off()
  pdf(paste0(str_name_in_list, "_DotPlot_Marker_Genes_RCAv2QCLabels.pdf"), width = 25, height = 15)
  print(DotPlot(seurat_RCA_all, features = list_features, group.by = "RCAv2_QC"))
  dev.off()
  
  seurat_RCA_all <- RunUMAP(seurat_RCA_all, dims = 1:20)
  
  pdf(paste0(str_name_in_list, "_UMAP_RCAv2QC_clusterlabels.pdf"))
  print(DimPlot(seurat_RCA_all, reduction = "umap", group.by = "seurat_clusters", raster = FALSE) + 
          ggtitle(paste0(str_name_in_list, "UMAP labelled by RCAv2 space clusters")))
  dev.off()
  
  pdf(paste0(str_name_in_list, "_UMAP_RCAv2QC_RCAv2labels.pdf"))
  print(DimPlot(seurat_RCA_all, reduction = "umap", group.by = "RCAv2_QC", raster = FALSE) + 
          ggtitle(paste0(str_name_in_list, "UMAP labelled by RCAv2 annotations")))
  dev.off()
  
  if (match(str_name_in_list, names(list_AIDA_libraries)) == 1){
    df_all_RCAv2_QC  <- df_RCAv2_cluster_to_RCAv2_QC
  }
  if (match(str_name_in_list, names(list_AIDA_libraries)) > 1){
    df_all_RCAv2_QC <- rbind(df_all_RCAv2_QC, df_RCAv2_cluster_to_RCAv2_QC)
  }
  
  ##### PRINTING: EXAMINATION OF DATA SUBSET #####
  print(str_name_in_list)
  print("Examination of dataset subset completed.")
  
}

##### Visualising cell type-specific QC metrics (pMito, NODG) #####

vec_unique_celltype <- unique(df_all_RCAv2_QC$RCAv2_QC)

for (str_celltype in vec_unique_celltype){
  
  df_celltype_of_interest <- df_all_RCAv2_QC[which(df_all_RCAv2_QC$RCAv2_QC == str_celltype), ]
  
  ggplot(df_celltype_of_interest, aes(x = NODG, y = pMito)) +
    geom_point(size = 0.1, stroke = 0) +
    theme_bw() + 
    ggtitle(paste0(str_celltype," n = ", dim(df_celltype_of_interest)[1])) +
    stat_density_2d()
  ggsave(paste0(str_celltype, "_CellTypeSpecific_QC_metrics.png"),
         width = 3, height = 3, device = "png")
  
  ggplot(df_celltype_of_interest, aes(x = NODG, y = pMito)) +
    geom_point(size = 0.1, stroke = 0) +
    theme_bw() + 
    ggtitle(paste0(str_celltype," n = ", dim(df_celltype_of_interest)[1])) +
    stat_density_2d() + 
    ylim(0, 0.2)
  ggsave(paste0(str_celltype, "_CellTypeSpecific_QC_metrics_Zoomin.png"),
         width = 3, height = 3, device = "png")
  
}

##### Cell type-specific QC filters (pMito, NODG) #####

vec_QC_celltype <- c("Myeloid",
                     "B",
                     "T",
                     "NK",
                     "pDC",
                     "CD34_HSPC",
                     "Plasma_Cell",
                     "Platelet")
vec_QC_NODG_min <- c(500, 1000, 1100, 1100, 1700, 1000, 1000, 300)
vec_QC_NODG_max <- c(5500, 3000, 2900, 2900, 4500, 5000, 6500, 1200)
vec_QC_pMito_min <- c(0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001)
vec_QC_pMito_max <- c(0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.125, 0.125)
df_QC_celltype_filters <- data.frame(CellType = vec_QC_celltype,
                                     NODG_min = vec_QC_NODG_min,
                                     NODG_max = vec_QC_NODG_max,
                                     pMito_min = vec_QC_pMito_min,
                                     pMito_max = vec_QC_pMito_max)
rownames(df_QC_celltype_filters) <- df_QC_celltype_filters$CellType

##### Filtering based on cell type-specific QC thresholds (pMito, NODG) #####

vec_barcodes_to_keep <- c()

for (str_celltype in vec_unique_celltype){
  
  print(str_celltype)
  df_celltype_of_interest <- df_all_RCAv2_QC[which(df_all_RCAv2_QC$RCAv2_QC == str_celltype), ]
  vec_barcodes_QC <- rownames(df_celltype_of_interest)[which((df_celltype_of_interest$NODG >= df_QC_celltype_filters[str_celltype, "NODG_min"]) & 
                                                               (df_celltype_of_interest$NODG <= df_QC_celltype_filters[str_celltype, "NODG_max"]) &
                                                               (df_celltype_of_interest$pMito >= df_QC_celltype_filters[str_celltype, "pMito_min"]) &
                                                               (df_celltype_of_interest$pMito <= df_QC_celltype_filters[str_celltype, "pMito_max"]))]
  vec_barcodes_to_keep <- c(vec_barcodes_to_keep, vec_barcodes_QC)
  
  df_celltype_of_interest$QC_Filter <- "Exclude"
  df_celltype_of_interest[vec_barcodes_QC, "QC_Filter"] <- "Retain"
  
  ggplot(df_celltype_of_interest, aes(x = NODG, y = pMito, colour = QC_Filter)) +
    geom_point(size = 0.1, stroke = 0) +
    theme_bw() + 
    ggtitle(paste0(str_celltype, ": retain n = ", length(vec_barcodes_QC), 
                   " (", round(length(vec_barcodes_QC)/dim(df_celltype_of_interest)[1]*100, 1), "%)")) +
    stat_density_2d() + 
    ylim(0, 0.2) + theme(legend.position = "none") 
  ggsave(paste0(str_celltype, "_CellTypeSpecific_QC_metrics_Filter_Zoomin.png"),
         width = 5, height = 5, device = "png")
  
}

##### Save data frame of AIDA Phase 2 cells with QC_Filter status #####

df_all_RCAv2_QC$QC_Filter <- "Exclude"
df_all_RCAv2_QC[vec_barcodes_to_keep, "QC_Filter"] <- "Retain"
write.table(df_all_RCAv2_QC, file = "AIDA_Phase2_DataFrame_RCAv2_QC_Filter_Status.txt", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
