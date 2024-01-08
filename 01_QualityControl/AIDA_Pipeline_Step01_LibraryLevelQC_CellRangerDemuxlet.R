### Standard pipeline for processing of AIDA Cell Ranger + Demuxlet datasets
### Last tested on 13 June 2023 on clean EC2 instance: 
### RStudio Server version 1.3.1073, R 4.0.2, Ubuntu 18.04 LTS

##### Installation of required R packages #####

# install.packages("Seurat")
# install.packages("tidyverse")
# install.packages("remotes")
# library(remotes)
### Production run on RStudio Server version 1.3.1073, R 4.0.2, Ubuntu 18.04 LTS
### required an older version of "Hmisc" for installing WGCNA and RCAv2.
# install_version("Hmisc", version = "4.7.0")
# install_github("prabhakarlab/RCAv2")
# install_github("chris-mcginnis-ucsf/DoubletFinder")
# install.packages("gridExtra")
# install.packages("RColorBrewer")

##### Preparing environment for dataset processing ##### 

library(RCAv2)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(DoubletFinder)
library(Seurat)

### Enter the following variables for the library / libraries of interest
### Ensure that input files are re-named: barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz

##### USER INPUT REQUIRED: FIELD 1 library / libraries of interest #####
### Enter a vector of library names to loop through

vec_str_name <- c("IN_NIB_B001_L001",
                  "IN_NIB_B001_L002",
                  "IN_NIB_B002_L001",
                  "IN_NIB_B002_L002")
# str_name <- vec_str_name[1]

##### USER INPUT REQUIRED: FIELD 2 RCAv2 file locations #####
### Input files containing RCAv2 immune cell types of interest

df_global_proj_immune <- read.table("/mnt/sdd/AIDA_Pipeline_v2/files_global/RCAv2_rownames_of_global_projection_immune_cells.txt", 
                                    header = FALSE)
df_RCAv2_to_annotations <- read.table(file = "/mnt/sdd/AIDA_Pipeline_v2/files_global/RCAv2_Names_to_Annotations_Immune_Global_Novershtern_Monaco.txt", sep = "\t", header = TRUE)
rownames(df_RCAv2_to_annotations) <- df_RCAv2_to_annotations$RCAv2_Names

##### USER INPUT REQUIRED: FIELD 3 consented genotype file location #####
### Input file containing list of genotypes consented for AIDA Phase 1

df_consented_genotypes <- read.table("/mnt/sdd/AIDA_Pipeline_v2/files_global/AIDA_Phase1_ApprovedDonorGenotypes.txt", 
                                     sep = "\t", header = TRUE)

##### USER INPUT REQUIRED: FIELD 4 chr Y gene names file locations #####
### Input files containing lists of chromosome Y gene names

df_chrY_nonPAR_genes <- read.table(file = "/mnt/sdd/AIDA_Pipeline_v2/files_global/list_chrY_nonPAR_genes.txt", 
                                   header = FALSE) 

df_chrY_PAR_genes <- read.table(file = "/mnt/sdd/AIDA_Pipeline_v2/files_global/list_chrY_PAR_genes.txt", 
                                header = FALSE) 

### Loop through the following code for each library

for (str_name in vec_str_name){
  
  ##### USER INPUT REQUIRED: FIELD 5 data and working directories #####
  str_link_dataset <- paste0("/mnt/sdd/AIDA_Pipeline_v2/data/", str_name, "_5GEX/")
  setwd(str_link_dataset)
  df_demuxlet <- read.table(paste0(str_link_dataset, str_name, "_using_CellRanger_BAM_include_introns_true.best"),
                            sep = "\t", header = TRUE)
  
  ##### USER INPUT REQUIRED: FIELD 6 number of samples per batch #####
  # Enter number of samples in this batch
  
  int_number_of_samples <- 16
  
  # Generate an RCA object from a single scRNA-seq library
  
  rca_PBMC <- createRCAObjectFrom10X(paste0(str_link_dataset, "outs/raw_feature_bc_matrix/"))
  
  ##### Examining library for basic QC parameters ##### 
  
  # Plot distribution of number of detected genes (NODG)
  
  matrix_PBMC <- rca_PBMC$raw.data
  vec_NODG <- Matrix::colSums(matrix_PBMC > 0)
  
  # Compute nUMI vector
  
  vec_nUMI <- Matrix::colSums(matrix_PBMC)
  
  # Identify mitochondrial genes, and compute % mitochondrial reads (pMito)
  
  vec_name_mito_genes = grep(pattern = "^MT-", x = rownames(matrix_PBMC), value = TRUE)
  vec_pMito <- Matrix::colSums(matrix_PBMC[vec_name_mito_genes, , drop = FALSE])/Matrix::colSums(matrix_PBMC)
  
  # Plot density plot of pMito against NODG
  
  df_all_NODG_UMI_pMito <- data.frame(NODG = vec_NODG, 
                                      nUMI = vec_nUMI,
                                      pMito = vec_pMito)
  
  ggplot(df_all_NODG_UMI_pMito, aes(x = NODG, y = pMito)) + 
    geom_point(size = 0.1, colour = "grey") +
    geom_density_2d() + 
    theme_bw(11) + 
    xlim(0, 6000) + ylim(0, 0.25) +
    ggtitle(paste0(str_name, " n = ", ncol(rca_PBMC$raw.data)))
  
  ggsave(paste0(str_name, "_Library_QC_Plot01_pMito_NODG.png"), width = 5, height = 5)
  
  # Retain cells with at least 300 detected genes
  # Use a 0.1% filter for genes (retaining genes expressed by at least 0.1% of cells)
  
  rca_PBMC <- dataFilter(rca_PBMC,
                         nGene.thresholds = c(300, Inf),
                         nUMI.thresholds = c(0, Inf),
                         percent.mito.thresholds = c(0, 1),
                         min.cell.exp = 0.001*ncol(rca_PBMC$raw.data), 
                         plot = FALSE)
  
  df_all_NODG_UMI_pMito <- df_all_NODG_UMI_pMito[colnames(rca_PBMC$raw.data), ]
  
  # Plot density plot of pMito against NODG for cells with NODG >= 300
  
  ggplot(df_all_NODG_UMI_pMito, aes(x = NODG, y = pMito)) + 
    geom_point(size = 0.1, colour = "grey") +
    geom_density_2d() + 
    theme_bw(11) + 
    xlim(0, 5000) + ylim(0, 0.1) +
    ggtitle(paste0(str_name, " n = ", ncol(rca_PBMC$raw.data)))
  
  ggsave(paste0(str_name, "_Library_QC_Plot02_pMito_NODG_300.png"), width = 5, height = 5)
  
  ###### Projecting dataset against all RCAv2 immune panels ######
  
  # Retain raw data for export of pMito values for downstream RCAv2 analyses
  
  rca_PBMC_rawdata <- rca_PBMC$raw.data
  
  # Remove "^MT-|^RPS|^RPL|^HSP" genes before RCAv2 projection
  
  vec_gene_names <- rownames(rca_PBMC$raw.data)
  vec_gene_names_keep <- vec_gene_names[grep(pattern = "^MT-|^RPS|^RPL|^HSP", 
                                             x = vec_gene_names, invert = TRUE)]
  rca_PBMC$raw.data <- rca_PBMC$raw.data[vec_gene_names_keep, ]
  
  # Normalise data
  
  rca_PBMC <- dataLogNormalise(rca_PBMC)
  
  # Obtain RCAv2 projection results against immune component in RCAv2 global panel
  
  rca_PBMC <- dataProject(rca_PBMC, method = "GlobalPanel",
                          corMeth = "pearson", scale = TRUE)
  df_global_proj <- as.data.frame(rca_PBMC$projection.data)
  df_global_proj <- df_global_proj[df_global_proj_immune$V1, ]
  
  # Obtain RCAv2 projection results against Novershtern and Monaco panels
  
  rca_PBMC <- dataProjectMultiPanel(rca_PBMC, 
                                    method = list("NovershternPanel", "MonacoPanel"),
                                    scale = TRUE, corMeth = "pearson")
  df_Novershtern_Monaco_proj <- as.data.frame(rca_PBMC$projection.data)
  
  # Combine immune component in RCAv2 global panel, Novershtern, and Monaco panels
  
  df_all_projection <- rbind(df_global_proj, df_Novershtern_Monaco_proj)
  df_all_projection <- as.matrix(df_all_projection)
  df_all_projection <- as(df_all_projection, "dgCMatrix")
  
  # Assign combined RCAv2 projection result to RCA object
  
  rca_PBMC$projection.data <- df_all_projection
  
  # Estimate the most probable cell type label for each cell
  
  rca_PBMC <- estimateCellTypeFromProjection(rca_PBMC, confidence = NULL)
  
  ##### Performing PCA and clustering in RCAv2 projection space #####
  
  seurat_rca_PBMC <- Seurat::CreateSeuratObject(rca_PBMC$raw.data)
  
  # Code for PCA of RCAv2 projection space adapted from Seurat RunPCA function, 
  # which by default computes the PCA on the cell x gene matrix, 
  # hence the transposed projection data (cell x RCAv2 cell type) being used in irlba
  
  npcs <- 30
  npcs <- min(npcs, nrow(rca_PBMC$projection.data) - 1)
  pca.results <- irlba::irlba(A = t(rca_PBMC$projection.data), nv = npcs)
  feature.loadings <- pca.results$v
  sdev <- pca.results$d/sqrt(max(1, ncol(rca_PBMC$projection.data) - 1))
  projection <- pca.results$u %*% diag(pca.results$d)
  
  rownames(x = feature.loadings) <- rownames(x = rca_PBMC$projection.data)
  colnames(x = feature.loadings) <- paste0("PC_", 1:npcs)
  rownames(x = projection) <- colnames(x = rca_PBMC$projection.data)
  colnames(x = projection) <- colnames(x = feature.loadings)
  
  # Save the sum of (variances per cell across all RCAv2 cell type projections) under "misc"
  
  total.variance <- sum(apply(X = rca_PBMC$projection.data, MARGIN = 2, FUN=var))
  seurat_rca_PBMC@reductions[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = projection,
    loadings = feature.loadings,
    assay = "RNA",
    stdev = sdev,
    key = "PC_",
    misc = list(total.variance = total.variance))
  
  png(paste0(str_name, "_RCAv2_space_PCA_ElbowPlot.png"))
  print(ElbowPlot(object = seurat_rca_PBMC, ndims = 30))
  dev.off()
  
  # Perform clustering in RCAv2 projection space
  
  seurat_rca_PBMC <- Seurat::FindNeighbors(object = seurat_rca_PBMC, dims = 1:20)
  seurat_rca_PBMC <- Seurat::FindClusters(seurat_rca_PBMC, resolution = 1)
  
  # Annotate each cluster from RCAv2 projection space by majority vote of RCAv2 label 
  
  df_metadata <- seurat_rca_PBMC[[]]
  df_cell_type_estimate <- data.frame(row.names = colnames(rca_PBMC$raw.data),
                                      RCA.proj.cell = unlist(rca_PBMC$cell.Type.Estimate.per.cell))
  df_merged_cluster_celltypeestimate <- merge(df_metadata, df_cell_type_estimate, 
                                              by = "row.names")
  df_merged_cluster_celltypeestimate <- df_merged_cluster_celltypeestimate[, c("seurat_clusters", "RCA.proj.cell")]
  
  df_merged_cluster_celltypeestimate$count <- 1
  df_merged_cluster_celltypeestimate <- aggregate(count ~ ., df_merged_cluster_celltypeestimate, FUN = sum)
  df_merged_cluster_celltypeestimate <- df_merged_cluster_celltypeestimate %>% group_by(seurat_clusters) %>% top_n(n = 5, wt = count)
  df_merged_cluster_celltypeestimate <- df_merged_cluster_celltypeestimate[order(df_merged_cluster_celltypeestimate$seurat_clusters), ]
  
  df_RCAv2_cluster_to_RCAv2_label <- FetchData(seurat_rca_PBMC, vars = c("seurat_clusters"))
  
  # Identify majority vote of RCAv2 annotation for each cluster (RCAv2 projection data clustering) 
  
  vec_RCAv2_cluster_annotations <- (df_merged_cluster_celltypeestimate %>% filter(count == max(count)))$RCA.proj.cell
  
  # Check through annotations to verify that per-cluster annotations look reasonable
  
  sink(file = paste(str_name, "RCAv2_Annotation_per_RCAv2_space_Cluster.txt", sep = "_"))
  print(str_name)
  for (i in 0:(length(vec_RCAv2_cluster_annotations)-1)){
    print(df_merged_cluster_celltypeestimate %>% filter(seurat_clusters == i))
  }
  sink()
  
  # Add RCAv2 majority vote label (per cluster based on RCAv2 projection data) to each cell
  
  df_RCAv2_cluster_to_RCAv2_label$RCAv2_Cluster_Label <- NA
  
  for (i in 0:(length(vec_RCAv2_cluster_annotations)-1)){
    df_RCAv2_cluster_to_RCAv2_label[which(df_RCAv2_cluster_to_RCAv2_label$seurat_clusters==i), "RCAv2_Cluster_Label"] <- df_RCAv2_to_annotations[vec_RCAv2_cluster_annotations[i+1], "Annotation_Names"]
  }
  
  # Save data frame of the RCAv2 majority-vote label of RCAv2 space clusters and QC metrics (pMito, NODG, nUMI) per cell #
  
  df_RCAv2_cluster_to_RCAv2_label$pMito <- df_all_NODG_UMI_pMito [rownames(df_RCAv2_cluster_to_RCAv2_label), "pMito"]
  df_RCAv2_cluster_to_RCAv2_label$NODG <- df_all_NODG_UMI_pMito [rownames(df_RCAv2_cluster_to_RCAv2_label), "NODG"]
  df_RCAv2_cluster_to_RCAv2_label$nUMI <- df_all_NODG_UMI_pMito [rownames(df_RCAv2_cluster_to_RCAv2_label), "nUMI"]
  df_RCAv2_cluster_to_RCAv2_label$barcode <- rownames(df_RCAv2_cluster_to_RCAv2_label)
  
  write.table(df_RCAv2_cluster_to_RCAv2_label, paste(str_name, "DataFrame_RCAv2_space_Cluster_Annotation_pMito_NODG_nUMI.txt", sep = "_"), 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # Plot marker gene DotPlots and UMAPs by cluster and by RCAv2 annotation label #
  
  seurat_rca_PBMC$RCAv2_Cluster_Label <- df_RCAv2_cluster_to_RCAv2_label[colnames(seurat_rca_PBMC), "RCAv2_Cluster_Label"]
  seurat_rca_PBMC <- NormalizeData(seurat_rca_PBMC, normalization.method = "LogNormalize", scale.factor = 10000)
  
  list_features <- list()
  list_features[["Monocytes, DCs"]] <- c("CD14", "FCGR3A", "CLEC9A", "CLEC10A", "SIGLEC6", "ITM2C")
  list_features[["Plasma B, B"]] <- c("MZB1", "MS4A1", "TCL1A", "TNFRSF17")
  list_features[["T"]] <- c("CD3D", "TRGV9", "GZMK", "GZMB", "CD4", "IL7R", 
                            "CCR7", "FOXP3", "CD8A", "TRAV1-2", "ITGB1", "KLRB1")
  list_features[["NK"]] <- c("NKG7", "GNLY", "NCAM1")
  list_features[["Others"]] <- c("PPBP", "HBB", "CD34")
  
  pdf(paste0(str_name, "_RCAv2_space_DotPlot_Marker_Genes_Clusters.pdf"), width = 25, height = 15)
  print(DotPlot(seurat_rca_PBMC, features = list_features, group.by = "seurat_clusters"))
  dev.off()
  pdf(paste0(str_name, "_RCAv2_space_DotPlot_Marker_Genes_RCAv2Labels.pdf"), width = 25, height = 15)
  print(DotPlot(seurat_rca_PBMC, features = list_features, group.by = "RCAv2_Cluster_Label"))
  dev.off()
  
  seurat_rca_PBMC <- RunUMAP(seurat_rca_PBMC, dims = 1:20)
  
  pdf(paste0(str_name, "_RCAv2_space_UMAP_clusterlabels.pdf"))
  print(DimPlot(seurat_rca_PBMC, reduction = "umap", group.by = "seurat_clusters") + 
          ggtitle(paste0(str_name, " UMAP \nlabelled by RCAv2 space clusters")))
  dev.off()
  
  pdf(paste0(str_name, "_RCAv2_space_UMAP_RCAv2labels.pdf"))
  print(DimPlot(seurat_rca_PBMC, reduction = "umap", group.by = "RCAv2_Cluster_Label") + 
          ggtitle(paste0(str_name, " UMAP \nlabelled by RCAv2 annotations")))
  dev.off()
  
  ##### Identifying genetic doublets from Demuxlet genetic demultiplexing #####
  
  df_demuxlet <- df_demuxlet[which(df_demuxlet$BARCODE %in% colnames(rca_PBMC$raw.data)), ]
  df_type <- as.data.frame(table(df_demuxlet$DROPLET.TYPE))
  
  write.table(df_type, paste0(str_name, "_Demuxlet_summary_AMB_DBL_SNG.txt"),
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  vec_demuxlet_singlet <- df_demuxlet[which(df_demuxlet$DROPLET.TYPE=="SNG"), "BARCODE"]
  
  # Estimate doublet rate: DBL + AMB droplets
  
  number_estimated_doublet_rate <- ((sum(df_demuxlet$DROPLET.TYPE == "DBL") + 
                                       sum(df_demuxlet$DROPLET.TYPE == "AMB"))/length(df_demuxlet$DROPLET.TYPE))*int_number_of_samples/(int_number_of_samples-1)
  
  ##### DoubletFinder #####
  
  # Create Seurat object from full gene-cell matrix (NODG > 300) from library
  # for running DoubletFinder
  
  data_PBMC <- Read10X(data.dir = paste0(str_link_dataset, "outs/raw_feature_bc_matrix/"))
  
  seurat_PBMC <- CreateSeuratObject(counts = data_PBMC, 
                                    project = str_name, 
                                    min.cells = 0, 
                                    min.features = 300)
  
  # seurat_test <- CreateSeuratObject(counts = rca_PBMC$raw.data, 
  #                                   project = "AIDA", 
  #                                   min.cells = 0, 
  #                                   min.features = 0)
  # Warning: Non-unique features (rownames) present in the input matrix, making unique
  
  ### Note that the following gene names appear twice - two ENSG identifiers
  # vec_gene_names[duplicated(vec_gene_names)]
  # [1] "TBCE"     "CYB561D2" "MATR3"    "HSPA14"   "TMSB15B" 
  ### We go with vec_gene_names_keep for retaining genes, to avoid double-counting genes
  
  seurat_PBMC$barcode_name <- paste(colnames(seurat_PBMC), "-", str_name, sep = "")
  # Later on before saving: rename cells to include library information
  # seurat_PBMC <- RenameCells(seurat_PBMC, new.names = seurat_PBMC$barcode_name)
  
  seurat_PBMC_subset <- subset(seurat_PBMC, 
                               features = vec_gene_names)
  
  # Add ground truth (GT) from Demuxlet genetic demultiplexing to Seurat object
  
  df_GT <- data.frame(barcode = colnames(seurat_PBMC_subset))
  rownames(df_GT) <- df_GT$barcode
  df_GT$GT <- "Doublet"
  df_GT[(df_GT$barcode %in% vec_demuxlet_singlet), "GT"] <- "Singlet"
  
  seurat_PBMC_subset$GT <- df_GT[colnames(seurat_PBMC_subset), "GT"]
  
  # Add RCAv2 projection annotation to Seurat object
  
  df_RCA_cell_estimate <- data.frame(barcode = colnames(rca_PBMC$projection.data),
                                     RCAv2_cell = unlist(rca_PBMC$cell.Type.Estimate.per.cell))
  rownames(df_RCA_cell_estimate) <- df_RCA_cell_estimate$barcode
  seurat_PBMC_subset$RCA.proj.cell <- df_RCA_cell_estimate[colnames(seurat_PBMC_subset), "RCAv2_cell"]
  
  # Normalise the data
  
  seurat_PBMC_subset <- NormalizeData(seurat_PBMC_subset, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  
  seurat_PBMC_subset <- FindVariableFeatures(seurat_PBMC_subset, selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  
  vec_all_genes_to_scale <- rownames(seurat_PBMC_subset)
  seurat_PBMC_subset <- ScaleData(seurat_PBMC_subset, features = vec_all_genes_to_scale)
  
  # Run PCA
  
  seurat_PBMC_subset <- RunPCA(seurat_PBMC_subset, features = VariableFeatures(object = seurat_PBMC_subset))
  png(paste0(str_name, "_Library_QC_ElbowPlot.png"))
  print(ElbowPlot(seurat_PBMC_subset, ndims = 50))
  dev.off()
  
  # Perform Louvain clustering
  
  seurat_PBMC_subset <- FindNeighbors(seurat_PBMC_subset, dims = 1:30)
  seurat_PBMC_subset <- FindClusters(seurat_PBMC_subset, resolution = 1)
  
  # Run UMAP
  
  seurat_PBMC_subset <- RunUMAP(seurat_PBMC_subset, dims = 1:30)
  
  pdf(paste0(str_name, "_Library_QC_Plot03_UMAP_Seurat_Clusters.pdf"))
  print(DimPlot(seurat_PBMC_subset, reduction = "umap", label = TRUE) + NoLegend() + 
          ggtitle(paste(str_name, " UMAP labelled by Seurat clusters", sep = "")))
  dev.off()
  
  # Add RCAv2 annotation information to each Seurat cluster (5GEX clustering)
  
  df_RCA_results <- FetchData(seurat_PBMC_subset, vars = c("seurat_clusters", "RCA.proj.cell"))
  
  df_RCA_results$count <- 1
  df_RCA_results <- aggregate(count ~ ., df_RCA_results, FUN = sum)
  df_RCA_results <- df_RCA_results %>% group_by(seurat_clusters) %>% top_n(n = 5, wt = count)
  df_RCA_results <- df_RCA_results[order(df_RCA_results$seurat_clusters), ]
  
  df_umap_cluster <- FetchData(seurat_PBMC_subset, vars = c("UMAP_1", "UMAP_2", "seurat_clusters"))
  
  # Identify majority vote of RCAv2 annotation for each Seurat cluster (GEX clustering) 
  
  vec_RCA_cluster_annotations <- (df_RCA_results %>% filter(count == max(count)))$RCA.proj.cell
  
  # Check through annotations to verify that per-cluster annotations look reasonable
  
  sink(file = paste(str_name, "RCAv2_Annotation_per_Seurat_cluster.txt", sep = "_"))
  print(str_name)
  for (i in 0:(length(vec_RCA_cluster_annotations)-1)){
    print(df_RCA_results %>% filter(seurat_clusters == i))
  }
  sink()
  
  df_umap_cluster$cell <- NA
  
  for (i in 0:(length(vec_RCA_cluster_annotations)-1)){
    df_umap_cluster[which(df_umap_cluster$seurat_clusters==i), "cell"] <- df_RCAv2_to_annotations[vec_RCA_cluster_annotations[i+1], "Annotation_Names"]
  }
  
  seurat_PBMC_subset$RCA.annotation <- df_umap_cluster[colnames(seurat_PBMC_subset), "cell"]
  
  # Plot UMAPs with RCAv2 annotations, and FeaturePlots of marker genes
  
  pdf(paste0(str_name, "_Library_QC_Plot04_UMAP_Seurat_RCAv2_Annotation_GEX.pdf"))
  print(DimPlot(seurat_PBMC_subset, reduction = "umap", group.by = "RCA.annotation", label = TRUE, label.size = 5) +
          NoLegend() + 
          ggtitle(paste(str_name, " UMAP labelled by RCAv2 annotations", sep = "")))
  dev.off()
  
  pdf(paste0(str_name, "_Library_QC_Plot05_UMAP_Seurat_Marker_Genes.pdf"), width = 30, height = 30)
  print(FeaturePlot(seurat_PBMC_subset, 
                    features = c("CD14", "FCGR3A", "CLEC9A", "CLEC10A", "SIGLEC6", "ITM2C", ### Monocytes, cDC, ASDC, pDC
                                 "MZB1", "MS4A1", "TCL1A", "TNFRSF17", ### Plasma B and B (naive and memory B)
                                 "CD3D", "TRGV9", ### "TRDV2", T and gdT 
                                 "GZMK", "GZMB", ### Central memory and effector / cytotoxic
                                 "CD4", "IL7R", "CCR7", "FOXP3", ### CD4+ T (naive, Treg)
                                 "CD8A", "TRAV1-2", ### CD8+ T, MAIT
                                 "ITGB1", "KLRB1", ### Memory T
                                 "NKG7", "GNLY", "NCAM1", ### NK, including CD56
                                 "PPBP", "HBB", "CD34"))) ### Platelet, RBC, HSPC
  dev.off()
  
  # pK identification (with ground truth "GT")
  
  sweep.res <- paramSweep_v3(seurat_PBMC_subset, PCs = 1:30, sct = FALSE)
  saveRDS(sweep.res, paste0(str_name, "_DoubletFinder.sweep.res.rds"))
  
  gt.calls <- seurat_PBMC_subset@meta.data[rownames(sweep.res[[1]]), "GT"]
  sweep.stats <- summarizeSweep(sweep.res, GT = TRUE, GT.calls = gt.calls)
  
  pdf(paste0(str_name, "_Library_QC_Plot06_DoubletFinder_BCmvn.pdf"), width = 10, height = 10)
  print(bcmvn <- find.pK(sweep.stats))
  print(var_pK <- as.numeric(as.character(bcmvn$pK)[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))
  print(title(paste(str_name, " Plot of BCmetric (blue) and Mean AUC (black) against pK")))
  dev.off()
  
  # Estimate of Homotypic Doublet proportion using RCAv2 space cluster annotations 
  # instead of 5GEX space cluster annotations
  
  annotations <- seurat_rca_PBMC$RCAv2_Cluster_Label
  # annotations <- seurat_PBMC_subset$RCA.annotation
  homotypic.prop <- modelHomotypic(annotations) 
  nExp_poi <- round(number_estimated_doublet_rate*length(colnames(seurat_PBMC_subset)))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # Run DoubletFinder with varying extents of classification stringency (with and without homotypic doublets)
  # Find column names of consecutive DoubletFinder runs
  
  seurat_PBMC_subset <- doubletFinder_v3(seurat_PBMC_subset, PCs = 1:30, pN = 0.25, 
                                         pK = var_pK, nExp = nExp_poi, 
                                         reuse.pANN = FALSE, sct = FALSE)
  
  var_pANN_colname <- colnames(seurat_PBMC_subset[[]])[grep("pANN", colnames(seurat_PBMC_subset[[]]))]
  str_DFclassification_colname <- colnames(seurat_PBMC_subset[[]])[grep("DF.classification", colnames(seurat_PBMC_subset[[]]))]
  
  seurat_PBMC_subset <- doubletFinder_v3(seurat_PBMC_subset, PCs = 1:30, pN = 0.25, 
                                         pK = var_pK, nExp = nExp_poi.adj, 
                                         reuse.pANN = var_pANN_colname, sct = FALSE)
  
  vec_DFclassification_both_colname <- colnames(seurat_PBMC_subset[[]])[grep("DF.classification", colnames(seurat_PBMC_subset[[]]))]
  var_DF_2nditeration <- vec_DFclassification_both_colname[vec_DFclassification_both_colname != str_DFclassification_colname]
  
  ##### Saving and plotting QC metrics, Singlet/Doublet status, metadata, and gene-cell matrices #####
  
  # Plot DoubletFinder singlets and doublets
  
  df_doublet <- FetchData(seurat_PBMC_subset, vars = c("UMAP_1", "UMAP_2", var_DF_2nditeration))
  
  pdf(paste0(str_name, "_Library_QC_Plot07_UMAP_DoubletFinder_Singlets_Doublets.pdf"))
  print(ggplot() + geom_point(data = df_doublet[which(df_doublet[, var_DF_2nditeration] == "Singlet"),],
                              mapping = aes(x = UMAP_1, UMAP_2), size = 0.5, stroke = 0, colour = "blue") + 
          geom_point(data = df_doublet[which(df_doublet[, var_DF_2nditeration] == "Doublet"),],
                     mapping = aes(x = UMAP_1, UMAP_2), size = 0.5, stroke = 0, colour = "red") +
          ggtitle("DoubletFinder singlets (blue) and doublets (red)"))
  dev.off()
  
  vec_final_singlet_barcode <- rownames(df_doublet[which(df_doublet[, var_DF_2nditeration] == "Singlet"), ])
  
  pdf(paste0(str_name, "_Library_QC_Plot08_Density_plot_pANN_all_vs_Singlets.pdf"))
  print(plot(density(seurat_PBMC_subset[[var_pANN_colname]][, 1]), ylim=c(0, 30), main = ""))
  print(lines(density(seurat_PBMC_subset[[var_pANN_colname]][vec_final_singlet_barcode, 1]), col = "red"))
  print(title(paste0(str_name, "\n pANN density of all cells (black) versus DoubletFinder Singlets (red)")))
  dev.off()  
  
  # Obtain vector of singlets identified via both 
  # Demuxlet genetic demultiplexing and DoubletFinder
  
  vec_final_singlet_barcode <- vec_final_singlet_barcode[vec_final_singlet_barcode %in% vec_demuxlet_singlet]
  
  df_doublet$barcode <- rownames(df_doublet)
  df_doublet$final.doublet <- "Doublet"
  df_doublet[which(df_doublet$barcode %in% vec_final_singlet_barcode), "final.doublet"] <- "Singlet"
  
  seurat_PBMC_subset$final.doublet_type <- df_doublet[colnames(seurat_PBMC_subset), "final.doublet"]
  
  pdf(paste0(str_name, "_Library_QC_Plot09_UMAP_Seurat_GeneticDoubletFinder_Singlets_Doublets.pdf"))
  print(ggplot() + geom_point(data = df_doublet[which(df_doublet$final.doublet == "Singlet"), ],
                              mapping = aes(x = UMAP_1, UMAP_2), size = 0.5, stroke = 0, colour = "blue") +
          geom_point(data = df_doublet[which(df_doublet$final.doublet=="Doublet"), ],
                     mapping = aes(x = UMAP_1, UMAP_2), size = 0.5, stroke = 0, colour = "red")  +
          ggtitle("UMAP of genetic+DoubletFinder singlets (blue) and any type of doublets (red)"))
  dev.off()
  
  # Highlight doublets identified by different software tools
  
  df_doublet$demuxlet <- "Doublet"
  df_doublet[which(df_doublet$barcode %in% vec_demuxlet_singlet), "demuxlet"] <- "Singlet"
  
  pdf(paste0(str_name, "_Library_QC_Plot10_UMAP_Seurat_Singlets_Doublets_Demuxlet_DoubletFinder.pdf"))
  print(ggplot() + 
          geom_point(data = df_doublet[which(df_doublet[, var_DF_2nditeration] == "Singlet" &
                                               df_doublet$demuxlet == "Singlet"), ],
                     mapping = aes(x = UMAP_1, UMAP_2), size = 0.5, stroke = 0, colour = "grey") +
          geom_point(data = df_doublet[which(df_doublet[, var_DF_2nditeration] == "Doublet" &
                                               df_doublet$demuxlet == "Singlet"), ],
                     mapping = aes(x = UMAP_1, UMAP_2), size = 0.5, stroke = 0, colour = "red") +
          geom_point(data = df_doublet[which(df_doublet[, var_DF_2nditeration] == "Singlet" &
                                               df_doublet$demuxlet == "Doublet"), ],
                     mapping = aes(x = UMAP_1, UMAP_2), size = 0.5, stroke = 0, colour = "orange") +
          geom_point(data = df_doublet[which(df_doublet[, var_DF_2nditeration] == "Doublet" &
                                               df_doublet$demuxlet == "Doublet"), ],
                     mapping = aes(x = UMAP_1, UMAP_2), size = 0.5, stroke = 0, colour = "blue") +
          ggtitle("Genetic Singlets and Doublet Finder Singlets (grey), 
\nGenetic Singlets and Doublet Finder Doublets (red), 
\nGenetic Doublets and Doublet Finder Singlets (orange), 
\nGenetic Doublets and Doublet Finder Doublets (blue)"))
  dev.off()
  
  # Plot QC metrics for Demuxlet and DoubletFinder singlets and doublets
  
  # DoubletFinder Doublet QC metrics
  
  df_all_NODG_UMI_pMito_DoubletFinder_Doublet <- df_all_NODG_UMI_pMito[df_doublet[which(df_doublet[, var_DF_2nditeration] == "Doublet"), "barcode"], ]
  
  ggplot(df_all_NODG_UMI_pMito_DoubletFinder_Doublet, aes(x = NODG, y = pMito)) + 
    geom_point(size = 0.1) +
    geom_density_2d() + theme_bw(11) + xlim(0, 6000) + ylim(0, 0.1) +
    ggtitle(paste0(str_name, " pMito and NODG \nfor DoubletFinder doublets"))
  ggsave(paste0(str_name, "_Library_QC_Plot11_pMito_NODG_DoubletFinder_Doublets.png"), width = 6, height = 5)
  
  # Demuxlet genetic demultiplexing doublet QC metrics
  
  df_all_NODG_UMI_pMito_DRAGEN_Doublet <- df_all_NODG_UMI_pMito[df_doublet[which(df_doublet$demuxlet == "Doublet"), "barcode"], ]
  
  ggplot(df_all_NODG_UMI_pMito_DRAGEN_Doublet, aes(x = NODG, y = pMito)) + 
    geom_point(size = 0.1) +
    geom_density_2d() + theme_bw(11) + xlim(0, 6000) + ylim(0, 0.1) +
    ggtitle(paste0(str_name, " pMito and NODG \nfor Demuxlet doublets"))
  ggsave(paste0(str_name, "_Library_QC_Plot12_pMito_NODG_Demuxlet_Doublets.png"), width = 6, height = 5)
  
  # QC metrics for Genetic+DoubletFinder singlets and any doublets
  
  df_all_NODG_UMI_pMito$Type <- "Doublet"
  df_all_NODG_UMI_pMito[vec_final_singlet_barcode, "Type"] <- "Singlet"
  
  pdf(paste0(str_name, "_Library_QC_Plot13_pMito_NODG_Singlets_Doublets.pdf"), width = 6, height = 5)
  print(ggplot(df_all_NODG_UMI_pMito, aes(x = NODG, y = pMito, colour = Type)) + 
          geom_point(size = 0.2, stroke = 0) +
          geom_density_2d() + theme_bw(11) + xlim(0, 6000) + ylim(0, 0.1) +
          ggtitle(paste0(str_name, " pMito and NODG for \nGenetic+DoubletFinder singlets and any doublets")))
  dev.off()
  
  # Plot number of Demuxlet+DoubletFinder singlets per donor
  
  df_demuxlet_by_donor <- df_demuxlet[which(df_demuxlet$DROPLET.TYPE == "SNG"), c("BARCODE", "SNG.BEST.GUESS")]
  colnames(df_demuxlet_by_donor) <- c("Barcode", "SampleIdentity")
  rownames(df_demuxlet_by_donor) <- df_demuxlet_by_donor$Barcode
  df_demuxlet_by_donor <- df_demuxlet_by_donor[vec_final_singlet_barcode, ]
  
  df_singlets_per_donor <- as.data.frame(table(df_demuxlet_by_donor$SampleIdentity))
  colnames(df_singlets_per_donor) <- c("Donor", "Freq")
  ggplot(df_singlets_per_donor, aes(x = Donor, y = Freq)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste("Bar chart of Demuxlet donor assignments for \n(Demuxlet and DoubletFinder) singlets in\n", str_name)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(paste0(str_name, "_Library_QC_Plot14_Singlets_by_donors.png"), units = "in", width = 16, height = 6.66)
  
  pdf(paste0(str_name, "_Library_QC_Plot15_UMAP_Seurat_Singlets_Highlighted.pdf"))
  print(DimPlot(seurat_PBMC_subset, reduction = "umap", group.by = "RCA.annotation", 
                cells.highlight = names(which(seurat_PBMC_subset$final.doublet_type == "Singlet")), 
                label = TRUE, label.size = 5) +  
          ggtitle(paste0(str_name, "\n UMAP with (Demuxlet and DoubletFinder) singlets \nhighlighted in red")) + 
          NoLegend())
  dev.off()
  
  pdf(paste0(str_name, "_Library_QC_Plot16_UMAP_Seurat_splitby_SingletsDoublets.pdf"))
  print(DimPlot(seurat_PBMC_subset, reduction = "umap", group.by = "RCA.annotation", 
                split.by = "final.doublet_type", 
                label = TRUE, label.size = 5) +
          ggtitle(paste0(str_name, "\n UMAP split by (Demuxlet and DoubletFinder) singlets \nversus any doublet type")) + 
          NoLegend())
  dev.off()
  
  ### Save data frame with singlet/doublet (Demuxlet, DoubletFinder) information: 
  ### NODG and pMito; cell numbers > NODG 300 (nrow indicative of capture rate); 
  ### Doublet rates (Demuxlet/Genetic, DoubletFinder, either)
  
  df_allcells_SingletDoublet_QC <- merge(df_doublet, df_all_NODG_UMI_pMito, 
                                         by = "row.names")
  df_allcells_SingletDoublet_QC$barcode_name <- paste(df_allcells_SingletDoublet_QC$barcode, "-", str_name, sep = "")
  
  write.table(df_allcells_SingletDoublet_QC, 
              file = paste0(str_name, "_DataFrame_AllCells_SingletDoublet_QC.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  ### Also save a summary data frame of cells, singlets, doublets:
  
  df_summary <- data.frame(Library = str_name, 
                           Number_Cells = nrow(df_allcells_SingletDoublet_QC),
                           Capture_Rate = nrow(df_allcells_SingletDoublet_QC)/40000*100,
                           Number_Both_Singlets = length(which(df_allcells_SingletDoublet_QC$final.doublet == "Singlet")),
                           Number_Both_Doublets = length(which(df_allcells_SingletDoublet_QC$final.doublet == "Doublet")),
                           Percent_DRAGEN_Doublets = length(which(df_allcells_SingletDoublet_QC$demuxlet == "Doublet"))/nrow(df_allcells_SingletDoublet_QC)*100,
                           Percent_DoubletFinder_Doublets = length(which(df_allcells_SingletDoublet_QC[, var_DF_2nditeration] == "Doublet"))/nrow(df_allcells_SingletDoublet_QC)*100,
                           Percent_Both_Doublets = length(which(df_allcells_SingletDoublet_QC$final.doublet == "Doublet"))/nrow(df_allcells_SingletDoublet_QC)*100
  )
  
  write.table(df_summary, 
              file = paste0(str_name, "_DataFrame_Summary_CellsSingletsDoublets.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  ### Save data frame of metadata from seurat_PBMC_subset
  ### seurat_PBMC_subset[[]]: Singlet/doublet (Demuxlet, DoubletFinder) information: 
  ### NODG and pMito; cell numbers > NODG 300 (nrow indicative of capture rate); 
  ### Doublet rates (Demuxlet/Genetic, DoubletFinder, either); RCAv2 information
  
  rownames(df_demuxlet) <- df_demuxlet$BARCODE
  df_demuxlet_singlet <- df_demuxlet[vec_final_singlet_barcode, ]
  df_demuxlet_singlet$Genotyping_ID <- df_demuxlet_singlet$SNG.BEST.GUESS
  seurat_PBMC_subset_singlet <- subset(x = seurat_PBMC_subset, 
                                       subset = final.doublet_type == "Singlet")
  seurat_PBMC_subset_singlet$Genotyping_ID <- df_demuxlet_singlet[colnames(seurat_PBMC_subset_singlet), "Genotyping_ID"]
  seurat_PBMC_subset_singlet$pMito <- df_all_NODG_UMI_pMito[colnames(seurat_PBMC_subset_singlet), "pMito"]
  seurat_PBMC_subset_singlet$NODG <- df_all_NODG_UMI_pMito[colnames(seurat_PBMC_subset_singlet), "NODG"]
  seurat_PBMC_subset_singlet$nUMI <- df_all_NODG_UMI_pMito[colnames(seurat_PBMC_subset_singlet), "nUMI"]
  seurat_PBMC_subset_singlet$barcode_only <- colnames(seurat_PBMC_subset_singlet)
  write.table(seurat_PBMC_subset_singlet[[]], 
              file = paste0(str_name, "_DataFrame_Singlets_SingletDoublet_QC_Annotations.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  ### Save data frames of names of genes that are expressed by >= 0.1% of cells in this library
  
  df_genes <- data.frame(Genes = row.names(seurat_PBMC_subset_singlet))
  write.table(df_genes, 
              file = paste0(str_name, "_DataFrame_Genes_passing_filter.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  ### Code for checking sex chromosome gene expression for data-metadata concordance
  
  ### Count chrY (non-PAR: likely male sex; PAR: both sexes) UMIs per cell
  
  vec_chrY_nonPAR_genes <- Matrix::colSums(matrix_PBMC[intersect(df_chrY_nonPAR_genes$V1, rownames(matrix_PBMC)), , drop = FALSE])
  vec_chrY_PAR_genes <- Matrix::colSums(matrix_PBMC[intersect(df_chrY_PAR_genes$V1, rownames(matrix_PBMC)), , drop = FALSE])
  
  ### Add chrY (non-PAR and PAR) UMIs per cell to Seurat object
  ### and aggregate these values per Genotyping_ID
  ### for elucidating scRNA-seq-inferred biological sex
  
  seurat_PBMC_subset_singlet$chrY_nonPAR_UMIs <- vec_chrY_nonPAR_genes[colnames(seurat_PBMC_subset_singlet)]
  seurat_PBMC_subset_singlet$chrY_PAR_UMIs <- vec_chrY_PAR_genes[colnames(seurat_PBMC_subset_singlet)]
  df_GenotypingID_chrY_UMIs <- seurat_PBMC_subset_singlet[[c("Genotyping_ID", "chrY_nonPAR_UMIs", "chrY_PAR_UMIs")]]
  
  df_merged_chrY_UMIs <- merge(aggregate(chrY_nonPAR_UMIs ~ Genotyping_ID, df_GenotypingID_chrY_UMIs, FUN = sum), 
                               aggregate(chrY_PAR_UMIs ~ Genotyping_ID, df_GenotypingID_chrY_UMIs, FUN = sum),
                               by = "Genotyping_ID")
  df_merged_chrY_UMIs$nonPAR_to_PAR <- as.numeric(df_merged_chrY_UMIs$chrY_nonPAR_UMIs/df_merged_chrY_UMIs$chrY_PAR_UMIs)
  df_merged_chrY_UMIs$scRNAseq_sex <- NA
  
  for (i in 1:length(df_merged_chrY_UMIs$Genotyping_ID)){
    if (df_merged_chrY_UMIs[i, "nonPAR_to_PAR"] >= 0.5){
      df_merged_chrY_UMIs[i, "scRNAseq_sex"] <- "Male"
    }
    if (df_merged_chrY_UMIs[i, "nonPAR_to_PAR"] < 0.5) {
      df_merged_chrY_UMIs[i, "scRNAseq_sex"] <- "Female"
    }
  }
  
  ### Save scRNA-seq inferred sex and numbers of singlets per donor
  ### to check for data-metadata concordance, successful genetic demultiplexing, 
  ### and to guard against possible sample swaps.
  ### Should have >100s of cells per donor.
  
  df_singlets_per_donor$Genotyping_ID <- as.character(df_singlets_per_donor$Donor)
  df_singlets_per_donor <- df_singlets_per_donor[, c("Genotyping_ID", "Freq")]
  colnames(df_singlets_per_donor) <- c("Genotyping_ID", "Singlets")
  df_merged_donor_scRNAseq_metadata <- merge(df_merged_chrY_UMIs, 
                                             df_singlets_per_donor,
                                             by = "Genotyping_ID")
  
  write.table(df_merged_donor_scRNAseq_metadata, 
              file = paste0(str_name, "_DataFrame_Singlets_and_scRNAseq_sex_by_GenotypingID.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  ### Code for identifying cells with more than 10 (HBA1 and HBB) (RBC) UMIs
  ### for exclusion from dataset and for writing these cells to an output file
  
  vec_HBA1_HBB <- Matrix::colSums(matrix_PBMC[c("HBA1", "HBB"), , drop = FALSE])
  vec_cells_to_remove <- names(vec_HBA1_HBB[vec_HBA1_HBB > 10])
  write.table(data.frame(Library = rep(str_name, length(vec_cells_to_remove)),
                         Cells_with_RBC_UMIs = vec_cells_to_remove, 
                         RBC_UMIs = vec_HBA1_HBB[vec_HBA1_HBB > 10]), 
              file = paste0(str_name, "_DataFrame_RBC_UMI_cells.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  ### Code for excluding cells corresponding to non-consented genotypes
  
  # Genotyping_IDs found in the current library
  vec_unique_genotypes <- unique(df_demuxlet_singlet$Genotyping_ID)
  vec_str_genotype <- c()
  for (str_genotype in vec_unique_genotypes){
    vec_str_genotype <- c(vec_str_genotype, 
                          paste(strsplit(str_genotype, split = "_")[[1]][-1], 
                                collapse = "_"))
  }
  
  # List of Genotyping_IDs consented for AIDA Phase 1
  vec_consented_genotypes <- df_consented_genotypes$Genotyping_ID
  
  # Get genotype IDs (including prefix) of samples to exclude from this library
  vec_genotypes_to_exclude <- vec_str_genotype[!(vec_str_genotype %in% vec_consented_genotypes)] ### Genotyping_ID to exclude
  vec_match_genotypes_to_exclude <- match(vec_genotypes_to_exclude, vec_str_genotype)
  vec_withprefix_genotypes_to_exclude <- vec_unique_genotypes[vec_match_genotypes_to_exclude]
  
  if (length(vec_genotypes_to_exclude) > 0){
    write.table(data.frame(Genotypes_to_exclude = vec_genotypes_to_exclude, 
                           Genotypes_to_exclude_with_prefix = vec_withprefix_genotypes_to_exclude),
                file = paste0(str_name, "_DataFrame_Genotypes_to_exclude.txt"),
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    for (str_prefix in vec_withprefix_genotypes_to_exclude){
      vec_cells_to_remove <- c(vec_cells_to_remove,
                               row.names(df_demuxlet_singlet[which(df_demuxlet_singlet$Genotyping_ID == str_prefix), ]))
    }
    vec_cells_to_remove <- unique(vec_cells_to_remove)
  }
  
  # Remove RBC UMI cells and cells from non-consented genotypes from vec_final_singlet_barcode
  
  vec_final_singlet_barcode <- vec_final_singlet_barcode[!(vec_final_singlet_barcode %in% vec_cells_to_remove)]
  
  # Save RCA object with Genetic+DoubletFinder singlets
  # after removing RBC UMI cells and cells corresponding to non-consented genotypes
  # All genes retained.
  
  rca_PBMC.new <- createRCAObject(rawData = seurat_PBMC@assays$RNA@counts[, vec_final_singlet_barcode])
  
  # Normalise data
  
  rca_PBMC.new <- dataLogNormalise(rca_PBMC.new)
  
  # Obtain RCAv2 projection results from earlier in this piece of code.
  # Assign combined RCAv2 projection result to RCA object.
  
  rca_PBMC.new$projection.data <- df_all_projection[, vec_final_singlet_barcode]
  
  # Rename RCA object barcodes to include library name in barcodes
  
  vec_RCA_barcode_name <- paste(vec_final_singlet_barcode, "-", str_name, sep = "")
  colnames(rca_PBMC.new$raw.data) <- vec_RCA_barcode_name
  colnames(rca_PBMC.new$data) <- vec_RCA_barcode_name
  colnames(rca_PBMC.new$projection.data) <- vec_RCA_barcode_name
  
  # Estimate the most probable cell type label for each cell and save RCA object
  
  rca_PBMC.new <- estimateCellTypeFromProjection(rca_PBMC.new, confidence = NULL)
  saveRDS(rca_PBMC.new, paste0(str_name, "_RCA_after_doublet_and_RBC_removal.rds"))
  
  # Save Seurat and RCA objects with all Demuxlet and DoubletFinder singlets
  # and all genes from library from consented genotypes
  # Before saving: rename cells to include country+institute+batch+library name
  
  seurat_PBMC_singlets <- seurat_PBMC[, vec_final_singlet_barcode]
  seurat_PBMC_singlets <- RenameCells(seurat_PBMC_singlets, new.names = seurat_PBMC_singlets$barcode_name)
  saveRDS(seurat_PBMC_singlets, file = paste0(str_name, "_Seurat_Consented_Singlets_AllGenes.RDS"))
  
  # Save Seurat object with all cells (NODG >= 300) and all genes from library
  # Before saving: rename cells to include country+institute+batch+library name
  
  seurat_PBMC <- RenameCells(seurat_PBMC, new.names = seurat_PBMC$barcode_name)
  saveRDS(seurat_PBMC, file = paste0(str_name, "_Seurat_AllCells_AllGenes.RDS"))
  
}
