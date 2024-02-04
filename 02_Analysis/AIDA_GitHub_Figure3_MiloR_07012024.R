### AIDA Phase 1 Data Freeze v2 Cell neighbourhood enrichment MiloR analysis: version 7 January 2024

##### Running miloR for cell neighbourhood enrichment analyses #####

### Installing MiloR and related packages
# install.packages("BiocManager")
# install.packages('Seurat')
# install.packages("tidyverse")
# BiocManager::install("miloR")
# BiocManager::install("scater")
# install.packages("statmod")

##### Preparing environment for analyses #####

library(miloR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(patchwork)

##### INPUT REQUIRED: READ IN CELL TYPE ANNOTATION HIERARCHY AS df_master_hierarchy #####
##### INPUT REQUIRED: READ IN AIDA METADATA AS df_all_metadata #####

df_master_hierarchy <- read.table("AIDA_Annotation_Subclustering_Final_AllHierarchy.txt", 
                                  header = TRUE, sep = "\t")
rownames(df_master_hierarchy) <- df_master_hierarchy$Annotation_Final

str_k <- 900
vec_majorcelltypes <- c("B", "Myeloid_pDC", "TNK_CD4", "TNK_nonCD4")
  
##### Consider only donors with at least 800 cells and country-ancestry combinations with at least 50 donors per combination #####

df_all_metadata <- read.table("df_metadata.txt", 
                              header = TRUE, sep = "\t")
df_all_metadata$count <- 1
rownames(df_all_metadata) <- df_all_metadata$barcode_name

df_total_cell_count <- df_all_metadata[, c("DCP_ID", "count")]
df_total_cell_count <- aggregate(count ~ ., df_total_cell_count, FUN = sum)
rownames(df_total_cell_count) <- df_total_cell_count$DCP_ID
vec_donors_to_exclude <- df_total_cell_count$DCP_ID[which(df_total_cell_count$count < 800)]
df_all_metadata <- df_all_metadata[!(df_all_metadata$DCP_ID %in% vec_donors_to_exclude), ]
df_all_metadata <- df_all_metadata[!(df_all_metadata$Country == "IN"), ]
df_all_metadata <- df_all_metadata[!(df_all_metadata$Ancestry == "European"), ]

##### Reading in cell subtype data, and including metadata and annotations (human diversity, cell type annotation) #####

for (str_Ancestry in unique(df_all_metadata$Ancestry)){
  df_all_metadata[ , str_Ancestry] <- 0
  df_all_metadata[ , str_Ancestry][which(df_all_metadata$Ancestry == str_Ancestry)] <- 1
}

for (str_celltype in c("Myeloid_pDC", "TNK_CD4", "TNK_nonCD4")){

seurat_AIDA_subset <- readRDS(paste0(str_celltype, "1.subclustering_reintegrated.rds"))
  
seurat_AIDA_subset[["Annotation_Final"]] <- NA
seurat_AIDA_subset[["Annotation_Level1"]] <- NA
seurat_AIDA_subset[["Annotation_Level2"]] <- NA
seurat_AIDA_subset[["Annotation_Level3"]] <- NA
seurat_AIDA_subset[["Annotation_Level4"]] <- NA
  
df_cluster_hierarchy <- read.table(paste0("AIDA_Annotation_Subclustering_Final_", str_celltype, ".txt"), 
                                     header = TRUE, sep = "\t")
rownames(df_cluster_hierarchy) <- df_cluster_hierarchy$Cluster
  
seurat_AIDA_subset$Annotation_Final <- df_cluster_hierarchy[as.character(seurat_AIDA_subset$seurat_clusters), "Annotation_Final"]
  
for (str_level_number in 1:4){
  seurat_AIDA_subset[[paste0("Annotation_Level", str_level_number)]] <- 
    df_master_hierarchy[seurat_AIDA_subset$Annotation_Final, paste0("Annotation_Level", str_level_number)]
}

vec_barcodes_to_keep <- intersect(rownames(df_all_metadata), colnames(seurat_AIDA_subset))  
df_barcodes_to_keep <- df_all_metadata[vec_barcodes_to_keep, ]
seurat_AIDA_subset <- seurat_AIDA_subset[, vec_barcodes_to_keep]

for (str_Ancestry in unique(seurat_AIDA_subset$Ancestry)){
  seurat_AIDA_subset[[str_Ancestry]] <- 0
  seurat_AIDA_subset[[str_Ancestry]] <- df_barcodes_to_keep[colnames(seurat_AIDA_subset), str_Ancestry]
}
  
seurat_AIDA_subset[["Batch"]] <- sapply(seurat_AIDA_subset$Library, function(x){strsplit(x, split = "_L")[[1]][1]})

seurat_AIDA_subset$Ages_19_to_32 <- 0
seurat_AIDA_subset$Ages_33_to_40 <- 0
seurat_AIDA_subset$Ages_41_to_49 <- 0
seurat_AIDA_subset$Ages_50_to_77 <- 0
seurat_AIDA_subset$Age <- as.numeric(seurat_AIDA_subset$Age)

seurat_AIDA_subset$Ages_19_to_32[which(seurat_AIDA_subset$Age >= 19 & seurat_AIDA_subset$Age < 33)] <- 1
seurat_AIDA_subset$Ages_33_to_40[which(seurat_AIDA_subset$Age >= 33 & seurat_AIDA_subset$Age < 41)] <- 1
seurat_AIDA_subset$Ages_41_to_49[which(seurat_AIDA_subset$Age >= 41 & seurat_AIDA_subset$Age < 50)] <- 1
seurat_AIDA_subset$Ages_50_to_77[which(seurat_AIDA_subset$Age >= 50 & seurat_AIDA_subset$Age <= 77)] <- 1

##### Running MiloR to assess influence of covariates on cell neighbourhood abundance ######

##### Preparing Milo object #####

DefaultAssay(seurat_AIDA_subset) <- "integrated"
Idents(object = seurat_AIDA_subset) <- "Annotation_Level3"

sce_AIDA_subset <- as.SingleCellExperiment(seurat_AIDA_subset)
colData(sce_AIDA_subset) <- DataFrame(seurat_AIDA_subset[[]])
# identical(DataFrame(seurat_AIDA_subset[[]]), colData(sce_AIDA_subset))
# [1] TRUE
milo_AIDA_subset <- Milo(sce_AIDA_subset)
# plotUMAP(milo_AIDA_subset)

# To save memory
rm(sce_AIDA_subset)

### Want average neighbourhood size of > 5N, where N is number of samples
### i.e. want average neighbourhood size of > 2810
### Trying a few k and identify peak neighbourhood size from plotNhoodSizeHist:
### Try k = 900 - gives peak Nhood size of ~3100

##### Performing MiloR neighbourhood analysis #####

milo_AIDA_subset <- buildGraph(milo_AIDA_subset, k = str_k, d = 30)

### For MiloR v1.5, used refinement_scheme = "graph" in makeNhoods function, 
### to avoid memory + computational time issues associated with miloR::calcNhoodDistance
### This takes ~1 hour to run, max memory ~25 GB, for 70k cells

milo_AIDA_subset <- makeNhoods(milo_AIDA_subset, prop = 0.1, k = str_k, d = 30, 
                               refined = TRUE, refinement_scheme = "graph")

plot_NhoodSizeHist <- plotNhoodSizeHist(milo_AIDA_subset)
ggsave(filename = paste("Rplot_MiloR_NhoodSizeHist", str_celltype, "k", str_k, ".pdf", sep = "_"), 
       plot = plot_NhoodSizeHist, width = 26.68, height = 11.24, units = "in")

milo_AIDA_subset <- countCells(milo_AIDA_subset,
                               meta.data = data.frame(DataFrame(seurat_AIDA_subset[[]])), 
                               sample = "DCP_ID")
# head(nhoodCounts(milo_AIDA_subset))

saveRDS(milo_AIDA_subset, file = paste("milo_AIDA_subset", str_celltype, "k", str_k, ".RDS", sep = "_"))

##### Design matrix for covariates of interest: sex, age, and ancestry, option of having batch covariate #####
traj_design <- 
  data.frame(colData(milo_AIDA_subset))[, c("DCP_ID", "Batch", "Sex",
                                            "Age", "Ages_19_to_32", "Ages_50_to_77",
                                            "Ancestry", unique(df_barcodes_to_keep$Ancestry))]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$DCP_ID

### Note that introducing batch as covariate alongside Ancestry would result in collinearity in covariates
### Will get a "Design matrix not of full rank" issue

### For MiloR v1.5, used fdr.weighting = "graph-overlap" in testNhoods function 
### to avoid memory + computational time issues

  ##### Running MiloR to assess influence of covariates on cell neighbourhood abundance ######
  
  ##### Preparing Milo object #####
  
  DefaultAssay(seurat_AIDA_subset) <- "integrated"
  Idents(object = seurat_AIDA_subset) <- "Annotation_Level3"
  
  milo_AIDA_subset <- readRDS(paste0(str_celltype, "/", paste("milo_AIDA_subset", str_celltype, "k", str_k, ".RDS", sep = "_")))
  
  str_covariate_comparison <- "Sex"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Ancestry + Age + Sex, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  
  ### Function for plotting beeswarm plot by cell subtypes
  
  function_beeswarm <- function(){
    da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level3")
    da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level4")
    
    plotDAbeeswarm(da_results, group.by = "Annotation_Level3") +
      ggtitle(paste0("Cell neighbourhood enrichment of ", str_celltype, " in ", str_covariate_comparison))
    ggsave(paste("Rplot_MiloR_Beeswarmby_Annotation_Level3", str_covariate_comparison, 
                 str_celltype, "k", str_k, ".png", sep = "_"), 
           width = 26.68, height = 11.24, units = "in")
    
    plotDAbeeswarm(da_results, group.by = "Annotation_Level4") +
      ggtitle(paste0("Cell neighbourhood enrichment of ", str_celltype, " in ", str_covariate_comparison))
    ggsave(paste("Rplot_MiloR_Beeswarmby_Annotation_Level4", str_covariate_comparison, 
                 str_celltype, "k", str_k, ".png", sep = "_"), 
           width = 26.68, height = 11.24, units = "in")
  }
  
  ### Function for plotting log2(mean fold-change of MiloR 2^log2FC) values on UMAP
  
  function_plot_log2_meanfoldchange <- function(){
    ### Order of Nhood in da_results is numerical ascending
    ### Converting log2 fold-change results to fold-change
    vec_FC <- 2^(da_results$logFC)
    matrix_cell_neighbourhood <- nhoods(milo_AIDA_subset)
    matrix_cell_neighbourhood_foldchange <- sweep(matrix_cell_neighbourhood, 
                                                  MARGIN = 2, 
                                                  STATS = vec_FC, 
                                                  FUN = `*`)
    rm(matrix_cell_neighbourhood)
    vec_log2_meanfoldchange <- log2(rowSums(matrix_cell_neighbourhood_foldchange) / 
                                      rowSums(matrix_cell_neighbourhood_foldchange > 0))
    df_UMAP <- FetchData(seurat_AIDA_subset, vars = c("UMAP_1","UMAP_2"))
    df_UMAP$log2_MeanFoldChange <- vec_log2_meanfoldchange[row.names(df_UMAP)]
    
    ### Capping plotted log2(mean fold change) at |2| 
    df_UMAP[(df_UMAP$log2_MeanFoldChange > 2),"log2_MeanFoldChange"] <- 2
    df_UMAP[(df_UMAP$log2_MeanFoldChange < -2),"log2_MeanFoldChange"] <- -2
    
    ggplot()+
      geom_point(df_UMAP, 
                 mapping = aes(x = UMAP_1, y = UMAP_2, 
                               colour = log2_MeanFoldChange), 
                 size = 0.5, stroke = 0.1) +
      scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                             values = c(1, 0.54, 0.5, 0.46, 0),  
                             limits = c(-2, 2)) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) + 
      ggtitle(paste0("Enrichment of ", str_celltype, " cells for ", str_covariate_comparison))
    ggsave(paste("RPlot_Figure_MiloR_UMAP_log2_MeanFoldChange", 
                 str_celltype, "k", str_k, str_covariate_comparison, 
                 ".png", sep = "_"), device = "png")
  }
  
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  ### Sanity checking output through plotting Male and Female enrichment
  
  traj_design$Male <- ifelse(traj_design$Sex == 'Male', 1, 0)
  colData(milo_AIDA_subset)$Male <- ifelse(colData(milo_AIDA_subset)$Sex == 'Male', 1, 0)
  traj_design$Female <- ifelse(traj_design$Sex == 'Female', 1, 0)
  colData(milo_AIDA_subset)$Female <- ifelse(colData(milo_AIDA_subset)$Sex == 'Female', 1, 0)
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Age + Ancestry + Male, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  str_covariate_comparison <- "Male"
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Age + Ancestry + Female, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  str_covariate_comparison <- "Female"
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  ##### Testing for influence of Ancestry on cell neighbourhood abundance: SG_Chinese, SG_Malay, SG_Indian, Japanese, Korean, Thai #####
  ### For MiloR v1.5, used fdr.weighting = "graph-overlap" in testNhoods function 
  ### to avoid memory + computational time issues
  
  str_covariate_comparison <- "SG_Chinese"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Age + Sex + SG_Chinese, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  str_covariate_comparison <- "SG_Malay"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Age + Sex + SG_Malay, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  str_covariate_comparison <- "SG_Indian"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Age + Sex + SG_Indian, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  str_covariate_comparison <- "Japanese"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Age + Sex + Japanese, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  str_covariate_comparison <- "Korean"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Age + Sex + Korean, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  str_covariate_comparison <- "Thai"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Age + Sex + Thai, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
  ##### Testing for influence of Age on cell neighbourhood abundance: Ages_50_to_77 for full workflow, Ages_19_to_32 for log2_meanfoldchange #####
  ### For MiloR v1.5, used fdr.weighting = "graph-overlap" in testNhoods function 
  ### to avoid memory + computational time issues
  
  ### Run into issue with beeswarm plot for Ages_19_to_32, so just plot log2(meanfoldchange) for Ages_19_to_32
  
  str_covariate_comparison <- "Ages_19_to_32"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Ancestry + Sex + Ages_19_to_32, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  function_plot_log2_meanfoldchange()
  
  str_covariate_comparison <- "Ages_50_to_77"
  
  da_results <- testNhoods(milo_AIDA_subset, 
                           design = ~ Ancestry + Sex + Ages_50_to_77, 
                           design.df = traj_design, fdr.weighting = "graph-overlap")
  function_beeswarm()
  function_plot_log2_meanfoldchange()
  
}

##### Figure panel specific plots #####

### Male, NK+CD8:

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Ancestry + Male, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")
str_covariate_comparison <- "Male"

da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level4")

vec_cells_to_remove <- c(which(da_results$Annotation_Level4 == "CD8+_T"), 
                         which(da_results$Annotation_Level4 == "flagged_CD16+_NK_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_gdT_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_naive_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_NK_low_exp"), 
                         which(da_results$Annotation_Level4 == "flagged_T_CD8+_T_gdT_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_T_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "T_IFNhi"))

plotDAbeeswarm(da_results[-vec_cells_to_remove, ], group.by = "Annotation_Level4") +
  ggtitle(paste0("Cell neighbourhood enrichment of CD8+ T, gdT, ILC, and NK cells in donors ", str_covariate_comparison))

ggsave(paste("Rplot_Figure_MiloR_Beeswarmby_Annotation_Level4", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"), 
       width = 13.34, height = 11.24, units = "in")

ggplot()+
  geom_point(df_UMAP, 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0),  
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of CD8+ T, gdT, ILC, and NK cells in ", str_covariate_comparison, " donors"))

### 50 and above, NK+CD8:

str_covariate_comparison <- "Ages_50_to_77"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Ancestry + Sex + Ages_50_to_77, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

ggplot()+
  geom_point(df_UMAP, 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0), 
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of CD8+ T, gdT, ILC, and NK cells in donors of ", str_covariate_comparison))

da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level4")

vec_cells_to_remove <- c(which(da_results$Annotation_Level4 == "CD8+_T"), 
                         which(da_results$Annotation_Level4 == "flagged_CD16+_NK_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_gdT_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_naive_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_NK_low_exp"), 
                         which(da_results$Annotation_Level4 == "flagged_T_CD8+_T_gdT_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_T_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "T_IFNhi"))

plotDAbeeswarm(da_results[-vec_cells_to_remove, ], group.by = "Annotation_Level4") +
  ggtitle(paste0("Cell neighbourhood enrichment of CD8+ T, gdT, ILC, and NK cells in donors ", str_covariate_comparison))

ggsave(paste("Rplot_Figure_MiloR_Beeswarmby_Annotation_Level4", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"), 
       width = 13.34, height = 11.24, units = "in")

### SG_Malay, NK+CD8:

str_covariate_comparison <- "SG_Malay"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Sex + SG_Malay, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level4")

vec_cells_to_remove <- c(which(da_results$Annotation_Level4 == "CD8+_T"), 
                         which(da_results$Annotation_Level4 == "flagged_CD16+_NK_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_gdT_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_naive_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_NK_low_exp"), 
                         which(da_results$Annotation_Level4 == "flagged_T_CD8+_T_gdT_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_T_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "T_IFNhi"))

plotDAbeeswarm(da_results[-vec_cells_to_remove, ], group.by = "Annotation_Level4") +
  ggtitle(paste0("Cell neighbourhood enrichment of CD8+ T, gdT, ILC, and NK cells in ", str_covariate_comparison, " donors"))

ggsave(paste("Rplot_Figure_MiloR_Beeswarmby_Annotation_Level4", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"), 
       width = 13.34, height = 11.24, units = "in")

### Run log2(mean_fold_change) plotting code

ggplot()+
  geom_point(df_UMAP, 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0),  
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of CD8+ T, gdT, ILC, and NK cells in ", str_covariate_comparison, " donors"))

### 50 and above, CD4:

str_covariate_comparison <- "Ages_50_to_77"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Ancestry + Sex + Ages_50_to_77, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

vec_cells_to_keep <- c(names(which(seurat_AIDA_subset$Annotation_Level2 == "CD4+_T")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "dnT")))

ggplot()+
  geom_point(df_UMAP[vec_cells_to_keep, ], 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0), 
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of CD4+ T cells in donors of ", str_covariate_comparison))

### Plotting DotPlot of DEG marking out least enriched neighbourhood

seurat_AIDA_subset_corecelltypes <- seurat_AIDA_subset[, vec_cells_to_keep]

### UMAP of least enriched cell neighbourhood

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "integrated"

seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "CD4+_T"))] <- "CD4+_T_unknown"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "CD4+_T_memory"))] <- "CD4+_T_memory"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "flagged_CD4+_T_low_exp"))] <- "CD4+_T_unknown"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "flagged_CD4+_T_platelet_gene"))] <- "CD4+_T_unknown"

Idents(object = seurat_AIDA_subset_corecelltypes) <- "Annotation_Level3"

int_least_enriched <- which(da_results$logFC == min(da_results$logFC))
vec_cells_interest <- names(which(nhoods(milo_AIDA_subset)[, int_least_enriched] > 0))
vec_cells_interest <- vec_cells_interest[vec_cells_interest %in% colnames(seurat_AIDA_subset_corecelltypes)]

DimPlot(seurat_AIDA_subset_corecelltypes, reduction = "umap", 
        cells.highlight = vec_cells_interest, sizes.highlight = 0.8, 
        label = TRUE, raster = FALSE, repel = TRUE, label.size = 8) + NoLegend()

ggsave(paste("Rplot_MiloR_UMAP_least_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".png", sep = "_"))

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "RNA"
da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level3")

### Identifying cell subtype implicated in least enriched neighbourhood
vec_celltype <- names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level3 == da_results$Annotation_Level3[int_least_enriched]))
vec_other <- vec_celltype [! vec_celltype %in% vec_cells_interest]
df_markergenes <- FindMarkers(seurat_AIDA_subset_corecelltypes, ident.1 = vec_cells_interest, ident.2 = vec_other)

vec_high_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = avg_log2FC)
vec_low_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = -avg_log2FC)

### Identifying annotated cell type of least enriched neighbourhood
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_cells_interest) <- paste0("Most depleted ", da_results$Annotation_Level3[int_least_enriched])
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_other) <- paste0("Other ", da_results$Annotation_Level3[int_least_enriched])
DotPlot(seurat_AIDA_subset_corecelltypes, features = c(rownames(vec_high_genes), rownames(vec_low_genes))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste("Rplot_MiloR_DotPlot_least_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"))

da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level4")

vec_cells_to_remove <- c(which(da_results$Annotation_Level4 == "CD4+_T"), 
                         which(da_results$Annotation_Level4 == "CD8+_T"), 
                         which(da_results$Annotation_Level4 == "CD8+_T_GZMBhi"), 
                         which(da_results$Annotation_Level4 == "CD8+_T_GZMKhi"), 
                         which(da_results$Annotation_Level4 == "CD8+_T_naive"), 
                         which(da_results$Annotation_Level4 == "flagged_CD4+_T_low_exp"), 
                         which(da_results$Annotation_Level4 == "flagged_CD4+_T_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_naive_CCR7lo"), 
                         which(da_results$Annotation_Level4 == "flagged_T_monocyte_NK_gene"))

plotDAbeeswarm(da_results[-vec_cells_to_remove, ], group.by = "Annotation_Level4") +
  ggtitle(paste0("Cell neighbourhood enrichment of CD4+ T cells in donors ", str_covariate_comparison))

ggsave(paste("Rplot_Figure_MiloR_Beeswarmby_Annotation_Level4", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"), 
       width = 13.34, height = 11.24, units = "in")

### Female, CD4:

str_covariate_comparison <- "Female"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Ancestry + Female, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level4")

vec_cells_to_remove <- c(which(da_results$Annotation_Level4 == "CD4+_T"), 
                         which(da_results$Annotation_Level4 == "CD8+_T"), 
                         which(da_results$Annotation_Level4 == "CD8+_T_GZMBhi"), 
                         which(da_results$Annotation_Level4 == "CD8+_T_GZMKhi"), 
                         which(da_results$Annotation_Level4 == "CD8+_T_naive"), 
                         which(da_results$Annotation_Level4 == "flagged_CD4+_T_low_exp"), 
                         which(da_results$Annotation_Level4 == "flagged_CD4+_T_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_naive_CCR7lo"), 
                         which(da_results$Annotation_Level4 == "flagged_T_monocyte_NK_gene"))

plotDAbeeswarm(da_results[-vec_cells_to_remove, ], group.by = "Annotation_Level4") +
  ggtitle(paste0("Cell neighbourhood enrichment of CD4+ T cells in ", str_covariate_comparison, " donors"))

ggsave(paste("Rplot_Figure_MiloR_Beeswarmby_Annotation_Level4", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"), 
       width = 13.34, height = 11.24, units = "in")

### SG_Malay, CD4:

str_covariate_comparison <- "SG_Malay"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Sex + SG_Malay, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

vec_cells_to_keep <- c(names(which(seurat_AIDA_subset$Annotation_Level2 == "CD4+_T")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "dnT")))

ggplot()+
  geom_point(df_UMAP[vec_cells_to_keep, ], 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0), 
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of CD4+ T cells in ", str_covariate_comparison, " donors"))

### Plotting DotPlot of DEG marking out most enriched neighbourhood

seurat_AIDA_subset_corecelltypes <- seurat_AIDA_subset[, vec_cells_to_keep]

### UMAP of most enriched cell neighbourhood

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "integrated"

seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "CD4+_T"))] <- "CD4+_T_unknown"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "CD4+_T_memory"))] <- "CD4+_T_memory"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "flagged_CD4+_T_low_exp"))] <- "CD4+_T_unknown"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "flagged_CD4+_T_platelet_gene"))] <- "CD4+_T_unknown"

Idents(object = seurat_AIDA_subset_corecelltypes) <- "Annotation_Level3"

int_most_enriched <- which(da_results$logFC == max(da_results$logFC))
vec_cells_interest <- names(which(nhoods(milo_AIDA_subset)[, int_most_enriched] > 0))
vec_cells_interest <- vec_cells_interest[vec_cells_interest %in% colnames(seurat_AIDA_subset_corecelltypes)]

DimPlot(seurat_AIDA_subset_corecelltypes, reduction = "umap", 
        cells.highlight = vec_cells_interest, sizes.highlight = 0.8, 
        label = TRUE, raster = FALSE, repel = TRUE, label.size = 8) + NoLegend()

ggsave(paste("Rplot_MiloR_UMAP_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".png", sep = "_"))

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "RNA"
da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level3")

### Identifying cell subtype implicated in most enriched neighbourhood
vec_celltype <- names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level3 == da_results$Annotation_Level3[int_most_enriched]))
vec_other <- vec_celltype [! vec_celltype %in% vec_cells_interest]
df_markergenes <- FindMarkers(seurat_AIDA_subset_corecelltypes, ident.1 = vec_cells_interest, ident.2 = vec_other)

vec_high_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = avg_log2FC)
vec_low_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = -avg_log2FC)

### Identifying annotated cell type of most enriched neighbourhood
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_cells_interest) <- paste0("Most enriched ", da_results$Annotation_Level3[int_most_enriched])
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_other) <- paste0("Other ", da_results$Annotation_Level3[int_most_enriched])
DotPlot(seurat_AIDA_subset_corecelltypes, features = c(rownames(vec_high_genes), rownames(vec_low_genes))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste("Rplot_MiloR_DotPlot_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"))

### SG_Indian, CD4:

str_covariate_comparison <- "SG_Indian"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Sex + SG_Indian, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

vec_cells_to_keep <- c(names(which(seurat_AIDA_subset$Annotation_Level2 == "CD4+_T")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "dnT")))

ggplot()+
  geom_point(df_UMAP[vec_cells_to_keep, ], 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0), 
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of CD4+ T cells in ", str_covariate_comparison, " donors"))

### Plotting DotPlot of DEG marking out most enriched neighbourhood

seurat_AIDA_subset_corecelltypes <- seurat_AIDA_subset[, vec_cells_to_keep]

### UMAP of most enriched cell neighbourhood

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "integrated"

seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "CD4+_T"))] <- "CD4+_T_unknown"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "CD4+_T_memory"))] <- "CD4+_T_memory"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "flagged_CD4+_T_low_exp"))] <- "CD4+_T_unknown"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level4 == "flagged_CD4+_T_platelet_gene"))] <- "CD4+_T_unknown"

Idents(object = seurat_AIDA_subset_corecelltypes) <- "Annotation_Level3"

int_most_enriched <- which(da_results$logFC == max(da_results$logFC))
vec_cells_interest <- names(which(nhoods(milo_AIDA_subset)[, int_most_enriched] > 0))
vec_cells_interest <- vec_cells_interest[vec_cells_interest %in% colnames(seurat_AIDA_subset_corecelltypes)]

DimPlot(seurat_AIDA_subset_corecelltypes, reduction = "umap", 
        cells.highlight = vec_cells_interest, sizes.highlight = 0.8, 
        label = TRUE, raster = FALSE, repel = TRUE, label.size = 8) + NoLegend()

ggsave(paste("Rplot_MiloR_UMAP_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".png", sep = "_"))

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "RNA"
da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level3")

### Identifying cell subtype implicated in most enriched neighbourhood
vec_celltype <- names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level3 == da_results$Annotation_Level3[int_most_enriched]))
vec_other <- vec_celltype [! vec_celltype %in% vec_cells_interest]
df_markergenes <- FindMarkers(seurat_AIDA_subset_corecelltypes, ident.1 = vec_cells_interest, ident.2 = vec_other)

vec_high_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = avg_log2FC)
vec_low_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = -avg_log2FC)

### Identifying annotated cell type of most enriched neighbourhood
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_cells_interest) <- paste0("Most enriched ", da_results$Annotation_Level3[int_most_enriched])
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_other) <- paste0("Other ", da_results$Annotation_Level3[int_most_enriched])
DotPlot(seurat_AIDA_subset_corecelltypes, features = c(rownames(vec_high_genes), rownames(vec_low_genes))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste("Rplot_MiloR_DotPlot_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"))

### Myeloid_pDC, 50 to 77

str_covariate_comparison <- "Ages_50_to_77"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Ancestry + Sex + Ages_50_to_77, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

ggplot()+
  geom_point(df_UMAP, 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0), 
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of pDC and myeloid cells in donors of ", str_covariate_comparison))

da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level4")

table(da_results$Annotation_Level4)

vec_cells_to_remove <- c(which(da_results$Annotation_Level4 == "flagged_CD14+_Monocyte_cDC2_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD14+_Monocyte_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD14+_Monocyte_TNK_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD16+_Monocyte_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD16+_Monocyte_TNK_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_myeloid_CD14+_Monocyte_cDC2_CD16+_Monocyte_gene"))

plotDAbeeswarm(da_results[-vec_cells_to_remove, ], group.by = "Annotation_Level4") +
  ggtitle(paste0("Cell neighbourhood enrichment of pDC and myeloid cells in donors ", str_covariate_comparison))

ggsave(paste("Rplot_Figure_MiloR_Beeswarmby_Annotation_Level4", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"), 
       width = 13.34, height = 11.24, units = "in")

### SG_Malay, Myeloid_pDC

str_covariate_comparison <- "SG_Malay"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Sex + SG_Malay, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

ggplot()+
  geom_point(df_UMAP, 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0), 
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of pDC and myeloid cells in ", str_covariate_comparison, " donors"))

### Plotting DotPlot of DEG marking out most enriched neighbourhood

vec_cells_to_keep <- c(names(which(seurat_AIDA_subset$Annotation_Level2 == "cDC")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "DC")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "Monocyte")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "pDC")))

seurat_AIDA_subset_corecelltypes <- seurat_AIDA_subset#[, vec_cells_to_keep]

### UMAP of most enriched cell neighbourhood

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "integrated"
seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level2 == "Myeloid"))] <- "Myeloid_unknown"

Idents(object = seurat_AIDA_subset_corecelltypes) <- "Annotation_Level3"

int_most_enriched <- which(da_results$logFC == max(da_results$logFC))
vec_cells_interest <- names(which(nhoods(milo_AIDA_subset)[, int_most_enriched] > 0))
vec_cells_interest <- vec_cells_interest[vec_cells_interest %in% colnames(seurat_AIDA_subset_corecelltypes)]

DimPlot(seurat_AIDA_subset_corecelltypes, reduction = "umap", 
        cells.highlight = vec_cells_interest, sizes.highlight = 0.8, 
        label = TRUE, raster = FALSE, repel = TRUE, label.size = 8) + NoLegend()

ggsave(paste("Rplot_MiloR_UMAP_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".png", sep = "_"))

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "RNA"
da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level3")

### Identifying cell subtype implicated in most enriched neighbourhood
vec_celltype <- names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level3 == da_results$Annotation_Level3[int_most_enriched]))
vec_other <- vec_celltype [! vec_celltype %in% vec_cells_interest]
df_markergenes <- FindMarkers(seurat_AIDA_subset_corecelltypes, ident.1 = vec_cells_interest, ident.2 = vec_other)

vec_high_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = avg_log2FC)
vec_low_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = -avg_log2FC)

vec_cells_to_keep <- c(names(which(seurat_AIDA_subset$Annotation_Level2 == "cDC")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "DC")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "Monocyte")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "pDC")))

### Identifying annotated cell type of most enriched neighbourhood
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_cells_interest) <- paste0("Most enriched ", da_results$Annotation_Level3[int_most_enriched])
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_other) <- paste0("Other ", da_results$Annotation_Level3[int_most_enriched])
DotPlot(seurat_AIDA_subset_corecelltypes[, vec_cells_to_keep], features = c(rownames(vec_high_genes), rownames(vec_low_genes))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste("Rplot_MiloR_DotPlot_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"))

### SG_Indian, B

str_covariate_comparison <- "SG_Indian"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Sex + SG_Indian, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

ggplot()+
  geom_point(df_UMAP, 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0), 
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of B cells in ", str_covariate_comparison, " donors"))

### Plotting DotPlot of DEG marking out most enriched neighbourhood

seurat_AIDA_subset_corecelltypes <- seurat_AIDA_subset#[, vec_cells_to_keep]

### UMAP of most enriched cell neighbourhood

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "integrated"
Idents(object = seurat_AIDA_subset_corecelltypes) <- "Annotation_Level3"

int_most_enriched <- which(da_results$logFC == max(da_results$logFC))
vec_cells_interest <- names(which(nhoods(milo_AIDA_subset)[, int_most_enriched] > 0))
vec_cells_interest <- vec_cells_interest[vec_cells_interest %in% colnames(seurat_AIDA_subset_corecelltypes)]

seurat_AIDA_subset_corecelltypes$Annotation_Level3[names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level3 == "B"))] <- "B_unknown"
Idents(seurat_AIDA_subset_corecelltypes) <- seurat_AIDA_subset_corecelltypes$Annotation_Level3

DimPlot(seurat_AIDA_subset_corecelltypes, reduction = "umap", 
        cells.highlight = vec_cells_interest, sizes.highlight = 0.8, 
        label = TRUE, raster = FALSE, repel = TRUE, label.size = 8) + NoLegend()

ggsave(paste("Rplot_MiloR_UMAP_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".png", sep = "_"))

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "RNA"
da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level3")

### Identifying cell subtype implicated in most enriched neighbourhood
vec_celltype <- names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level3 == da_results$Annotation_Level3[int_most_enriched]))
vec_other <- vec_celltype [! vec_celltype %in% vec_cells_interest]
df_markergenes <- FindMarkers(seurat_AIDA_subset_corecelltypes, ident.1 = vec_cells_interest, ident.2 = vec_other)

vec_high_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = avg_log2FC)
vec_low_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = -avg_log2FC)

vec_cells_to_keep <- c(names(which(seurat_AIDA_subset$Annotation_Level2 == "atypical_B")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "memory_B")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "naive_B")))

### Identifying annotated cell type of most enriched neighbourhood
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_cells_interest) <- paste0("Most enriched ", da_results$Annotation_Level3[int_most_enriched])
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_other) <- paste0("Other ", da_results$Annotation_Level3[int_most_enriched])
DotPlot(seurat_AIDA_subset_corecelltypes[, vec_cells_to_keep], features = c(rownames(vec_high_genes), rownames(vec_low_genes))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste("Rplot_MiloR_DotPlot_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"))

### Male, B

str_covariate_comparison <- "Male"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Ancestry + Male, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

ggplot()+
  geom_point(df_UMAP, 
             mapping = aes(x = UMAP_1, y = UMAP_2, 
                           colour = log2_MeanFoldChange), 
             size = 0.5, stroke = 0.1) +
  scale_colour_gradientn(colours = c("darkorange","#e8ded3", "#E8E8E8", "#d5e0e8", "darkblue"),
                         values = c(1, 0.54, 0.5, 0.46, 0), 
                         limits = c(-2, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggtitle(paste0("Enrichment of B cells in ", str_covariate_comparison, " donors"))

### Plotting DotPlot of DEG marking out most enriched neighbourhood

seurat_AIDA_subset_corecelltypes <- seurat_AIDA_subset#[, vec_cells_to_keep]

### UMAP of most enriched cell neighbourhood

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "integrated"
Idents(object = seurat_AIDA_subset_corecelltypes) <- "Annotation_Level3"

int_most_enriched <- which(da_results$logFC == max(da_results$logFC))
vec_cells_interest <- names(which(nhoods(milo_AIDA_subset)[, int_most_enriched] > 0))
vec_cells_interest <- vec_cells_interest[vec_cells_interest %in% colnames(seurat_AIDA_subset_corecelltypes)]

DimPlot(seurat_AIDA_subset_corecelltypes, reduction = "umap", 
        cells.highlight = vec_cells_interest, sizes.highlight = 0.8, 
        label = TRUE, raster = FALSE, repel = TRUE, label.size = 8) + NoLegend()

ggsave(paste("Rplot_MiloR_UMAP_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".png", sep = "_"))

DefaultAssay(seurat_AIDA_subset_corecelltypes) <- "RNA"
da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level3")

vec_cells_to_remove <- c(which(da_results$Annotation_Level3 == "B"))

plotDAbeeswarm(da_results[-vec_cells_to_remove, ], group.by = "Annotation_Level3") +
  ggtitle(paste0("Cell neighbourhood enrichment of B cells in ", str_covariate_comparison, " donors"))
ggsave(paste("Rplot_Figure_MiloR_Beeswarmby_Annotation_Level3", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"), 
       width = 13.34, height = 11.24, units = "in")

### Identifying cell subtype implicated in most enriched neighbourhood
vec_celltype <- names(which(seurat_AIDA_subset_corecelltypes$Annotation_Level3 == da_results$Annotation_Level3[int_most_enriched]))
vec_other <- vec_celltype [! vec_celltype %in% vec_cells_interest]
df_markergenes <- FindMarkers(seurat_AIDA_subset_corecelltypes, ident.1 = vec_cells_interest, ident.2 = vec_other)

vec_high_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = avg_log2FC)
vec_low_genes <- df_markergenes[which(df_markergenes$p_val_adj < 0.05), ] %>%
  slice_max(n = 5, order_by = -avg_log2FC)

vec_cells_to_keep <- c(names(which(seurat_AIDA_subset$Annotation_Level2 == "atypical_B")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "memory_B")),
                       names(which(seurat_AIDA_subset$Annotation_Level2 == "naive_B")))

### Identifying annotated cell type of most enriched neighbourhood
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_cells_interest) <- paste0("Most enriched ", da_results$Annotation_Level3[int_most_enriched])
Idents(object = seurat_AIDA_subset_corecelltypes, cells = vec_other) <- paste0("Other ", da_results$Annotation_Level3[int_most_enriched])
DotPlot(seurat_AIDA_subset_corecelltypes[, vec_cells_to_keep], features = c(rownames(vec_high_genes), rownames(vec_low_genes))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste("Rplot_MiloR_DotPlot_most_enriched", str_covariate_comparison, 
             str_celltype, "k", str_k, ".pdf", sep = "_"))

### Male, CD8+NK

str_covariate_comparison <- "Male"

da_results <- testNhoods(milo_AIDA_subset, 
                         design = ~ Age + Ancestry + Male, 
                         design.df = traj_design, fdr.weighting = "graph-overlap")

da_results <- annotateNhoods(milo_AIDA_subset, da_results, coldata_col = "Annotation_Level4")

vec_cells_to_remove <- c(which(da_results$Annotation_Level4 == "CD8+_T"), 
                         which(da_results$Annotation_Level4 == "flagged_CD16+_NK_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_gdT_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_CD8+_T_naive_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_NK_low_exp"), 
                         which(da_results$Annotation_Level4 == "flagged_T_CD8+_T_gdT_gene"), 
                         which(da_results$Annotation_Level4 == "flagged_T_platelet_gene"), 
                         which(da_results$Annotation_Level4 == "T_IFNhi"))

plotDAbeeswarm(da_results[-vec_cells_to_remove, ], group.by = "Annotation_Level4") +
  ggtitle(paste0("Cell neighbourhood enrichment of CD8+ T, gdT, ILC, and NK cells in ", str_covariate_comparison, " donors"))

ggsave(paste("Rplot_Figure_MiloR_Beeswarmby_Annotation_Level4", 
             str_covariate_comparison, str_celltype, "k", str_k, ".png", sep = "_"), 
       width = 26.68, height = 11.24, units = "in")
