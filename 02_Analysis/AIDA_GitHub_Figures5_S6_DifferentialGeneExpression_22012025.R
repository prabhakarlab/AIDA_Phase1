### AIDA Phase 1 Data Freeze v2 Differential gene expression analysis: version 22 January 2025

##### INPUT REQUIRED: READ IN FULL ANNOTATION TABLE AND METADATA AS df_AIDA_metadata #####

### Consider only donors with at least 800 cells and population groups with at least 50 donors

df_AIDA_metadata$count <- 1
rownames(df_AIDA_metadata) <- df_AIDA_metadata$barcode_name

df_total_cell_count <- df_AIDA_metadata[, c("DCP_ID", "count")]
df_total_cell_count <- aggregate(count ~ ., df_total_cell_count, FUN = sum)
rownames(df_total_cell_count) <- df_total_cell_count$DCP_ID
vec_donors_to_exclude <- df_total_cell_count$DCP_ID[which(df_total_cell_count$count < 800)]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$DCP_ID %in% vec_donors_to_exclude), ]
df_AIDA_metadata <- df_AIDA_metadata[(df_AIDA_metadata$Country == "SG"), ]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$ethnicity == "European"), ]

library(Seurat)
library(edgeR)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
require(DOSE)
vec_colours <- c("#F8766D", "#B79F00", "#619CFF")

cyt.go.genes <- as.list(org.Hs.egGO2ALLEGS)

##### INPUT REQUIRED: READ IN THE RELEVANT SEURAT OBJECT FOR EACH CELL POPULATION AS seurat_majorcelltype #####

vec_barcodes_to_keep <- colnames(seurat_majorcelltype)[colnames(seurat_majorcelltype) %in% rownames(df_AIDA_metadata)]
seurat_majorcelltype <- seurat_majorcelltype[, vec_barcodes_to_keep]
seurat_majorcelltype$Annotation_Level3 <- df_AIDA_metadata[colnames(seurat_majorcelltype), "Annotation_Level3"]

# for (str_subtype in c("atypical_B", "memory_B_IGHMhi", "memory_B_IGHMlo", "naive_B"))
# for (str_subtype in c("CD14+_Monocyte", "CD16+_Monocyte", "cDC2", "pDC"))
for (str_subtype in c("CD16+_NK", "CD4+_T_cm", "CD4+_T_cyt", "CD4+_T_em", "CD4+_T_naive", 
                      "CD56+_NK", "CD8+_T_GZMBhi", "CD8+_T_GZMKhi", "CD8+_T_naive", 
                      "gdT_GZMBhi", "gdT_GZMKhi", "MAIT", "Treg")){
  
  seurat_subset <- seurat_majorcelltype[, which(seurat_majorcelltype$Annotation_Level3 == str_subtype)]
  df_design <- seurat_subset[[]]
  df_design$Batch <- sapply(df_design$Library, FUN = function(x){strsplit(x, split = "_L00")[[1]][1]})
  df_design <- df_design[, c("DCP_ID", "Age", "Sex", "ethnicity", "Batch")]
  df_design <- unique(df_design)
  rownames(df_design) <- df_design$DCP_ID
  df_cellcount <- data.frame(table(seurat_subset[[]]$DCP_ID))
  df_cellcount <- df_cellcount[which(df_cellcount$Freq >= 10), ]
  df_design <- df_design[df_design$DCP_ID %in% df_cellcount$Var1, ]
  df_design$SG_Chinese <- "Not_SG_Chinese"
  df_design$SG_Chinese[which(df_design$ethnicity == "SG_Chinese")] <- "SG_Chinese"
  df_design$Not_SG_Chinese <- 1
  df_design$Not_SG_Chinese[which(df_design$ethnicity == "SG_Chinese")] <- 0
  
  df_design$SG_Malay <- "Not_SG_Malay"
  df_design$SG_Malay[which(df_design$ethnicity == "SG_Malay")] <- "SG_Malay"
  df_design$Not_SG_Malay <- 1
  df_design$Not_SG_Malay[which(df_design$ethnicity == "SG_Malay")] <- 0
  
  df_design$SG_Indian <- "Not_SG_Indian"
  df_design$SG_Indian[which(df_design$ethnicity == "SG_Indian")] <- "SG_Indian"
  df_design$Not_SG_Indian <- 1
  df_design$Not_SG_Indian[which(df_design$ethnicity == "SG_Indian")] <- 0
  
  seurat_subset
  
  ##### Pseudobulk of count matrix from gene-cell to gene-donor #####
  
  matrix_counts <- AggregateExpression(seurat_subset, 
                                       assays = "RNA", 
                                       group.by = "DCP_ID", 
                                       slot = "counts")$RNA
  matrix_counts <- t(matrix_counts)
  matrix_counts <- matrix_counts[rownames(df_design), ]
  
  vec_donors_per_gene <- colSums(matrix_counts > 0)
  vec_UMIs_per_gene <- colSums(matrix_counts)
  
  vec_genes_to_keep <- intersect(names(which(vec_donors_per_gene >= 0.1*dim(df_design)[1])), 
                                 names(which(vec_UMIs_per_gene >= dim(df_design)[1])))
  matrix_counts <- matrix_counts[, vec_genes_to_keep]
  
  ##### edgeR analyses #####
  
  y <- DGEList(counts = t(matrix_counts))
  y <- calcNormFactors(y)
  
  for (str_ethnicity in c("SG_Chinese", "SG_Malay", "SG_Indian")){
    df_design[, str_ethnicity] <- factor(df_design[, str_ethnicity], 
                                         levels=c(paste0("Not_", str_ethnicity), 
                                                  str_ethnicity))
    design <- model.matrix(formula(paste0("~Sex + Age + Batch + ", str_ethnicity)), 
                           data=df_design)
    
    y <- estimateDisp(y, design)
    
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit)
    topTags(lrt)
    write.table(lrt$table, 
                file = paste("edgeR", str_subtype, str_ethnicity, ".txt", sep = "_"), 
                col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    
    group <- factor(y$design[, "SexMale"])
    colours <- c("red", "blue")
    points <- c(15, 16)
    ### This shows the segregation of profiles primarily by sex
    png(paste("Rplot_MDS_Sex", str_ethnicity, str_subtype, ".png", sep = "_"))
    plotMDS(y, col=colours[group], pch=points[group])
    title("MDS plot of sex")
    dev.off()
    
    group <- factor(y$design[, paste0(str_ethnicity, str_ethnicity)])
    colours <- c("black", "pink")
    points <- c(15, 16)
    png(paste("Rplot_MDS_ethnicity", str_ethnicity, str_subtype, ".png", sep = "_"))
    plotMDS(y, col=colours[group], pch=points[group])
    title("MDS plot of ethnicity")
    dev.off()
    
    # MA plot
    png(paste("Rplot_MA_ethnicity", str_ethnicity, str_subtype, ".png", sep = "_"))
    plot(lrt$table$logCPM, lrt$table$logFC)
    dev.off()
      }
    }

##### Preparing object for GO term analysis #####
{
  df_convert <- select(org.Hs.eg.db, keys = rownames((y$counts)), 
                       columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
  which(duplicated(df_convert$SYMBOL))
  vec_duplicates <- c()
  for (str_index in which(duplicated(df_convert$SYMBOL))){
    vec_duplicates <- c(vec_duplicates, 
                        which(df_convert$SYMBOL == df_convert[str_index, "SYMBOL"])[1])
  }
  df_convert <- df_convert[-vec_duplicates, ]
  y_GO <- y
  y_GO$counts <- y_GO$counts[rownames(y_GO)[!is.na(df_convert$ENTREZID)], ]
  rownames(y_GO$counts) <- df_convert$ENTREZID[!is.na(df_convert$ENTREZID)]
  
  y_GO <-  calcNormFactors(y_GO)
  y_GO <- estimateDisp(y_GO, design)
  
  fit_GO <- glmFit(y_GO, design)
  lrt_GO <- glmLRT(fit_GO)
  topTags(lrt_GO)
  
  go <- goana(lrt_GO, species = "Hs")
  df_topGO <- topGO(go, ont = "BP", sort = "Up", n = 30, truncate = 30)
  df_topGO
}

### Examining DEGs

vec_filenames <- grep("^edgeR_", list.files(), value = TRUE)

for (str_filename in vec_filenames){
  print(str_filename)
  df_edgeR <- read.table(str_filename, header = TRUE, sep = "\t")
  df_edgeR$gene <- row.names(df_edgeR)
  df_edgeR$FDR <- p.adjust(df_edgeR$PValue, method = "fdr")
  df_edgeR <- df_edgeR[which(df_edgeR$FDR < 0.05), ]
  if (dim(df_edgeR)[1] > 0){
    
    str_to_split <- strsplit(str_filename, split = "_.txt")[[1]][1]
    str_to_split <- strsplit(str_to_split, split = "edgeR_")[[1]][2]
    df_edgeR$ethnicity <- strsplit(str_to_split, split = "_")[[1]][length(strsplit(str_to_split, split = "_")[[1]])]
    df_edgeR$subtype <- paste(strsplit(str_to_split, split = "_")[[1]][1:length(strsplit(str_to_split, split = "_")[[1]])-1], collapse = "_")
    
    if (match(str_filename, vec_filenames) == 1){
      df_strong_DEG <- df_edgeR } else {
        df_strong_DEG <- rbind(df_strong_DEG, df_edgeR)
      }
  }
}

write.table(df_strong_DEG, 
            "AIDA_AtlasManuscript_SupplementaryTable_DEG_FDR005.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

for (str_ethnicity in c("SG_Chinese", "SG_Indian", "SG_Malay")){
  #  for (str_celltype in c("CD16+_NK", "CD14+_Monocyte", "CD4+_T_cm", "naive_B")){
  for (str_celltype in c("atypical_B", "memory_B_IGHMhi", "memory_B_IGHMlo", "naive_B", 
                         "CD14+_Monocyte", "CD16+_Monocyte", "cDC2", "pDC",
                         "CD16+_NK", "CD4+_T_cm", "CD4+_T_cyt", "CD4+_T_em", "CD4+_T_naive", 
                         "CD56+_NK", "CD8+_T_GZMBhi", "CD8+_T_GZMKhi", "CD8+_T_naive", 
                         "gdT_GZMBhi", "gdT_GZMKhi", "MAIT", "Treg")){
    print(str_ethnicity)
    print(str_celltype)
    
    df_test <- read.table(paste("edgeR", 
                                str_celltype, str_ethnicity, ".txt", sep = "_"), header = TRUE, sep = "\t")
    
    df_test$logFC <- as.numeric(df_test$logFC)
    df_test$PValue <- as.numeric(df_test$PValue)
    df_test$FDR <- p.adjust(df_test$PValue, method = "fdr")
    df_test$FDR_threshold <- ifelse(df_test$FDR < 0.05, "Yes", "No")
    df_test$Sign <- ifelse(df_test$logFC > 0, 1, -1)
    
    ggplot(data = df_test, 
           aes(x = logFC, y = -log10(PValue), colour = FDR_threshold), 
           label = rownames(df_test)) + geom_point(size = 0.8, shape = 20) + 
      theme_minimal() + theme(legend.position = "none") + 
      scale_colour_manual(values = c("Black", "#F8766D")) + 
      geom_label_repel(aes(label = ifelse(((abs(logFC) >= 0.5) & (FDR_threshold == "Yes")), 
                                          rownames(df_test), "")), max.overlaps = Inf) +
      ggtitle(paste0("Differentially expressed genes in ", str_celltype,  
                     ": \n", str_ethnicity, " versus other Singapore ethnicities"))
    ggsave(paste("Volcano_GSEA_output/Rplot_Fig5_VolcanoPlot", 
                 str_ethnicity, str_celltype, ".pdf", sep = "_"), 
           height = 6.36, width = 11.96, units = "in")

    ##### GO term analysis #####
    
    # P-value version (ranking p-values by most significant by using -log10 transformation, then multiply by fold-change sign)
    vec_gene_list <- -log10(df_test$PValue)*df_test$Sign
    names(vec_gene_list) <- rownames(df_test)
    vec_gene_list <- sort(vec_gene_list, decreasing = TRUE)
    
    set.seed(1234)
    gse <- gseGO(geneList = vec_gene_list, 
                 ont = "BP", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "fdr", 
                 seed = TRUE)
    
    if (dim(gse@result)[1] > 0){
      dotplot(gse, showCategory = 10, split = ".sign") + facet_grid( . ~ .sign) + 
        ggtitle(paste0("GSEA (GO biological process) of gene expression in ", str_celltype,  
                       ": \n", str_ethnicity, " versus other Singapore ethnicities"))
      ggsave(paste("Volcano_GSEA_output/Rplot_Fig5_GSEA_log10pvalue", str_ethnicity, str_celltype, ".pdf", sep = "_"), 
             height = 9.54, width = 11.96, units = "in")
    }
  }
}

### Adjusting plot height for SG_Indian MAIT, SG_Chinese CD16+_NK, SG_Malay CD14+_Monocyte GSEA

str_ethnicity <- "SG_Indian"
str_celltype <- "MAIT"

str_ethnicity <- "SG_Chinese"
str_celltype <- "CD16+_NK"

str_ethnicity <- "SG_Malay"
str_celltype <- "CD14+_Monocyte"

dotplot(gse, showCategory = 5, split = ".sign") + facet_grid( . ~ .sign) + 
  ggtitle(paste0("GSEA (GO biological process) of gene expression in ", str_celltype,  
                 ": \n", str_ethnicity, " versus other Singapore ethnicities"))
ggsave(paste("Rplot_Fig5_GSEA_log10pvalue", str_ethnicity, str_celltype, ".pdf", sep = "_"), 
       height = 6.36, width = 11.96, units = "in")

df_test <- gse@result
df_test$Number_Core_Enrichment <- sapply(df_test$core_enrichment, 
                                         FUN = function(x){length(unlist(strsplit(x, split = "/")))})
df_test$GeneRatio <- df_test$Number_Core_Enrichment/df_test$setSize
df_test <- df_test[order(df_test$pvalue),]
df_test[1:20, c("Description", "enrichmentScore", "pvalue")]
df_test <- df_test[order(-df_test$qvalue),]
df_test[1:10, c("Description", "enrichmentScore", "pvalue")]

df_test <- df_test[order(-df_test$GeneRatio),]
df_test[1:20, c("Description", "enrichmentScore", "pvalue", "GeneRatio")]

### Adjusting highlighted points for CD14+_Monocyte SG_Malay and CD4+_T_naive SG_Indian:

str_ethnicity <- "SG_Malay"
str_celltype <- "CD14+_Monocyte"

str_ethnicity <- "SG_Chinese"
str_celltype <- "CD16+_NK"

str_ethnicity <- "SG_Indian"
str_celltype <- "MAIT"

print(str_ethnicity)
print(str_celltype)

df_test <- read.table(paste("edgeR", 
                            str_celltype, str_ethnicity, ".txt", sep = "_"), header = TRUE, sep = "\t")

df_test$logFC <- as.numeric(df_test$logFC)
df_test$PValue <- as.numeric(df_test$PValue)
df_test$FDR <- p.adjust(df_test$PValue, method = "fdr")
df_test$FDR_threshold <- ifelse(df_test$FDR < 0.05, "Yes", "No")
df_test$Sign <- ifelse(df_test$logFC > 0, 1, -1)

ggplot(data = df_test, 
       aes(x = logFC, y = -log10(PValue), colour = FDR_threshold), 
       label = rownames(df_test)) + geom_point(size = 0.8, shape = 20) + 
  theme_minimal() + theme(legend.position = "none") + 
  scale_colour_manual(values = c("Black", "#F8766D")) + 
  geom_label_repel(aes(label = ifelse(((abs(logFC) >= 0.0) & (FDR_threshold == "Yes")), 
                                      rownames(df_test), "")), max.overlaps = Inf) +
  xlab("Log2 fold-change") + ylab("-log10(p-value)") + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18)) + 
  ggtitle(paste0("Differentially expressed genes in ", str_celltype,  
                 ": ", str_ethnicity, " versus other Singapore ethnicities"))

ggplot(data = df_test, 
       aes(x = logFC, y = -log10(PValue), colour = FDR_threshold), 
       label = rownames(df_test)) + geom_point(size = 0.8, shape = 20) + 
  theme_minimal() + theme(legend.position = "none") + 
  scale_colour_manual(values = c("Black", "#F8766D")) + 
  geom_label_repel(aes(label = ifelse(((abs(logFC) >= 0.75) & (FDR_threshold == "Yes")), 
                                      rownames(df_test), "")), max.overlaps = Inf) +
  xlab("Log2 fold-change") + ylab("-log10(p-value)") + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18)) + 
  ggtitle(paste0("Differentially expressed genes in ", str_celltype,  
                 ": ", str_ethnicity, " versus other Singapore ethnicities"))

ggsave(paste("Rplot_Fig5_VolcanoPlot", 
             str_ethnicity, str_celltype, ".pdf", sep = "_"), 
       height = 6.36, width = 11.96, units = "in")

### Getting UTS2 values from each ethnicity

vec_UTS2_ethnicity <- c()
vec_UTS2_celltype <- c()

for (str_ethnicity in c("SG_Chinese", "SG_Indian", "SG_Malay")){
  for (str_celltype in c("atypical_B", "memory_B_IGHMhi", "memory_B_IGHMlo", "naive_B", 
                         "CD14+_Monocyte", "CD16+_Monocyte", "cDC2", "pDC",
                         "CD16+_NK", "CD4+_T_cm", "CD4+_T_cyt", "CD4+_T_em", "CD4+_T_naive", 
                         "CD56+_NK", "CD8+_T_GZMBhi", "CD8+_T_GZMKhi", "CD8+_T_naive", 
                         "gdT_GZMBhi", "gdT_GZMKhi", "MAIT", "Treg")){
    print(str_ethnicity)
    print(str_celltype)
    
    df_test <- read.table(paste("edgeR", 
                                str_celltype, str_ethnicity, ".txt", sep = "_"), header = TRUE, sep = "\t")
    print(df_test[which(rownames(df_test) == "UTS2"), ])
    
    if (dim(df_test[which(rownames(df_test) == "UTS2"), ])[1] > 0){
      if ((str_ethnicity == "SG_Chinese") & (str_celltype == "memory_B_IGHMlo")){
        print("Yes")
        df_UTS2 <- df_test[which(rownames(df_test) == "UTS2"), ]
        vec_UTS2_ethnicity <- c(vec_UTS2_ethnicity, str_ethnicity)
        vec_UTS2_celltype <- c(vec_UTS2_celltype, str_celltype)
      } else {
        df_UTS2 <- rbind(df_UTS2, df_test[which(rownames(df_test) == "UTS2"), ])
        vec_UTS2_ethnicity <- c(vec_UTS2_ethnicity, str_ethnicity)
        vec_UTS2_celltype <- c(vec_UTS2_celltype, str_celltype)
      }
    }
  }
}

df_UTS2$ethnicity <- vec_UTS2_ethnicity
df_UTS2$Subtype <- vec_UTS2_celltype

df_UTS2$Subtype <- factor(df_UTS2$Subtype, 
                          levels = c("naive_B", "memory_B_IGHMhi",
                                     "memory_B_IGHMlo", "atypical_B", 
                                     "CD14+_Monocyte", "CD16+_Monocyte", 
                                     "cDC2", "pDC", "CD4+_T_naive", 
                                     "CD4+_T_cm", "CD4+_T_em", "CD4+_T_cyt", 
                                     "Treg", "CD8+_T_naive", "CD8+_T_GZMKhi", 
                                     "CD8+_T_GZMBhi", "MAIT", "gdT_GZMKhi", 
                                     "gdT_GZMBhi", "CD16+_NK", "CD56+_NK"))

ggplot(df_UTS2, aes(y = logFC, x = logCPM, colour = ethnicity, shape = Subtype)) + 
  geom_point(size = 4.9) + 
  scale_shape_manual(values = c(15:17, 1:4, 5:14)) + 
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF")) + 
  scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF")) + 
  ggtitle("Log2 fold-change against log2 (counts per million) of UTS2 \nby ethnicity and cell subtype") + 
  ylab("Log2 fold-change") + xlab("Log2 (Counts per million)") + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18), 
        legend.text = element_text(size = 10)) + 
  guides(col = guide_legend(ncol = 3))

ggsave(paste("Rplot_Fig4_UTS2.pdf"), 
       height = 6.36, width = 11.96, units = "in")

### Getting TIFA values from each ethnicity

vec_TIFA_ethnicity <- c()
vec_TIFA_celltype <- c()

for (str_ethnicity in c("SG_Chinese", "SG_Indian", "SG_Malay")){
  for (str_celltype in c("atypical_B", "memory_B_IGHMhi", "memory_B_IGHMlo", "naive_B", 
                         "CD14+_Monocyte", "CD16+_Monocyte", "cDC2", "pDC",
                         "CD16+_NK", "CD4+_T_cm", "CD4+_T_cyt", "CD4+_T_em", "CD4+_T_naive", 
                         "CD56+_NK", "CD8+_T_GZMBhi", "CD8+_T_GZMKhi", "CD8+_T_naive", 
                         "gdT_GZMBhi", "gdT_GZMKhi", "MAIT", "Treg")){
    print(str_ethnicity)
    print(str_celltype)
    
    df_test <- read.table(paste("edgeR", 
                                str_celltype, str_ethnicity, ".txt", sep = "_"), header = TRUE, sep = "\t")
    print(df_test[which(rownames(df_test) == "TIFA"), ])
    
    if (dim(df_test[which(rownames(df_test) == "TIFA"), ])[1] > 0){
      if ((str_ethnicity == "SG_Chinese") & (str_celltype == "atypical_B")){
        print("Yes")
        df_TIFA <- df_test[which(rownames(df_test) == "TIFA"), ]
        vec_TIFA_ethnicity <- c(vec_TIFA_ethnicity, str_ethnicity)
        vec_TIFA_celltype <- c(vec_TIFA_celltype, str_celltype)
      } else {
        df_TIFA <- rbind(df_TIFA, df_test[which(rownames(df_test) == "TIFA"), ])
        vec_TIFA_ethnicity <- c(vec_TIFA_ethnicity, str_ethnicity)
        vec_TIFA_celltype <- c(vec_TIFA_celltype, str_celltype)
      }
    }
  }
}

df_TIFA$ethnicity <- vec_TIFA_ethnicity
df_TIFA$Subtype <- vec_TIFA_celltype

df_TIFA$Subtype <- factor(df_TIFA$Subtype, 
                          levels = c("naive_B", "memory_B_IGHMhi",
                                     "memory_B_IGHMlo", "atypical_B", 
                                     "CD14+_Monocyte", "CD16+_Monocyte", 
                                     "cDC2", "pDC", "CD4+_T_naive", 
                                     "CD4+_T_cm", "CD4+_T_em", "CD4+_T_cyt", 
                                     "Treg", "CD8+_T_naive", "CD8+_T_GZMKhi", 
                                     "CD8+_T_GZMBhi", "MAIT", "gdT_GZMKhi", 
                                     "gdT_GZMBhi", "CD16+_NK", "CD56+_NK"))

ggplot(df_TIFA, aes(y = logFC, x = logCPM, colour = ethnicity, shape = Subtype)) + 
  geom_point(size = 4.9) + 
  scale_shape_manual(values = c(15, 18, 16, 20, 17, 0, 1, 97, 2:14)) + 
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF")) + 
  scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF")) + 
  ggtitle("Log2 fold-change against log2 (counts per million) of \nTIFA by ethnicity and cell subtype") + 
  ylab("Log2 fold-change") + xlab("Log2 (Counts per million)") + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18), 
        legend.text = element_text(size = 10)) + 
  guides(col = guide_legend(ncol = 3))

ggsave(paste("Rplot_Fig4_TIFA.pdf"), 
       height = 6.36, width = 11.96, units = "in")

ggplot(df_TIFA[which(df_TIFA$Subtype %in% c("naive_B", "memory_B_IGHMhi",
                                            "memory_B_IGHMlo", "atypical_B")), ], 
       aes(y = logFC, x = logCPM, colour = ethnicity, shape = Subtype)) + 
  geom_point(size = 3) + 
  scale_shape_manual(values = c(15, 18, 16, 20)) + 
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF")) + 
  scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF")) + 
  ggtitle("Log2 fold-change against log2 (counts per million) of TIFA \nby ethnicity and cell subtype") + 
  ylab("Log2 fold-change") + xlab("Log2 (Counts per million)")

ggsave("Rplot_Fig4_TIFA_Bcells.pdf", 
       height = 3.18, width = 5.98, units = "in")

##### edgeR analyses: pairwise between Singapore population groups #####

### Consider only donors with at least 800 cells and ethnicities with at least 40 donors per ethnicity

df_total_cell_count <- df_AIDA_metadata[, c("DCP_ID", "count")]
df_total_cell_count <- aggregate(count ~ ., df_total_cell_count, FUN = sum)
rownames(df_total_cell_count) <- df_total_cell_count$DCP_ID
vec_donors_to_exclude <- df_total_cell_count$DCP_ID[which(df_total_cell_count$count < 800)]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$DCP_ID %in% vec_donors_to_exclude), ]
df_AIDA_metadata <- df_AIDA_metadata[(df_AIDA_metadata$Country == "SG"), ]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$Ethnicity == "European"), ]

##### INPUT REQUIRED: READ IN THE RELEVANT SEURAT OBJECT FOR EACH CELL POPULATION AS seurat_majorcelltype #####

vec_barcodes_to_keep <- colnames(seurat_majorcelltype)[colnames(seurat_majorcelltype) %in% rownames(df_AIDA_metadata)]
seurat_majorcelltype <- seurat_majorcelltype[, vec_barcodes_to_keep]
seurat_majorcelltype$Annotation_Level3 <- df_AIDA_metadata[colnames(seurat_majorcelltype), "Annotation_Level3"]

list_pairwise <- list()
list_pairwise[[1]] <- c("Chinese", "Malay")
list_pairwise[[2]] <- c("Chinese", "Indian")
list_pairwise[[3]] <- c("Malay", "Indian")

# for (str_subtype in c("atypical_B", "memory_B_IGHMhi", "memory_B_IGHMlo", "naive_B"))
# for (str_subtype in c("CD14+_Monocyte", "CD16+_Monocyte", "cDC2", "pDC"))
for (str_subtype in c("CD16+_NK", "CD4+_T_cm", "CD4+_T_cyt", "CD4+_T_em", "CD4+_T_naive", 
                      "CD56+_NK", "CD8+_T_GZMBhi", "CD8+_T_GZMKhi", "CD8+_T_naive", 
                      "gdT_GZMBhi", "gdT_GZMKhi", "MAIT", "Treg")){
  
  for (vec_pair in list_pairwise){
    
    seurat_subset <- seurat_majorcelltype[, which(seurat_majorcelltype$Annotation_Level3 == str_subtype)]
    
    print(str_subtype)
    print(vec_pair)
    seurat_subset <- seurat_subset[, c(which(seurat_subset$Ethnicity == vec_pair[1]), 
                                       which(seurat_subset$Ethnicity == vec_pair[2]))]
    
    df_design <- seurat_subset[[]]
    df_design$Batch <- sapply(df_design$Library, FUN = function(x){strsplit(x, split = "_L00")[[1]][1]})
    df_design <- df_design[, c("DCP_ID", "Age", "Sex", "Ethnicity", "Batch")]
    df_design <- unique(df_design)
    rownames(df_design) <- df_design$DCP_ID
    df_cellcount <- data.frame(table(seurat_subset[[]]$DCP_ID))
    df_cellcount <- df_cellcount[which(df_cellcount$Freq >= 10), ]
    df_design <- df_design[df_design$DCP_ID %in% df_cellcount$Var1, ]
    df_design$Chinese <- "Not_Chinese"
    df_design$Chinese[which(df_design$Ethnicity == "Chinese")] <- "Chinese"
    df_design$Not_Chinese <- 1
    df_design$Not_Chinese[which(df_design$Ethnicity == "Chinese")] <- 0
    
    df_design$Malay <- "Not_Malay"
    df_design$Malay[which(df_design$Ethnicity == "Malay")] <- "Malay"
    df_design$Not_Malay <- 1
    df_design$Not_Malay[which(df_design$Ethnicity == "Malay")] <- 0
    
    df_design$Indian <- "Not_Indian"
    df_design$Indian[which(df_design$Ethnicity == "Indian")] <- "Indian"
    df_design$Not_Indian <- 1
    df_design$Not_Indian[which(df_design$Ethnicity == "Indian")] <- 0
    
    seurat_subset
    
    ##### Pseudobulk of count matrix from gene-cell to gene-donor #####
    
    matrix_counts <- AggregateExpression(seurat_subset, 
                                         assays = "RNA", 
                                         group.by = "DCP_ID", 
                                         slot = "counts")$RNA
    matrix_counts <- t(matrix_counts)
    matrix_counts <- matrix_counts[rownames(df_design), ]
    
    vec_donors_per_gene <- colSums(matrix_counts > 0)
    vec_UMIs_per_gene <- colSums(matrix_counts)
    
    vec_genes_to_keep <- intersect(names(which(vec_donors_per_gene >= 0.1*dim(df_design)[1])), 
                                   names(which(vec_UMIs_per_gene >= dim(df_design)[1])))
    matrix_counts <- matrix_counts[, vec_genes_to_keep]
    
    y <- DGEList(counts = t(matrix_counts))
    y <- calcNormFactors(y)
    vec_topDEG <- c()
    
    str_ethnicity <- vec_pair[1]
    
    df_design[, str_ethnicity] <- factor(df_design[, str_ethnicity], 
                                         levels=c(paste0("Not_", str_ethnicity), 
                                                  str_ethnicity))
    design <- model.matrix(formula(paste0("~Sex + Age + Batch + ", str_ethnicity)), 
                           data=df_design)
    
    y <- estimateDisp(y, design)
    
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit)
    topTags(lrt)
    vec_topDEG <- c(vec_topDEG, rownames(lrt$table %>% slice_max(logFC, n = 10)))
    write.table(lrt$table, 
                file = paste("edgeR", "pairwise", vec_pair[1], vec_pair[2], str_subtype, ".txt", sep = "_"), 
                col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
  }
}

### Examining DEGs

for (str_pairwise in c("Chinese_Malay", "Chinese_Indian", "Malay_Indian")){
  for (str_celltype in c("atypical_B", "memory_B_IGHMhi", "memory_B_IGHMlo", "naive_B", 
                         "CD14+_Monocyte", "CD16+_Monocyte", "cDC2", 
                         "CD16+_NK", "CD4+_T_cm", "CD4+_T_cyt", "CD4+_T_em", "CD4+_T_naive", 
                         "CD56+_NK", "CD8+_T_GZMBhi", "CD8+_T_GZMKhi", "CD8+_T_naive", 
                         "gdT_GZMBhi", "gdT_GZMKhi", "MAIT", "Treg")){
    
    df_edgeR <- read.table(paste("edgeR", "pairwise", str_pairwise, str_celltype, ".txt", sep = "_"), 
                           header = TRUE, sep = "\t")
    df_edgeR$gene <- row.names(df_edgeR)
    df_edgeR$FDR <- p.adjust(df_edgeR$PValue, method = "fdr")
    df_edgeR <- df_edgeR[which(df_edgeR$FDR < 0.05), ]
    
    if (dim(df_edgeR)[1] > 0){
      
      df_edgeR$comparison <- str_pairwise
      df_edgeR$subtype <- str_celltype
      
      if ((str_pairwise == "Chinese_Malay") & (str_celltype == "atypical_B")){
        df_strong_DEG <- df_edgeR } else {
          df_strong_DEG <- rbind(df_strong_DEG, df_edgeR)
        }
    }
  }
}

##### Concordance analysis with Singapore Integrative Omics whole blood microarray dataset (iOmics) #####

##### USER INPUT: READ LSM TABLE FROM iOmics AS df_SIOS #####

df_SIOS <- read.table(file = "SIOS_Top280ProbeSets.txt", header = TRUE, sep = "\t")
df_SIOS <- df_SIOS[-which(df_SIOS$Gene == "FBXO9"), ]
row.names(df_SIOS) <- df_SIOS$Gene
df_SIOS$log2FC_Chinese_Malay <- df_SIOS$lsm_chinese - df_SIOS$lsm_malay
df_SIOS$log2FC_Malay_Indian <- df_SIOS$lsm_malay - df_SIOS$lsm_indian

df_pseudobulk_Chinese_Indian <- read.table(file = "edgeR_pairwise_pseudobulk_PBMC_Chinese_Indian_.txt", header = TRUE, sep = "\t")
df_pseudobulk_Chinese_Malay <- read.table(file = "edgeR_pairwise_pseudobulk_PBMC_Chinese_Malay_.txt", header = TRUE, sep = "\t")
df_pseudobulk_Malay_Indian <- read.table(file = "edgeR_pairwise_pseudobulk_PBMC_Malay_Indian_.txt", header = TRUE, sep = "\t")

df_pseudobulk_Chinese_Indian$Comparison <- "Chinese_Indian"
df_pseudobulk_Chinese_Malay$Comparison <- "Chinese_Malay"
df_pseudobulk_Malay_Indian$Comparison <- "Malay_Indian"

df_pseudobulk_Chinese_Indian$Gene <- row.names(df_pseudobulk_Chinese_Indian)
df_pseudobulk_Chinese_Malay$Gene <- row.names(df_pseudobulk_Chinese_Malay)
df_pseudobulk_Malay_Indian$Gene <- row.names(df_pseudobulk_Malay_Indian)

df_merge_pseudobulk <- rbind(df_pseudobulk_Chinese_Indian, df_pseudobulk_Chinese_Malay, df_pseudobulk_Malay_Indian)
df_merge_pseudobulk$FDR <- p.adjust(df_merge_pseudobulk$PValue, method = "fdr")

df_SIOS_Chinese_Indian <- df_SIOS[, c("Gene", "log2FC_Chinese_Indian")]
df_SIOS_Chinese_Indian$Comparison <- "Chinese_Indian"
df_SIOS_Chinese_Malay <- df_SIOS[, c("Gene", "log2FC_Chinese_Malay")]
df_SIOS_Chinese_Malay$Comparison <- "Chinese_Malay"
df_SIOS_Malay_Indian <- df_SIOS[, c("Gene", "log2FC_Malay_Indian")]
df_SIOS_Malay_Indian$Comparison <- "Malay_Indian"

colnames(df_SIOS_Chinese_Indian)[grep("log2FC", colnames(df_SIOS_Chinese_Indian))] <- "log2FC_SIOS"
colnames(df_SIOS_Chinese_Malay)[grep("log2FC", colnames(df_SIOS_Chinese_Malay))] <- "log2FC_SIOS"
colnames(df_SIOS_Malay_Indian)[grep("log2FC", colnames(df_SIOS_Malay_Indian))] <- "log2FC_SIOS"

df_merged_SIOS <- rbind(df_SIOS_Chinese_Indian, df_SIOS_Chinese_Malay, df_SIOS_Malay_Indian)

df_merged <- merge(df_merged_SIOS, df_merge_pseudobulk, by = c("Gene", "Comparison"))

library(ggplot2)

ggplot(df_merged, aes(x = log2FC_SIOS, y = logFC, colour = Comparison)) + geom_point(size = 0.8) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  xlab(label = "SIOS log2FC") + ylab("AIDA log2FC") + 
  ggtitle("Comparisons of AIDA log2FC values (all)\nagainst top SIOS probe sets") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 15)) +
  theme(legend.title=element_text(size = 12), 
        legend.text=element_text(size = 12))

cor.test(df_merged$log2FC_SIOS, df_merged$logFC)

length(which((df_merged$log2FC_SIOS > 0) & (df_merged$logFC < 0)))
length(which((df_merged$log2FC_SIOS < 0) & (df_merged$logFC > 0)))

df_merge_005 <- df_merged[(df_merged$FDR < 0.05), ]

ggplot(df_merge_005, aes(log2FC_SIOS, logFC, colour = Comparison)) + geom_point(size = 0.8) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  xlab(label = "SIOS log2FC") + ylab("AIDA log2FC") + 
  ggtitle("Comparisons of AIDA log2FC values (FDR<0.05)\nagainst top SIOS probe sets") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 15)) +
  theme(legend.title=element_text(size = 12), 
        legend.text=element_text(size = 12))

cor.test(df_merge_005$log2FC_SIOS, df_merge_005$logFC)

length(which((df_merge_005$log2FC_SIOS > 0) & (df_merge_005$logFC < 0)))
length(which((df_merge_005$log2FC_SIOS < 0) & (df_merge_005$logFC > 0)))