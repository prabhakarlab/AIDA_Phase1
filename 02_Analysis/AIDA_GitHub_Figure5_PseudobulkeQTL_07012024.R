### AIDA Phase 1 Data Freeze v1 pseudobulk eQTL analysis: version 7 January 2024

##### INPUT REQUIRED: READ IN COVARIATES TABLE AS df_merged_covariates AND METADATA AS df_AIDA_metadata #####
##### INPUT REQUIRED: READ IN SEURAT OBJECT RDS RELEVANT TO CELL POPULATION OF INTEREST AS seurat_annotated #####
##### INPUT REQUIRED: READ IN CELL ANNOTATIONS AS df_cellannotations #####

### Wrangle covariates (10 gene expression PCs)

library(Seurat)
vec_barcodes_all <- Cells(seurat_annotated)

rownames(df_cellannotations) <- df_cellannotations$barcode

vec_barcodes_flagged <- 
  df_cellannotations$barcode[which(df_cellannotations$Annotation_Level3 == "flagged_cluster")]

vec_barcodes_flagged <- c(vec_barcodes_flagged,
                          df_cellannotations$barcode[which(df_cellannotations$Annotation_Level3 == "flagged_platelet_sum")])

vec_barcodes_flagged <- c(vec_barcodes_flagged,
                          df_cellannotations$barcode[which(df_cellannotations$Annotation_Level3 == "flagged_RBC_UMI")])

### Removing cells from flagged clusters

vec_barcodes_to_keep <- vec_barcodes_all[-which(vec_barcodes_all %in% vec_barcodes_flagged)]
seurat_annotated_to_keep <- seurat_annotated[,vec_barcodes_to_keep]
rm(seurat_annotated)

rownames(df_AIDA_metadata) <- df_AIDA_metadata$barcode
df_AIDA_metadata_filtered_Asian <- df_AIDA_metadata[which(df_AIDA_metadata$DCP_ID %in% rownames(df_merged_covariates)),]

vec_barcodes_to_keep_infilteredonors <- Cells(seurat_annotated_to_keep)[which(Cells(seurat_annotated_to_keep) %in% df_AIDA_metadata_filtered_Asian$barcode)]

seurat_annotated_to_keep <- subset(seurat_annotated_to_keep, 
                                   cells = vec_barcodes_to_keep_infilteredonors)

seurat_annotated_to_keep[["DCP_ID"]] <- df_AIDA_metadata_filtered_Asian[Cells(seurat_annotated_to_keep),
                                                                        "DCP_ID"]
seurat_annotated_to_keep[["Annotation_Level3"]] <- df_cellannotations[Cells(seurat_annotated_to_keep),
                                                                      "Annotation_Level3"]
seurat_annotated_to_keep[["Annotation_Level2"]] <- df_cellannotations[Cells(seurat_annotated_to_keep),
                                                                      "Annotation_Level2"]

DefaultAssay(seurat_annotated_to_keep) <- "RNA"
unique(seurat_annotated_to_keep[["Annotation_Level3"]])$Annotation_Level3
unique(seurat_annotated_to_keep[["Annotation_Level2"]])$Annotation_Level2

###### USER INPUT: Set str_subtype to determine which subtype to find eQTLs for
### Loop through all vec_celltypes_of_interest

### str_subtype = unique(seurat_annotated_to_keep[["Annotation_Level2"]])$Annotation_Level2[2]

vec_celltypes_of_interest = c("IGHMhi_memory_B", "IGHMlo_memory_B", "atypical_B")

for (str_subtype in vec_celltypes_of_interest){
  seurat_annotated_to_keep_subtype <- subset(seurat_annotated_to_keep, 
                                             subset = Annotation_Level2 == str_subtype)
  
  ### Implementing eQTL pipeline
  
  library(dplyr)
  library("MatrixEQTL")
  base.dir = find.package("MatrixEQTL")
  useModel = modelLINEAR
  
  ### Which genes (rows) in 
  ### GetAssayData(object = seurat_annotated_to_keep, slot = "counts")
  ### have more than 1% of cells expressing the gene of interest?
  ### Keep these genes for QTL analyses
  
  matrix_AIDA_Asian <- GetAssayData(object = seurat_annotated_to_keep_subtype, slot = "counts")
  vec_genes_numberofcells <- apply(matrix_AIDA_Asian, MARGIN = 1, function (x) length(which(x > 0)))
  vec_genes_to_keep <- vec_genes_numberofcells[vec_genes_numberofcells > dim(matrix_AIDA_Asian)[2]*0.01]
  
  ### For donors with >=10 cells: 
  ### assigning vector of normalised gene expression to object with donor_ID
  ### Normalisation using total UMI count within cell
  
  df_donors_of_interest <- as.data.frame(table(seurat_annotated_to_keep_subtype[["DCP_ID"]]))
  colnames(df_donors_of_interest) <- c("DCP_ID", "Freq")
  vec_donors_of_interest <- as.character(df_donors_of_interest$DCP_ID[which(df_donors_of_interest$Freq >= 10)])
  
  ### Removing related donor samples
  
  vec_donors_to_remove <- c()
  
  ### For related pairs, check if match length == 2; if so, keep donor with the higher genotyping call rate
  
  list_related_match_check <- list()
  
  list_related_match_check[[1]] <- match(c("SG_HEL_H180", "SG_HEL_H07a"), vec_donors_of_interest)
  list_related_match_check[[2]] <- match(c("SG_HEL_H152", "SG_HEL_H09a"), vec_donors_of_interest)
  list_related_match_check[[3]] <- match(c("SG_HEL_H185", "SG_HEL_H217"), vec_donors_of_interest)
  list_related_match_check[[4]] <- match(c("SG_HEL_H276", "SG_HEL_H313"), vec_donors_of_interest)
  list_related_match_check[[5]] <- match(c("SG_HEL_H230", "SG_HEL_H234"), vec_donors_of_interest)
  list_related_match_check[[6]] <- match(c("SG_HEL_H135", "SG_HEL_H134"), vec_donors_of_interest)
  list_related_match_check[[7]] <- match(c("KR_SGI_H003", "KR_SGI_H002"), vec_donors_of_interest)
  list_related_match_check[[8]] <- match(c("KR_SGI_H046", "KR_SGI_H009"), vec_donors_of_interest)
  list_related_match_check[[9]] <- match(c("KR_SGI_H078", "KR_SGI_H017"), vec_donors_of_interest)
  list_related_match_check[[10]] <- match(c("KR_SGI_H025", "KR_SGI_H026"), vec_donors_of_interest)
  list_related_match_check[[11]] <- match(c("KR_SGI_H056", "KR_SGI_H057"), vec_donors_of_interest)
  list_related_match_check[[12]] <- match(c("KR_SGI_H071", "KR_SGI_H069"), vec_donors_of_interest)
  list_related_match_check[[13]] <- match(c("KR_SGI_H074", "KR_SGI_H158"), vec_donors_of_interest)
  list_related_match_check[[14]] <- match(c("KR_SGI_H081", "KR_SGI_H079"), vec_donors_of_interest)
  list_related_match_check[[15]] <- match(c("KR_SGI_H095", "KR_SGI_H093"), vec_donors_of_interest)
  
  for (i in 1:15){
    vec_match_check <- list_related_match_check[[i]]
    if (sum(!is.na(vec_match_check)) == 2){
      vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[2])
    }
  }
  
  ### For related trios, check if match length > 1; if so
  ### Check if match length == 3; if so, keep the first donor noted above
  ### Check if match length == 2; if so, check through all three possible pairs to keep higher genotyping call rate donor
  ### In terms of genotyping call rates: "SG_HEL_H379" > "SG_HEL_H378" > "SG_HEL_H380"
  ### In terms of genotyping call rates: "KR_SGI_H163" > "KR_SGI_H162" > "KR_SGI_H161"
  
  vec_match_check <- match(c("SG_HEL_H380", "SG_HEL_H378", "SG_HEL_H379"), vec_donors_of_interest)
  if (sum(!is.na(vec_match_check)) > 1){
    if (sum(!is.na(vec_match_check)) == 3){
      vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[1], vec_match_check[2])
    }
    if (sum(!is.na(vec_match_check)) == 2){
      if (sum(!is.na(match(c("SG_HEL_H380", "SG_HEL_H378"), vec_donors_of_interest))) == 2){
        vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[1])
      }
      if (sum(!is.na(match(c("SG_HEL_H380", "SG_HEL_H379"), vec_donors_of_interest))) == 2){
        vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[1])
      }
      if (sum(!is.na(match(c("SG_HEL_H378", "SG_HEL_H379"), vec_donors_of_interest))) == 2){
        vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[2])
      }
    }
  }
  
  vec_match_check <- match(c("KR_SGI_H161", "KR_SGI_H162", "KR_SGI_H163"), vec_donors_of_interest)
  if (sum(!is.na(vec_match_check)) > 1){
    if (sum(!is.na(vec_match_check)) == 3){
      vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[1], vec_match_check[2])
    }
    if (sum(!is.na(vec_match_check)) == 2){
      if (sum(!is.na(match(c("KR_SGI_H161", "KR_SGI_H162"), vec_donors_of_interest))) == 2){
        vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[1])
      }
      if (sum(!is.na(match(c("KR_SGI_H161", "KR_SGI_H163"), vec_donors_of_interest))) == 2){
        vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[1])
      }
      if (sum(!is.na(match(c("KR_SGI_H162", "KR_SGI_H163"), vec_donors_of_interest))) == 2){
        vec_donors_to_remove <- c(vec_donors_to_remove, vec_match_check[2])
      }
    }
  }
  
  ### Retaining non-related donors
  
  vec_donors_of_interest <- vec_donors_of_interest[-vec_donors_to_remove]
  
  for (donor in vec_donors_of_interest){
    print(donor)
    eQTL_subset <- subset(x = seurat_annotated_to_keep_subtype, subset = DCP_ID == donor)
    eQTL_matrix <- GetAssayData(object = eQTL_subset, slot = "counts")
    eQTL_matrix <- eQTL_matrix[names(vec_genes_to_keep),]
    eQTL_normalised_matrix <- sweep(eQTL_matrix, MARGIN = 2, colSums(eQTL_matrix), "/")
    
    ### Mean of all cells in a single donor
    eQTL_vec_per_donor <- apply(eQTL_normalised_matrix, MARGIN = 1, mean)
    ### Log1p transformation of cell means in each donor
    assign(donor, log1p(eQTL_vec_per_donor*10000))
    
  }
  
  ### Getting normalised gene expression data
  df_normalised_gene_exp <- data.frame(Gene = names(vec_genes_to_keep))
  
  ### Putting together a data frame of the normalised data from donors of interest 
  for (donor in vec_donors_of_interest){
    df_normalised_gene_exp <- cbind(df_normalised_gene_exp, get(donor))
  }
  
  colnames(df_normalised_gene_exp) <- c("Gene", as.character(vec_donors_of_interest))
  rownames(df_normalised_gene_exp) <- df_normalised_gene_exp$Gene
  ### Checking distribution of log1p(values)
  ### Log1p values look more like normal distribution than original matrix
  # hist(as.numeric(df_normalised_gene_exp[,c(2:dim(df_normalised_gene_exp)[2])][1,]), breaks = 20)
  # hist(log1p(as.numeric(df_normalised_gene_exp[,c(2:dim(df_normalised_gene_exp)[2])][1,])*10000))
  
  ### Remove AC240274.1, AC004556.3, and AL592183.1 from df_normalised_gene_exp
  df_normalised_gene_exp <- df_normalised_gene_exp[-c(which(df_normalised_gene_exp$Gene == "AC240274.1"),
                                                      which(df_normalised_gene_exp$Gene == "AC004556.3"),
                                                      which(df_normalised_gene_exp$Gene == "AL592183.1")), ]
  
  ### Adjust colname(1) of df_normalised_gene_exp
  colnames(df_normalised_gene_exp)[1] <- "DCP_ID"
  filename_df_normalised_gene_exp <-
    paste("/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_", str_subtype, "_normgeneexp.txt", sep = "")
  write.table(df_normalised_gene_exp, file = filename_df_normalised_gene_exp, sep = "\t", quote = FALSE, row.names = FALSE)
  
  ### Checking for normality
  # qqnorm(df_normalised_gene_exp[1,c(2:dim(df_normalised_gene_exp)[2])], pch = 1, frame = FALSE)
  # qqline(df_normalised_gene_exp[1,c(2:dim(df_normalised_gene_exp)[2])], col = "steelblue", lwd = 2)
  
  ### Perform principal component analysis of gene expression values
  
  PCA_df_normalised_gene_exp <- prcomp(df_normalised_gene_exp[,c(2:dim(df_normalised_gene_exp)[2])])
  df_PCA_rotation_normalised_gene_exp <- as.data.frame(PCA_df_normalised_gene_exp$rotation[,c(1:10)])
  colnames(df_PCA_rotation_normalised_gene_exp) <- c("GeneExp_PC1",
                                                     "GeneExp_PC2",
                                                     "GeneExp_PC3",
                                                     "GeneExp_PC4",
                                                     "GeneExp_PC5",
                                                     "GeneExp_PC6",
                                                     "GeneExp_PC7",
                                                     "GeneExp_PC8",
                                                     "GeneExp_PC9",
                                                     "GeneExp_PC10")
  
  ### Wrangle covariates data frame, and 
  ### sort according to gene expression file DCP_ID column order
  
  df_covariates <- merge(df_PCA_rotation_normalised_gene_exp, df_merged_covariates, 
                         by = "row.names")
  rownames(df_covariates) <- df_covariates$Row.names
  df_covariates$Row.names <- NULL
  
  ### Checking distribution of ancestries, and that df_expression column names
  ### are in same order as df_covariate row names
  rownames(df_covariates) == colnames(df_normalised_gene_exp[2:dim(df_normalised_gene_exp)[2]])
  ### All true
  
  ### In case: write gene names to a file
  # filename_genename <-
  #   paste("/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_", str_subtype, "_genename.txt", sep = "")
  # write.table(as.data.frame(df_normalised_gene_exp[,1]),
  #             file = filename_genename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ### Need to transpose covariates file and convert strings to categorical numerical variables
  ### Write covariates data frame to file
  
  df_transposed_covariates <- t(df_covariates)
  filename_covariates <-
    paste("/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_", str_subtype, "_covariates.txt", sep = "")
  ### Used col.names = NA to include a "id + \t" at top-left corner
  write.table(df_transposed_covariates, 
              file = filename_covariates, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  
  covariates_file_name = filename_covariates
  
  ##### INPUT REQUIRED: LOCATION OF TAB-DELIMITED SNPS FILE AS SNP_file_name #####
  ##### INPUT REQUIRED: LOCATION OF POSITIONS OF SNPS FILE AS eQTL_df_snpos #####
  ##### INPUT REQUIRED: LOCATION OF POSITIONS OF GENES FILE AS eQTL_df_genepos #####
  
  SNP_file_name = "/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_VCF_allchr_filtered.txt"
  
  expression_file_name = filename_df_normalised_gene_exp
  
  eQTL_df_snpos <- read.table(file = "/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_snppos_for_matrixQTL_allchr_filtered.txt",
                              sep = "\t", header = TRUE)
  
  ### Wrote master genepos file for re-arranging in R, for input into MatrixeQTL
  ### To re-arrange rows using df_normalised_gene_exp[,1] vector:
  
  eQTL_df_genepos <- read.table(file = "/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_AIDA_GENCODEv32_genepos.txt",
                                sep = "\t", header = TRUE)
  rownames(eQTL_df_genepos) <- eQTL_df_genepos$Gene
  eQTL_df_genepos <- eQTL_df_genepos[df_normalised_gene_exp[,1],]
  
  filename_genepos <- paste("/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_", str_subtype, "_genepos.txt", sep = "")
  write.table(eQTL_df_genepos, file = filename_genepos, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  eQTL_df_genepos <- read.table(file = filename_genepos,
                                sep = "\t", header = TRUE)
  
  output_file_name_cis = paste("eQTL_cis_output_", str_subtype, ".txt", sep = "")
  output_file_name_tra = paste("eQTL_trans_output_", str_subtype, ".txt", sep = "")
  
  ### Set as 1 to record all cis-eQTL p-values
  pvOutputThreshold_cis = 1
  ### Much higher p-value threshold for trans-eQTLs
  pvOutputThreshold_tra = 1e-9
  errorCovariance = numeric()
  
  cisDist = 1e6
  
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"      # the TAB character
  snps$fileOmitCharacters = "NA" # denote missing values;
  snps$fileSkipRows = 1          # one row of column labels
  snps$fileSkipColumns = 1       # one column of row labels
  snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  snps$LoadFile( SNP_file_name )
  
  vec_order_donors_in_snps <- c()
  for (donor in as.character(vec_donors_of_interest)){
    vec_order_donors_in_snps <- c(vec_order_donors_in_snps, which(colnames(snps) == donor))
  }
  
  snps$ColumnSubsample(vec_order_donors_in_snps)
  
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"      # the TAB character
  gene$fileOmitCharacters = "NA" # denote missing values;
  gene$fileSkipRows = 1          # one row of column labels
  gene$fileSkipColumns = 1       # one column of row labels
  gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  gene$LoadFile( expression_file_name )
  
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = "\t"      # the TAB character
  cvrt$fileOmitCharacters = "NA" # denote missing values;
  cvrt$fileSkipRows = 1          # one row of column labels
  cvrt$fileSkipColumns = 1       # one column of row labels
  cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  cvrt$LoadFile( covariates_file_name )
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = eQTL_df_snpos,
    genepos = eQTL_df_genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);
  
  ## Plot the Q-Q plot of local and distant p-values
  
  png(filename = paste("eQTL_QQplot_", str_subtype, ".png", sep = ""))
  plot(me)
  dev.off()
}
