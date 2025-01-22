### AIDA Phase 1 Data Freeze v1 Colocalisation (RScript) analysis: version 22 January 2025

##### Preparing environment for colocalisation analyses #####

# install.packages("remotes")
# library(remotes)
# install_github("chr1swallace/coloc@main", build_vignettes = TRUE)
# install_github("boxiangliu/locuscomparer")
library(coloc)
library(locuscomparer)
library(ggplot2)

params <- commandArgs(trailingOnly = TRUE)
str_celltype <- params[1]

##### INPUT REQUIRED: READ IN GWAS INFORMATION TABLE AS df_GWAS_metadata #####
##### INPUT REQUIRED: READ IN AIDA ALLELE FREQUENCIES TABLE AS df_eQTL_allelefreq #####
##### INPUT REQUIRED: READ IN 1000 GENOMES ALLELE FREQUENCIES TABLE AS df_1kG #####

df_GWAS_metadata <- read.table("AIDA_Colocalisation_Compiled_GWAS_SummaryStatistics.txt", 
                               sep = "\t", header = TRUE)
rownames(df_GWAS_metadata) <- df_GWAS_metadata$Study
df_GWAS_metadata <- df_GWAS_metadata[which(df_GWAS_metadata$type == "cc"), ]

df_eQTL_allelefreq <- read.table(file = "AIDA_eQTL_variants_allelefreq.txt", 
                                 header = FALSE, sep = "\t")
colnames(df_eQTL_allelefreq) <- c("SNP", "AltAlleleFreq")
rownames(df_eQTL_allelefreq) <- df_eQTL_allelefreq$SNP

df_1kG <- read.table(file = "Formatted_1000GENOMES-phase_3.txt", header = TRUE, sep = "\t")

##### Preparing eQTL files for cell type of interest #####

  df_eQTL <- read.table(file = paste0("/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_cis_output_", 
                                      str_celltype, ".txt"),
                        header = TRUE, sep = "\t")
  df_covariates <- read.table(file = paste0("/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL_", 
                                            str_celltype, "_covariates.txt"),
                              header = TRUE, sep = "\t")
  
  df_SNP <- read.table(paste("/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eigenMT/eQTL", str_celltype, "genotypes.txt", sep = "_"),
                       header = TRUE, row.names = 1)
  
  df_expression <- read.table(file = paste("/mnt/volume1/AIDA_5GEX_DataFreeze_v1/eQTL/eQTL", str_celltype, "normgeneexp.txt", sep = "_"), header = TRUE, sep = "\t", row.names = 1)
  vec_donors_of_interest <- colnames(df_expression)
  
  ##### Looping through each set of GWAS summary statistics #####
  
  for (str_GWAS in rownames(df_GWAS_metadata)){
    
    vec_coloc_celltype <- c()
    vec_coloc_gene <- c()
    vec_coloc_GWAS <- c()
    vec_coloc_PP_H4_abf <- c()
    
    df_GWAS <- read.table(file = paste0("Formatted_", str_GWAS, ".txt"), 
                          header = TRUE, sep = "\t")
    rownames(df_GWAS) <- df_GWAS$SNPID_hg38
    vec_unique_eQTL_gene <- unique(df_eQTL$gene)
    
    ##### Looping through each gene in each set of eQTL summary statistics #####
    
    for (str_gene in vec_unique_eQTL_gene){
      
      df_eQTL_subset <- df_eQTL[which(df_eQTL$gene == str_gene), ]
      # Only set row names here, as multiple SNPs for each gene
      rownames(df_eQTL_subset) <- df_eQTL_subset$SNP
      # Find intersection of loci ID of eQTL with loci ID of a GWAS summary statistic
      vec_intersecting_SNPs <- intersect(rownames(df_eQTL_subset), rownames(df_GWAS))
      
      if (length(vec_intersecting_SNPs) >= 100){
        
        df_GWAS_subset <- df_GWAS[vec_intersecting_SNPs, ]
        df_eQTL_subset <- df_eQTL_subset[vec_intersecting_SNPs, ]
        
        df_eQTL_subset$beta_se <- df_eQTL_subset$beta / df_eQTL_subset$t.stat
        df_eQTL_subset$beta_var <- (df_eQTL_subset$beta_se)^2
        
        df_GWAS_subset$Var <- (df_GWAS_subset$SE)^2
        df_eQTL_subset$pos <- sapply(df_eQTL_subset$SNP, FUN = function(x){strsplit(x, split = "_")[[1]][2]})
        df_GWAS_subset$pos <- sapply(df_GWAS_subset$SNPID_hg38, FUN = function(x){strsplit(x, split = "_")[[1]][2]})
        
        df_eQTL_subset$AltAlleleFreq <- df_eQTL_allelefreq[rownames(df_eQTL_subset), 
                                                           "AltAlleleFreq"]
        
        ##### Objects for colocalisation analyses #####
        
        list_eQTL <- list()
        list_eQTL[["beta"]] <- df_eQTL_subset$beta
        list_eQTL[["varbeta"]] <- df_eQTL_subset$beta_var
        list_eQTL[["pvalues"]] <- df_eQTL_subset$p.value
        list_eQTL[["MAF"]] <- df_eQTL_subset$AltAlleleFreq
        list_eQTL[["snp"]] <- df_eQTL_subset$SNP
        list_eQTL[["N"]] <- (length(colnames(df_covariates)) - 1)
        list_eQTL[["type"]] <- "quant"
        list_eQTL[["position"]] <- df_eQTL_subset$pos
        
        list_GWAS <- list()
        list_GWAS[["beta"]] <- df_GWAS_subset$Beta
        list_GWAS[["varbeta"]] <- df_GWAS_subset$Var
        list_GWAS[["pvalues"]] <- df_GWAS_subset$Pvalue
        # Assumed same AltAlleleFreq as df_eQTL_subset, in case not present in df_GWAS
        list_GWAS[["MAF"]] <- df_eQTL_subset$AltAlleleFreq
        list_GWAS[["snp"]] <- df_GWAS_subset$SNPID_hg38
        list_GWAS[["N"]] <- df_GWAS_metadata[str_GWAS, "N"]
        list_GWAS[["type"]] <- df_GWAS_metadata[str_GWAS, "type"]
        if (df_GWAS_metadata[str_GWAS, "type"] == "cc"){
          list_GWAS[["s"]] <- df_GWAS_metadata[str_GWAS, "s"]
        }
        list_GWAS[["position"]] <- df_GWAS_subset$pos
        
        my.res <- coloc.abf(dataset1 = list_eQTL,
                            dataset2 = list_GWAS)
        
        vec_coloc_celltype <- c(vec_coloc_celltype, str_celltype)
        vec_coloc_gene <- c(vec_coloc_gene, str_gene)
        vec_coloc_GWAS <- c(vec_coloc_GWAS, str_GWAS)
        vec_coloc_PP_H4_abf <- c(vec_coloc_PP_H4_abf, 
                                 as.numeric(my.res$summary["PP.H4.abf"]))
        
        if (as.numeric(my.res$summary["PP.H4.abf"]) > 0.9){
          print(str_celltype)
          print(str_gene)
          print(str_GWAS)
          print(as.numeric(my.res$summary["PP.H4.abf"]))
        }
      }
    }
    
    df_colocalisation_summary <- data.frame(CellType = vec_coloc_celltype, 
                                            Gene = vec_coloc_gene, 
                                            GWAS = vec_coloc_GWAS, 
                                            coloc_PP_H4_abf = vec_coloc_PP_H4_abf)
    
    write.table(df_colocalisation_summary, 
                file = paste("df_colocalisation", str_celltype, str_GWAS, 
                             "_summary.txt", sep = "_"),  
                col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
##### After completing colocalisation for a particular GWAS summary statistic, plot examples that had coloc_PP_H4_abf > 0.9 #####

df_priority <- df_colocalisation_summary[which(df_colocalisation_summary$coloc_PP_H4_abf > 0.9), ]
vec_AA <- c()
vec_AB <- c()
vec_BB <- c()
vec_Wilcoxon <- c()

##### Looping through each prioritised gene for this set of GWAS summary statistics #####

for (str_gene in unique(df_priority$Gene)){
  
  df_eQTL_subset <- df_eQTL[which(df_eQTL$gene == str_gene), ]
  # Only set row names here, as multiple SNPs to each gene
  rownames(df_eQTL_subset) <- df_eQTL_subset$SNP
  # Find intersection of loci ID of eQTL with loci ID of a GWAS summary statistic
  vec_intersecting_SNPs <- intersect(rownames(df_eQTL_subset), rownames(df_GWAS))
  
  ### Convert SNP_IDs to rsIDs
  
  vec_intersecting_SNPs <- vec_intersecting_SNPs[vec_intersecting_SNPs %in% df_1kG$str_SNPID]
  
  df_GWAS_subset <- df_GWAS[vec_intersecting_SNPs, ]
  df_GWAS_subset$rsid <- df_1kG$str_rsID[match(vec_intersecting_SNPs, df_1kG$str_SNPID)]
  df_eQTL_subset <- df_eQTL_subset[vec_intersecting_SNPs, ]
  df_eQTL_subset$rsid <- df_GWAS_subset$rsid
  
  ### Need to format names of columns
  
  colnames(df_GWAS_subset)[which(colnames(df_GWAS_subset) == "Pvalue")] <- "pval"
  colnames(df_eQTL_subset)[which(colnames(df_eQTL_subset) == "p.value")] <- "pval"
  filename_df_GWAS_subset <- paste("df_GWAS_subset", 
                                   str_celltype, str_gene, str_GWAS, ".txt", sep = "_")
  write.table(df_GWAS_subset[, c("rsid", "pval")], 
              file = filename_df_GWAS_subset, 
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  filename_df_eQTL_subset <- paste("df_eQTL_subset", 
                                   str_celltype, str_gene, str_GWAS, ".txt", sep = "_")
  write.table(df_eQTL_subset[, c("rsid", "pval")], 
              file = filename_df_eQTL_subset, 
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  pdf(paste("Rplot_locuscomparer", 
            str_celltype, str_gene, str_GWAS, ".pdf", sep = "_"), 
      width = 10, height = 10)
  print(locuscompare(in_fn1 = filename_df_GWAS_subset, 
                     in_fn2 = filename_df_eQTL_subset, 
                     title1 = str_GWAS, title2 = paste0(str_gene, "_eQTL"), 
                     population = "EAS", genome = "hg38"))
  dev.off()
  
  ##### Get information about lead SNP #####
  
  d1 = read_metal(filename_df_GWAS_subset)
  d2 = read_metal(filename_df_eQTL_subset)
  merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
  print(df_1kG[which(df_1kG$str_rsID == get_lead_snp(merged)), ])
  
  ##### Plotting cis-eQTLs effects: histogram and QQplot of gene expression, cis-eQTL plot #####
  
  str_SNP <- df_1kG[which(df_1kG$str_rsID == get_lead_snp(merged)), ]$str_SNPID
  str_genotype_A <- strsplit(str_SNP, split = "_")[[1]][3]
  str_genotype_B <- strsplit(str_SNP, split = "_")[[1]][4]
  
  png(paste("Rplot_Histogram", 
            str_celltype, str_gene, str_GWAS, ".png", sep = "_"), 
      width = 2000, height = 1000, units = "px")
  hist((as.numeric(df_expression[str_gene,])), breaks = 20)
  dev.off()
  
  vec_to_plot <- as.numeric(df_expression[str_gene,])
  
  png(paste("Rplot_QQplot", 
            str_celltype, str_gene, str_GWAS, ".png", sep = "_"), 
      width = 2000, height = 1000, units = "px")
  print(qqnorm(vec_to_plot))
  print(qqline(vec_to_plot))
  dev.off()
  
  df_plot <- data.frame("Genotype" = as.character(df_SNP[str_SNP, c(vec_donors_of_interest)]),
                        "Expression" = as.numeric(df_expression[str_gene, c(vec_donors_of_interest)]))
  df_plot[which(df_plot$Genotype == "0"), "Genotype"] <- paste(str_genotype_A, str_genotype_A, sep = "")
  df_plot[which(df_plot$Genotype == "1"), "Genotype"] <- paste(str_genotype_A, str_genotype_B, sep = "")
  df_plot[which(df_plot$Genotype == "2"), "Genotype"] <- paste(str_genotype_B, str_genotype_B, sep = "")
  df_plot$Genotype <- factor(df_plot$Genotype, 
                             levels = c(paste(str_genotype_A, str_genotype_A, sep = ""), 
                                        paste(str_genotype_A, str_genotype_B, sep = ""), 
                                        paste(str_genotype_B, str_genotype_B, sep = "")))
  
  ggplot(df_plot, aes(Genotype, Expression)) +
    geom_boxplot() +
    ggtitle(paste(str_SNP, ", a cis-eQTL for ", str_gene, sep = ""))
  ggsave(paste("Rplot_eQTL", 
               str_celltype, str_gene, str_GWAS, ".png", sep = "_"))
  
  ##### Save information about numbers of donors per genotype, and Wilcoxon rank-sum p-value #####
  
  vec_AA <- c(vec_AA, length(which(df_plot$Genotype == paste(str_genotype_A, str_genotype_A, sep = ""))))
  vec_AB <- c(vec_AB, length(which(df_plot$Genotype == paste(str_genotype_A, str_genotype_B, sep = ""))))
  vec_BB <- c(vec_BB, length(which(df_plot$Genotype == paste(str_genotype_B, str_genotype_B, sep = ""))))
  vec_Wilcoxon <- c(vec_Wilcoxon, wilcox.test((df_plot$Expression[which(df_plot$Genotype == paste(str_genotype_A, str_genotype_A, sep = ""))]), 
                                              (df_plot$Expression[which(df_plot$Genotype == paste(str_genotype_B, str_genotype_B, sep = ""))]))$p.value)
  
  if (match(str_gene, unique(df_priority$Gene)) == 1){
    df_1kG_allele_information <- df_1kG[which(df_1kG$str_rsID == get_lead_snp(merged)), ]
  } else {
    df_1kG_allele_information <- rbind(df_1kG_allele_information, 
                                       df_1kG[which(df_1kG$str_rsID == get_lead_snp(merged)), ])
  }
}

##### Save variant-gene-allele frequency information #####

df_1kG_allele_information$Gene <- unique(df_priority$Gene)
df_1kG_allele_information$vec_AA <- vec_AA
df_1kG_allele_information$vec_AB <- vec_AB
df_1kG_allele_information$vec_BB <- vec_BB
df_1kG_allele_information$vec_Wilcoxon <- vec_Wilcoxon

write.table(df_1kG_allele_information, 
            file = paste("df_1kG_allele_information", 
                         str_celltype, str_gene, str_GWAS, ".txt", sep = "_"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  }