### AIDA Phase 1 Data Freeze v2 Cell Type Proportions analysis: version 7 January 2024

library(ggplot2)

##### INPUT REQUIRED: READ IN FULL ANNOTATION TABLE AND METADATA AS df_AIDA_metadata #####
  
### Filtering donors for cell counts >= 800
  
df_AIDA_metadata$count <- 1
rownames(df_AIDA_metadata) <- df_AIDA_metadata$barcode_name
  
### Consider only donors with at least 800 cells and population groups with at least 50 donors
  
df_total_cell_count <- df_AIDA_metadata[, c("DCP_ID", "count")]
df_total_cell_count <- aggregate(count ~ ., df_total_cell_count, FUN = sum)
rownames(df_total_cell_count) <- df_total_cell_count$DCP_ID
vec_donors_to_exclude <- df_total_cell_count$DCP_ID[which(df_total_cell_count$count < 800)]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$DCP_ID %in% vec_donors_to_exclude), ]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$Country == "IN"), ]

##### Annotation Level 1 against all PBMCs #####

vec_unique_celltypes <- unique(df_AIDA_metadata$Annotation_Level1)

### Figure Panels for cell type proportions

for (i in 1:length(vec_unique_celltypes)) {
  str_celltype <- vec_unique_celltypes[i]
  vec_rows_celltype <- which(df_AIDA_metadata$Annotation_Level1 == str_celltype)
  df_celltype <- df_AIDA_metadata[vec_rows_celltype, c("DCP_ID", "count")]
  df_celltype <- aggregate(count ~ ., df_celltype, FUN = sum)
  colnames(df_celltype) <- c("DCP_ID", vec_unique_celltypes[i])
  rownames(df_celltype) <- df_celltype$DCP_ID
  df_celltype$count <- df_total_cell_count[rownames(df_celltype), "count"]
  
  df_celltype_metadata <- df_AIDA_metadata[match(df_celltype$DCP_ID, df_AIDA_metadata$DCP_ID),
                                           c("Genotyping_ID", "Library", "Country",
                                             "DCP_ID", "Age", "Sex", "Ancestry")]
  rownames(df_celltype_metadata) <- df_celltype_metadata$DCP_ID
  
  df_merged <- merge(df_celltype_metadata, df_celltype, by = "DCP_ID")
  df_merged$Age <- as.numeric(df_merged$Age)
  df_merged$Proportion <- df_merged[, vec_unique_celltypes[i]]/df_merged$count
  
  df_merged_without_European <- df_merged[-which(df_merged$Ancestry == "European"),]
  df_merged_without_European$Ancestry <- factor(df_merged_without_European$Ancestry, 
                                                levels = c("SG_Chinese", "SG_Indian", "SG_Malay", 
                                                           "Japanese", "Korean", "Thai"))
  
  for (str_ancestry in unique(df_merged_without_European$Ancestry)){
    df_merged_without_European[ , str_ancestry] <- 0
    df_merged_without_European[ , str_ancestry][which(df_merged_without_European$Ancestry == str_ancestry)] <- 1
  }
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Ancestry, colour = Ancestry)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "in each ancestry", collapse = " ")) + 
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 10)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 10))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "ancestry.pdf", sep = "_"), 
         height = 4, width = 5, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Sex)) + 
    geom_boxplot() +
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin females versus males", collapse = " ")) + 
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "sex.pdf", sep = "_"), 
         height = 4, width = 3, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Ancestry, fill = Sex)) + 
    geom_boxplot() +
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby sex and ancestry", collapse = " ")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 10))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "sex_ancestry.pdf", sep = "_"))
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age)) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    ggtitle(paste0("Log10(Proportion) of ", str_celltype, " by age")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 10))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "age.pdf", sep = "_"), 
         height = 4, width = 3, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age, colour = Ancestry)) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby age, trend by ancestry", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 10))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "age_by_Ancestry.pdf", sep = "_"))
  
  pdf(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "sex_ancestry_interaction.pdf", sep = "_"),
      width = 700, height = 700)
  par(mar=c(5.1, 6.1, 4.1, 2.1))
  interaction.plot(x.factor = df_merged_without_European$Sex,
                   trace.factor = df_merged_without_European$Ancestry,
                   response = log10(df_merged_without_European$Proportion),
                   fun = median,
                   col = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3"),#c("red","green", "cyan", "blue", "pink"),
                   xlab = "Sex",
                   ylab = paste("Median log10(proportion) \nof", str_celltype),
                   trace.label = "Ancestry",
                   lwd = 4,
                   lty = 5,
                   cex.lab = 1.2,
                   cex.axis = 1.2)
  title(main = paste("Interaction plot of Sex and Ancestry on median log10(proportion) of",
                     str_celltype))
  dev.off()

  test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion", str_celltype, ".txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + Japanese + SG_Chinese + SG_Malay + SG_Indian + Thai, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion_IndividualAncestryvsKorean", str_celltype, ".txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + Korean + Thai + SG_Chinese + SG_Malay + SG_Indian, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion_IndividualAncestryvsJapanese", str_celltype, ".txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
  
  for (str_individual_ancestry in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
    test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", str_individual_ancestry)), data = df_merged_without_European)
    
    sink(file = paste("Level1/LinearModel_Figure_Log10Proportion_IndividualAncestry_SG_Chinese", str_celltype, str_individual_ancestry, ".txt", sep = "_"))
    print(str_celltype)
    print(str_individual_ancestry)
    print(summary(test_model))
    sink()
  }
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry + Sex:Ancestry + Age:Ancestry, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion", str_celltype, "InteractionsAncestry.txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry + Sex:Age + Sex:Ancestry + Age:Ancestry, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion", str_celltype, "AllTwoWayInteractions.txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
}  

##### Annotation Level 3 against all PBMCs #####
{
  vec_unique_celltypes <- unique(df_AIDA_metadata$Annotation_Level3)
  
  ### Figure Panel for cell type proportions
  
  library(ggplot2)
  vec_individual_covariates <- c("Sex", "Age", "SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")
  vec_model <- c()
  vec_covariate <- c()
  vec_subtype <- c()
  vec_df <- c()
  vec_N <- c()
  vec_Rsquared <- c()
  vec_coefficient <- c()
  vec_tvalue <- c()
  vec_prtvalue <- c()
  
  vec_allancestries_individual_covariates <- c("Sex", "Age", "Ancestry")
  vec_allancestries_model <- c()
  vec_allancestries_covariate <- c()
  vec_allancestries_subtype <- c()
  vec_allancestries_df <- c()
  vec_allancestries_N <- c()
  vec_allancestries_Rsquared <- c()
  vec_allancestries_coefficient <- c()
  vec_allancestries_tvalue <- c()
  vec_allancestries_prtvalue <- c()
  
  vec_ancestry_model <- c()
  vec_ancestry_ancestry <- c()
  vec_ancestry_subtype <- c()
  vec_ancestry_df <- c()
  vec_ancestry_N <- c()
  vec_ancestry_Rsquared <- c()
  vec_ancestry_coefficient <- c()
  vec_ancestry_tvalue <- c()
  vec_ancestry_prtvalue <- c()
}

for (i in 1:length(vec_unique_celltypes)) {
  str_celltype <- vec_unique_celltypes[i]
  vec_rows_celltype <- which(df_AIDA_metadata$Annotation_Level3 == str_celltype)
  df_celltype <- df_AIDA_metadata[vec_rows_celltype, c("DCP_ID", "count")]
  df_celltype <- aggregate(count ~ ., df_celltype, FUN = sum)
  colnames(df_celltype) <- c("DCP_ID", vec_unique_celltypes[i])
  rownames(df_celltype) <- df_celltype$DCP_ID
  df_celltype$count <- df_total_cell_count[rownames(df_celltype), "count"]
  
  df_celltype_metadata <- df_AIDA_metadata[match(df_celltype$DCP_ID, df_AIDA_metadata$DCP_ID),
                                           c("Genotyping_ID", "Library", "Country",
                                             "DCP_ID", "Age", "Sex", "Ancestry")]
  rownames(df_celltype_metadata) <- df_celltype_metadata$DCP_ID
  
  df_merged <- merge(df_celltype_metadata, df_celltype, by = "DCP_ID")
  df_merged$Age <- as.numeric(df_merged$Age)
  df_merged$Proportion <- df_merged[, vec_unique_celltypes[i]]/df_merged$count
  
  df_merged_without_European <- df_merged[-which(df_merged$Ancestry == "European"),]
  df_merged_without_European$Ancestry <- factor(df_merged_without_European$Ancestry, 
                                                levels = c("SG_Chinese", "SG_Indian", "SG_Malay", 
                                                           "Japanese", "Korean", "Thai"))
  
  for (str_ancestry in unique(df_merged_without_European$Ancestry)){
    df_merged_without_European[ , str_ancestry] <- 0
    df_merged_without_European[ , str_ancestry][which(df_merged_without_European$Ancestry == str_ancestry)] <- 1
  }
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Ancestry, colour = Ancestry)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin each ancestry", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "ancestry.pdf", sep = "_"), 
         height = 4, width = 6, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Sex)) +
    geom_boxplot() +
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin each sex", collapse = " ")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15))

  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex.pdf", sep = "_"))

  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Ancestry, fill = Sex)) + 
    geom_boxplot() +
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby sex and ancestry", collapse = " ")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex_ancestry.pdf", sep = "_"), 
         height = 4, width = 6, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age)) +
    geom_point() +
    geom_smooth(method = "lm") +
    ggtitle(paste("Log10(Proportion) of", str_celltype, "by age", collapse = " ")) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 11)) +
    theme(legend.title=element_text(size = 12),
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "age.pdf", sep = "_"), 
         height = 4, width = 3, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age, colour = Ancestry)) +
    geom_point() +
    geom_smooth(method = "lm") + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby age, trend by ancestry", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "age_by_Ancestry.pdf", sep = "_"), 
         height = 4, width = 5, dpi = 300, units = "in")
  
  pdf(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex_ancestry_interaction.pdf", sep = "_"),
      width = 700, height = 700)
  par(mar=c(5.1, 6.1, 4.1, 2.1))
  interaction.plot(x.factor = df_merged_without_European$Sex,
                   trace.factor = df_merged_without_European$Ancestry,
                   response = log10(df_merged_without_European$Proportion),
                   fun = median,
                   col = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3"),#c("red","green", "cyan", "blue", "pink"),
                   xlab = "Sex",
                   ylab = paste("Median log10(proportion) \nof", str_celltype),
                   trace.label = "Ancestry",
                   lwd = 4,
                   lty = 5,
                   cex.lab = 1.2,
                   cex.axis = 1.2)
  title(main = paste("Interaction plot of Sex and Ancestry on median \nlog10(proportion) of",
                     str_celltype))
  dev.off()

  {
    for (str_individual_covariate in vec_individual_covariates){
      test_model <- lm(formula(paste0("log10(Proportion) ~ ", str_individual_covariate)), data = df_merged_without_European)
      
      sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualCovariate", str_celltype, str_individual_covariate, ".txt", sep = "_"))
      print(str_celltype)
      print(str_individual_covariate)
      print(summary(test_model))
      sink()
      vec_model <- c(vec_model, as.character(summary(test_model)$call)[2])
      vec_covariate <- c(vec_covariate, str_individual_covariate)
      vec_subtype <- c(vec_subtype, str_celltype)
      vec_df <- c(vec_df, test_model$df.residual)
      vec_N <- c(vec_N, dim(test_model$model)[1])
      vec_Rsquared <- c(vec_Rsquared, summary(test_model)$r.squared)
      if (str_individual_covariate == "Sex"){
        vec_coefficient <- c(vec_coefficient, summary(test_model)$coefficients[, "Estimate"]["SexMale"])
        vec_tvalue <- c(vec_tvalue, summary(test_model)$coefficients[, "t value"]["SexMale"])
        vec_prtvalue <- c(vec_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"]["SexMale"])
      } else{
        vec_coefficient <- c(vec_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_covariate])
        vec_tvalue <- c(vec_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_covariate])
        vec_prtvalue <- c(vec_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_covariate])
      }
    }
    
    for (str_individual_covariate in vec_allancestries_individual_covariates){
      test_model <- lm(formula(paste0("log10(Proportion) ~ ", str_individual_covariate)), data = df_merged_without_European)
      
      sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualCovariate_AllAncestries", str_celltype, str_individual_covariate, ".txt", sep = "_"))
      print(str_celltype)
      print(str_individual_covariate)
      print(summary(test_model))
      sink()
      vec_allancestries_model <- c(vec_allancestries_model, as.character(summary(test_model)$call)[2])
      vec_allancestries_covariate <- c(vec_allancestries_covariate, str_individual_covariate)
      vec_allancestries_subtype <- c(vec_allancestries_subtype, str_celltype)
      vec_allancestries_df <- c(vec_allancestries_df, test_model$df.residual)
      vec_allancestries_N <- c(vec_allancestries_N, dim(test_model$model)[1])
      vec_allancestries_Rsquared <- c(vec_allancestries_Rsquared, summary(test_model)$r.squared)
      if (str_individual_covariate == "Sex"){
        vec_allancestries_coefficient <- c(vec_allancestries_coefficient, summary(test_model)$coefficients[, "Estimate"]["SexMale"])
        vec_allancestries_tvalue <- c(vec_allancestries_tvalue, summary(test_model)$coefficients[, "t value"]["SexMale"])
        vec_allancestries_prtvalue <- c(vec_allancestries_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"]["SexMale"])
      } else{
        vec_allancestries_coefficient <- c(vec_allancestries_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_covariate])
        vec_allancestries_tvalue <- c(vec_allancestries_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_covariate])
        vec_allancestries_prtvalue <- c(vec_allancestries_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_covariate])
      }
    }
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Japanese + SG_Chinese + SG_Malay + SG_Indian + Thai, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualAncestryvsKorean", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Korean + Thai + SG_Chinese + SG_Malay + SG_Indian, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualAncestryvsJapanese", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    for (str_individual_ancestry in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
      test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", str_individual_ancestry)), data = df_merged_without_European)
      
      sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualAncestry", str_celltype, str_individual_ancestry, ".txt", sep = "_"))
      print(str_celltype)
      print(str_individual_ancestry)
      print(summary(test_model))
      sink()
      vec_ancestry_model <- c(vec_ancestry_model, as.character(summary(test_model)$call)[2])
      vec_ancestry_ancestry <- c(vec_ancestry_ancestry, str_individual_ancestry)
      vec_ancestry_subtype <- c(vec_ancestry_subtype, str_celltype)
      vec_ancestry_df <- c(vec_ancestry_df, test_model$df.residual)
      vec_ancestry_N <- c(vec_ancestry_N, dim(test_model$model)[1])
      vec_ancestry_Rsquared <- c(vec_ancestry_Rsquared, summary(test_model)$r.squared)
      vec_ancestry_coefficient <- c(vec_ancestry_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_ancestry])
      vec_ancestry_tvalue <- c(vec_ancestry_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_ancestry])
      vec_ancestry_prtvalue <- c(vec_ancestry_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_ancestry])
    }
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry + Sex:Ancestry + Age:Ancestry, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, "InteractionsAncestry.txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry + Sex:Age + Sex:Ancestry + Age:Ancestry, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, "AllTwoWayInteractions.txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
  }
}  

df_individual_covariates <- data.frame(Model = vec_model, 
                                       Covariate = vec_covariate, 
                                       Subtype = vec_subtype, 
                                       df = vec_df, 
                                       Number = vec_N, 
                                       Rsquared = vec_Rsquared, 
                                       Coefficient = vec_coefficient, 
                                       tvalue = vec_tvalue, 
                                       pvalue = vec_prtvalue)
write.table(df_individual_covariates, file = "df_individual_covariates.txt", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

df_individual_ancestries <- data.frame(Model = vec_ancestry_model, 
                                       Covariate = vec_ancestry_ancestry, 
                                       Subtype = vec_ancestry_subtype, 
                                       df = vec_ancestry_df, 
                                       Number = vec_ancestry_N, 
                                       Rsquared = vec_ancestry_Rsquared, 
                                       Coefficient = vec_ancestry_coefficient, 
                                       tvalue = vec_ancestry_tvalue, 
                                       pvalue = vec_ancestry_prtvalue)
write.table(df_individual_ancestries, file = "df_individual_ancestries.txt", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

df_all_ancestries <- data.frame(Model = vec_allancestries_model, 
                                Covariate = vec_allancestries_covariate, 
                                Subtype = vec_allancestries_subtype, 
                                df = vec_allancestries_df, 
                                Number = vec_allancestries_N, 
                                Rsquared = vec_allancestries_Rsquared, 
                                Coefficient = vec_allancestries_coefficient, 
                                tvalue = vec_allancestries_tvalue, 
                                pvalue = vec_allancestries_prtvalue)
write.table(df_all_ancestries, file = "df_all_ancestries.txt", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

df_individual_covariates$Covariate <- factor(df_individual_covariates$Covariate, 
                                             levels = vec_individual_covariates)
ggplot(df_individual_covariates, aes(x = Covariate, y = Rsquared)) +
  geom_boxplot() +
  #  scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
  ggtitle(paste("Variance explained by individual covariates", collapse = " ")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 15)) +
  theme(legend.title=element_text(size = 12), 
        legend.text=element_text(size = 12))

ggplot(df_individual_covariates[which(df_individual_covariates$Subtype == "CD16+_NK"), ], 
       aes(x = Covariate, y = Rsquared)) + 
  geom_boxplot(aes(fill = Covariate), alpha = .2) +
  geom_line(aes(group = Subtype)) + 
  geom_point(size = 2)

ggplot(df_all_ancestries, aes(x = Covariate, y = Rsquared)) +
  geom_boxplot() +
  #  scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
  ggtitle(paste("Variance explained by individual covariates", collapse = " ")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 15)) #+
# theme(legend.title=element_text(size = 12), 
#       legend.text=element_text(size = 12))

ggplot(df_all_ancestries, aes(x = Covariate, y = sqrt(Rsquared))) +
  geom_boxplot() +
  #  scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
  ggtitle(paste("Variance explained by individual covariates", collapse = " ")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 15)) #+
# theme(legend.title=element_text(size = 12), 
#       legend.text=element_text(size = 12))

ggplot(df_all_ancestries, 
       aes(x = Covariate, y = log10(Rsquared))) + 
  geom_boxplot(aes(fill = Covariate), alpha = .2) +
  geom_line(aes(group = Subtype)) + 
  geom_point(size = 2)

##### Interaction terms for individual ancestries #####

### For Annotation_Level3, naive_B is i = 14, CD16+_NK is i = 1, CD4+_T_naive is i = 2

i <- 14
i <- 1
i <- 2

### Ran the Level3 code

for (str_individual_ancestry in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
  test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", 
                                  str_individual_ancestry, 
                                  " + Sex:Age + Sex:", str_individual_ancestry, 
                                  " + Age:", str_individual_ancestry)), data = df_merged_without_European)
  
  sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualAncestry", str_celltype, str_individual_ancestry, "AllTwoWayInteractions.txt", sep = "_"))
  print(str_celltype)
  print(str_individual_ancestry)
  print(summary(test_model))
  sink()
}

##### Option of considering only CD4+ T cells as denominator for Treg #####

vec_TNK <- c(
  which(df_AIDA_metadata$Annotation_Level2 == "CD4+_T")
)

df_total_cell_count_filtered <- df_AIDA_metadata[vec_TNK, c("DCP_ID", "count")]
df_total_cell_count_filtered <- aggregate(count ~ ., df_total_cell_count_filtered, FUN = sum)
rownames(df_total_cell_count_filtered) <- df_total_cell_count_filtered$DCP_ID

vec_unique_celltypes <- unique(df_AIDA_metadata[vec_TNK, ]$Annotation_Level3)

### Figure Panel for cell type proportions

library(ggplot2)
vec_individual_covariates <- c("Sex", "Age", "SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")
vec_model <- c()
vec_covariate <- c()
vec_subtype <- c()
vec_df <- c()
vec_N <- c()
vec_Rsquared <- c()
vec_coefficient <- c()
vec_tvalue <- c()
vec_prtvalue <- c()

vec_ancestry_model <- c()
vec_ancestry_ancestry <- c()
vec_ancestry_subtype <- c()
vec_ancestry_df <- c()
vec_ancestry_N <- c()
vec_ancestry_Rsquared <- c()
vec_ancestry_coefficient <- c()
vec_ancestry_tvalue <- c()
vec_ancestry_prtvalue <- c()

# Set i <- 5 for Treg

for (i in 1:length(vec_unique_celltypes)) {
  str_celltype <- vec_unique_celltypes[i]
  vec_rows_celltype <- which(df_AIDA_metadata$Annotation_Level3 == str_celltype)
  df_celltype <- df_AIDA_metadata[vec_rows_celltype, c("DCP_ID", "count")]
  df_celltype <- aggregate(count ~ ., df_celltype, FUN = sum)
  colnames(df_celltype) <- c("DCP_ID", vec_unique_celltypes[i])
  rownames(df_celltype) <- df_celltype$DCP_ID
  df_celltype$count <- df_total_cell_count_filtered[rownames(df_celltype), "count"]
  
  df_celltype_metadata <- df_AIDA_metadata[match(df_celltype$DCP_ID, df_AIDA_metadata$DCP_ID),
                                           c("Genotyping_ID", "Library", "Country",
                                             "DCP_ID", "Age", "Sex", "Ancestry")]
  rownames(df_celltype_metadata) <- df_celltype_metadata$DCP_ID
  
  df_merged <- merge(df_celltype_metadata, df_celltype, by = "DCP_ID")
  df_merged$Age <- as.numeric(df_merged$Age)
  df_merged$Proportion <- df_merged[, vec_unique_celltypes[i]]/df_merged$count
  
  df_merged_without_European <- df_merged[-which(df_merged$Ancestry == "European"),]
  df_merged_without_European$Ancestry <- factor(df_merged_without_European$Ancestry, 
                                                levels = c("SG_Chinese", "SG_Indian", "SG_Malay", 
                                                           "Japanese", "Korean", "Thai"))
  
  for (str_ancestry in unique(df_merged_without_European$Ancestry)){
    df_merged_without_European[ , str_ancestry] <- 0
    df_merged_without_European[ , str_ancestry][which(df_merged_without_European$Ancestry == str_ancestry)] <- 1
  }
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Ancestry, colour = Ancestry)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin each ancestry out of all CD4+ T cells per donor", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion_Treg_CD4+_T_denominator_ancestry.pdf", sep = "_"), 
         height = 4, width = 6, dpi = 300, units = "in")
  
  for (str_individual_ancestry in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
    test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", 
                                    str_individual_ancestry, 
                                    " + Sex:Age + Sex:", str_individual_ancestry, 
                                    " + Age:", str_individual_ancestry)), data = df_merged_without_European)
    
    sink(file = paste("Treg_CD4+_T_denominator_LinearModel_Figure_Log10Proportion_IndividualAncestry", str_individual_ancestry, "AllTwoWayInteractions.txt", sep = "_"))
    print(str_individual_ancestry)
    print(summary(test_model))
    sink()
  }
  
  for (str_individual_ancestry in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
    test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", 
                                    str_individual_ancestry)), data = df_merged_without_European)
    
    sink(file = paste("Treg_CD4+_T_denominator_LinearModel_Figure_Log10Proportion_IndividualAncestry", str_individual_ancestry, ".txt", sep = "_"))
    print(str_individual_ancestry)
    print(summary(test_model))
    sink()
  }
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "ancestry.pdf", sep = "_"))
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Sex)) + 
    geom_boxplot() +
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin each sex", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex.pdf", sep = "_"))
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Ancestry, fill = Sex)) + 
    geom_boxplot() + #(aes(colour = Ancestry)) +
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby sex and ancestry", collapse = " ")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex_ancestry.pdf", sep = "_"))
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age)) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "by age", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "age.pdf", sep = "_"))
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age, colour = Ancestry)) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby age, trend by ancestry", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "age_by_Ancestry.pdf", sep = "_"))
  
  pdf(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex_ancestry_interaction.pdf", sep = "_"), 
      width = 700, height = 700)
  par(mar=c(5.1, 6.1, 4.1, 2.1))
  interaction.plot(x.factor = df_merged_without_European$Sex,
                   trace.factor = df_merged_without_European$Ancestry,
                   response = log10(df_merged_without_European$Proportion),
                   fun = median,
                   col = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3"),#c("red","green", "cyan", "blue", "pink"),
                   xlab = "Sex",
                   ylab = paste("Median log10(proportion) \nof", str_celltype), 
                   trace.label = "Ancestry",
                   lwd = 4,
                   lty = 5,
                   cex.lab = 1.2, 
                   cex.axis = 1.2)
  title(main = paste("Interaction plot of Sex and Ancestry on median log10(proportion) of", 
                     str_celltype))
  dev.off()
  
  {
    for (str_individual_covariate in vec_individual_covariates){
      test_model <- lm(formula(paste0("log10(Proportion) ~ ", str_individual_covariate)), data = df_merged_without_European)
      
      sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualCovariate", str_celltype, str_individual_covariate, ".txt", sep = "_"))
      print(str_celltype)
      print(str_individual_covariate)
      print(summary(test_model))
      sink()
      vec_model <- c(vec_model, as.character(summary(test_model)$call)[2])
      vec_covariate <- c(vec_covariate, str_individual_covariate)
      vec_subtype <- c(vec_subtype, str_celltype)
      vec_df <- c(vec_df, test_model$df.residual)
      vec_N <- c(vec_N, dim(test_model$model)[1])
      vec_Rsquared <- c(vec_Rsquared, summary(test_model)$r.squared)
      vec_coefficient <- c(vec_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_covariate])
      vec_tvalue <- c(vec_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_covariate])
      vec_prtvalue <- c(vec_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_covariate])
    }
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry, data = df_merged_without_European)
    
    sink(file = paste("Treg_CD4+_T_denominator_LinearModel_Figure_Log10Proportion", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Japanese + SG_Chinese + SG_Malay + SG_Indian + Thai, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualAncestryvsKorean", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Korean + Thai + SG_Chinese + SG_Malay + SG_Indian, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualAncestryvsJapanese", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    for (str_individual_ancestry in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
      test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", str_individual_ancestry)), data = df_merged_without_European)
      
      sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualAncestry", str_celltype, str_individual_ancestry, ".txt", sep = "_"))
      print(str_celltype)
      print(str_individual_ancestry)
      print(summary(test_model))
      sink()
      vec_ancestry_model <- c(vec_ancestry_model, as.character(summary(test_model)$call)[2])
      vec_ancestry_ancestry <- c(vec_ancestry_ancestry, str_individual_ancestry)
      vec_ancestry_subtype <- c(vec_ancestry_subtype, str_celltype)
      vec_ancestry_df <- c(vec_ancestry_df, test_model$df.residual)
      vec_ancestry_N <- c(vec_ancestry_N, dim(test_model$model)[1])
      vec_ancestry_Rsquared <- c(vec_ancestry_Rsquared, summary(test_model)$r.squared)
      vec_ancestry_coefficient <- c(vec_ancestry_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_ancestry])
      vec_ancestry_tvalue <- c(vec_ancestry_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_ancestry])
      vec_ancestry_prtvalue <- c(vec_ancestry_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_ancestry])
    }
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry + Sex:Ancestry + Age:Ancestry, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, "InteractionsAncestry.txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Ancestry + Sex:Age + Sex:Ancestry + Age:Ancestry, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, "AllTwoWayInteractions.txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
  }
}