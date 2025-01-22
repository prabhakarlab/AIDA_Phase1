### AIDA Phase 1 Data Freeze v2 Cell Type Proportions analysis: version 22 January 2025

library(ggplot2)
library(dplyr)

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
                                             "DCP_ID", "Age", "Sex", "ethnicity")]
  rownames(df_celltype_metadata) <- df_celltype_metadata$DCP_ID
  
  df_merged <- merge(df_celltype_metadata, df_celltype, by = "DCP_ID")
  df_merged$Age <- as.numeric(df_merged$Age)
  df_merged$Proportion <- df_merged[, vec_unique_celltypes[i]]/df_merged$count
  
  df_merged_without_European <- df_merged[-which(df_merged$ethnicity == "European"),]
  df_merged_without_European$ethnicity <- factor(df_merged_without_European$ethnicity, 
                                                levels = c("SG_Chinese", "SG_Indian", "SG_Malay", 
                                                           "Japanese", "Korean", "Thai"))
  
  for (str_ethnicity in unique(df_merged_without_European$ethnicity)){
    df_merged_without_European[ , str_ethnicity] <- 0
    df_merged_without_European[ , str_ethnicity][which(df_merged_without_European$ethnicity == str_ethnicity)] <- 1
  }
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = ethnicity, colour = ethnicity)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "in each ethnicity", collapse = " ")) + 
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 10)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 10))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "ethnicity.pdf", sep = "_"), 
         height = 4, width = 5, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Sex)) + 
    geom_boxplot() +
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin females versus males", collapse = " ")) + 
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "sex.pdf", sep = "_"), 
         height = 4, width = 3, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = ethnicity, fill = Sex)) + 
    geom_boxplot() +
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby sex and ethnicity", collapse = " ")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 10))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "sex_ethnicity.pdf", sep = "_"))
  
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
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age, colour = ethnicity)) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby age, trend by ethnicity", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 10))
  
  ggsave(filename = paste("Level1/RPlot_Figure_Log10Proportion", str_celltype, "age_by_ethnicity.pdf", sep = "_"))
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion", str_celltype, ".txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + Japanese + SG_Chinese + SG_Malay + SG_Indian + Thai, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion_IndividualethnicityvsKorean", str_celltype, ".txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + Korean + Thai + SG_Chinese + SG_Malay + SG_Indian, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion_IndividualethnicityvsJapanese", str_celltype, ".txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
  
  for (str_individual_ethnicity in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
    test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", str_individual_ethnicity)), data = df_merged_without_European)
    
    sink(file = paste("Level1/LinearModel_Figure_Log10Proportion_Individualethnicity_SG_Chinese", str_celltype, str_individual_ethnicity, ".txt", sep = "_"))
    print(str_celltype)
    print(str_individual_ethnicity)
    print(summary(test_model))
    sink()
  }
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity + Sex:ethnicity + Age:ethnicity, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion", str_celltype, "Interactionsethnicity.txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
  
  test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity + Sex:Age + Sex:ethnicity + Age:ethnicity, data = df_merged_without_European)
  
  sink(file = paste("Level1/LinearModel_Figure_Log10Proportion", str_celltype, "AllTwoWayInteractions.txt", sep = "_"))
  print(str_celltype)
  print(summary(test_model))
  sink()
}  

##### Annotation Level 3 against all PBMCs #####
{
  vec_unique_celltypes <- unique(df_AIDA_metadata$Annotation_Level3)
  
  ### Figure Panel for cell type proportions
  
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
  
  vec_allethnicities_individual_covariates <- c("Sex", "Age", "ethnicity")
  vec_allethnicities_model <- c()
  vec_allethnicities_covariate <- c()
  vec_allethnicities_subtype <- c()
  vec_allethnicities_df <- c()
  vec_allethnicities_N <- c()
  vec_allethnicities_Rsquared <- c()
  vec_allethnicities_coefficient <- c()
  vec_allethnicities_tvalue <- c()
  vec_allethnicities_prtvalue <- c()
  
  vec_ethnicity_model <- c()
  vec_ethnicity_ethnicity <- c()
  vec_ethnicity_subtype <- c()
  vec_ethnicity_df <- c()
  vec_ethnicity_N <- c()
  vec_ethnicity_Rsquared <- c()
  vec_ethnicity_coefficient <- c()
  vec_ethnicity_tvalue <- c()
  vec_ethnicity_prtvalue <- c()
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
                                             "DCP_ID", "Age", "Sex", "ethnicity")]
  rownames(df_celltype_metadata) <- df_celltype_metadata$DCP_ID
  
  df_merged <- merge(df_celltype_metadata, df_celltype, by = "DCP_ID")
  df_merged$Age <- as.numeric(df_merged$Age)
  df_merged$Proportion <- df_merged[, vec_unique_celltypes[i]]/df_merged$count
  
  df_merged_without_European <- df_merged[-which(df_merged$ethnicity == "European"),]
  df_merged_without_European$ethnicity <- factor(df_merged_without_European$ethnicity, 
                                                levels = c("SG_Chinese", "SG_Indian", "SG_Malay", 
                                                           "Japanese", "Korean", "Thai"))
  
  for (str_ethnicity in unique(df_merged_without_European$ethnicity)){
    df_merged_without_European[ , str_ethnicity] <- 0
    df_merged_without_European[ , str_ethnicity][which(df_merged_without_European$ethnicity == str_ethnicity)] <- 1
  }
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = ethnicity, colour = ethnicity)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin each ethnicity", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "ethnicity.pdf", sep = "_"), 
         height = 4, width = 6, dpi = 300, units = "in")
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Sex)) +
    geom_boxplot() +
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin each sex", collapse = " ")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15))

  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex.pdf", sep = "_"))

  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = ethnicity, fill = Sex)) + 
    geom_boxplot() +
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby sex and ethnicity", collapse = " ")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex_ethnicity.pdf", sep = "_"), 
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
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age, colour = ethnicity)) +
    geom_point() +
    geom_smooth(method = "lm") + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby age, trend by ethnicity", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "age_by_ethnicity.pdf", sep = "_"), 
         height = 4, width = 5, dpi = 300, units = "in")
  
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
    
    for (str_individual_covariate in vec_allethnicities_individual_covariates){
      test_model <- lm(formula(paste0("log10(Proportion) ~ ", str_individual_covariate)), data = df_merged_without_European)
      
      sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualCovariate_Allethnicities", str_celltype, str_individual_covariate, ".txt", sep = "_"))
      print(str_celltype)
      print(str_individual_covariate)
      print(summary(test_model))
      sink()
      vec_allethnicities_model <- c(vec_allethnicities_model, as.character(summary(test_model)$call)[2])
      vec_allethnicities_covariate <- c(vec_allethnicities_covariate, str_individual_covariate)
      vec_allethnicities_subtype <- c(vec_allethnicities_subtype, str_celltype)
      vec_allethnicities_df <- c(vec_allethnicities_df, test_model$df.residual)
      vec_allethnicities_N <- c(vec_allethnicities_N, dim(test_model$model)[1])
      vec_allethnicities_Rsquared <- c(vec_allethnicities_Rsquared, summary(test_model)$r.squared)
      if (str_individual_covariate == "Sex"){
        vec_allethnicities_coefficient <- c(vec_allethnicities_coefficient, summary(test_model)$coefficients[, "Estimate"]["SexMale"])
        vec_allethnicities_tvalue <- c(vec_allethnicities_tvalue, summary(test_model)$coefficients[, "t value"]["SexMale"])
        vec_allethnicities_prtvalue <- c(vec_allethnicities_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"]["SexMale"])
      } else{
        vec_allethnicities_coefficient <- c(vec_allethnicities_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_covariate])
        vec_allethnicities_tvalue <- c(vec_allethnicities_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_covariate])
        vec_allethnicities_prtvalue <- c(vec_allethnicities_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_covariate])
      }
    }
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Japanese + SG_Chinese + SG_Malay + SG_Indian + Thai, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualethnicityvsKorean", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Korean + Thai + SG_Chinese + SG_Malay + SG_Indian, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualethnicityvsJapanese", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    for (str_individual_ethnicity in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
      test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", str_individual_ethnicity)), data = df_merged_without_European)
      
      sink(file = paste("LinearModel_Figure_Log10Proportion_Individualethnicity", str_celltype, str_individual_ethnicity, ".txt", sep = "_"))
      print(str_celltype)
      print(str_individual_ethnicity)
      print(summary(test_model))
      sink()
      vec_ethnicity_model <- c(vec_ethnicity_model, as.character(summary(test_model)$call)[2])
      vec_ethnicity_ethnicity <- c(vec_ethnicity_ethnicity, str_individual_ethnicity)
      vec_ethnicity_subtype <- c(vec_ethnicity_subtype, str_celltype)
      vec_ethnicity_df <- c(vec_ethnicity_df, test_model$df.residual)
      vec_ethnicity_N <- c(vec_ethnicity_N, dim(test_model$model)[1])
      vec_ethnicity_Rsquared <- c(vec_ethnicity_Rsquared, summary(test_model)$r.squared)
      vec_ethnicity_coefficient <- c(vec_ethnicity_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_ethnicity])
      vec_ethnicity_tvalue <- c(vec_ethnicity_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_ethnicity])
      vec_ethnicity_prtvalue <- c(vec_ethnicity_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_ethnicity])
    }
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity + Sex:ethnicity + Age:ethnicity, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, "Interactionsethnicity.txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity + Sex:Age + Sex:ethnicity + Age:ethnicity, data = df_merged_without_European)
    
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

df_individual_ethnicities <- data.frame(Model = vec_ethnicity_model, 
                                       Covariate = vec_ethnicity_ethnicity, 
                                       Subtype = vec_ethnicity_subtype, 
                                       df = vec_ethnicity_df, 
                                       Number = vec_ethnicity_N, 
                                       Rsquared = vec_ethnicity_Rsquared, 
                                       Coefficient = vec_ethnicity_coefficient, 
                                       tvalue = vec_ethnicity_tvalue, 
                                       pvalue = vec_ethnicity_prtvalue)
write.table(df_individual_ethnicities, file = "df_individual_ethnicities.txt", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

df_all_ethnicities <- data.frame(Model = vec_allethnicities_model, 
                                Covariate = vec_allethnicities_covariate, 
                                Subtype = vec_allethnicities_subtype, 
                                df = vec_allethnicities_df, 
                                Number = vec_allethnicities_N, 
                                Rsquared = vec_allethnicities_Rsquared, 
                                Coefficient = vec_allethnicities_coefficient, 
                                tvalue = vec_allethnicities_tvalue, 
                                pvalue = vec_allethnicities_prtvalue)
write.table(df_all_ethnicities, file = "df_all_ethnicities.txt", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

df_individual_covariates$Covariate <- factor(df_individual_covariates$Covariate, 
                                             levels = vec_individual_covariates)
ggplot(df_individual_covariates, aes(x = Covariate, y = Rsquared)) +
  geom_boxplot() +
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

ggplot(df_all_ethnicities, aes(x = Covariate, y = Rsquared)) +
  geom_boxplot() +
  ggtitle(paste("Variance explained by individual covariates", collapse = " ")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 15))

ggplot(df_all_ethnicities, aes(x = Covariate, y = sqrt(Rsquared))) +
  geom_boxplot() +
  ggtitle(paste("Variance explained by individual covariates", collapse = " ")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 15))

ggplot(df_all_ethnicities, 
       aes(x = Covariate, y = log10(Rsquared))) + 
  geom_boxplot(aes(fill = Covariate), alpha = .2) +
  geom_line(aes(group = Subtype)) + 
  geom_point(size = 2)

##### Interaction terms for individual ethnicities #####

### For Annotation_Level3, naive_B is i = 14, CD16+_NK is i = 1, CD4+_T_naive is i = 2

i <- 14
i <- 1
i <- 2

### Ran the Level3 code

for (str_individual_ethnicity in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
  test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", 
                                  str_individual_ethnicity, 
                                  " + Sex:Age + Sex:", str_individual_ethnicity, 
                                  " + Age:", str_individual_ethnicity)), data = df_merged_without_European)
  
  sink(file = paste("LinearModel_Figure_Log10Proportion_Individualethnicity", str_celltype, str_individual_ethnicity, "AllTwoWayInteractions.txt", sep = "_"))
  print(str_celltype)
  print(str_individual_ethnicity)
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

vec_ethnicity_model <- c()
vec_ethnicity_ethnicity <- c()
vec_ethnicity_subtype <- c()
vec_ethnicity_df <- c()
vec_ethnicity_N <- c()
vec_ethnicity_Rsquared <- c()
vec_ethnicity_coefficient <- c()
vec_ethnicity_tvalue <- c()
vec_ethnicity_prtvalue <- c()

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
                                             "DCP_ID", "Age", "Sex", "ethnicity")]
  rownames(df_celltype_metadata) <- df_celltype_metadata$DCP_ID
  
  df_merged <- merge(df_celltype_metadata, df_celltype, by = "DCP_ID")
  df_merged$Age <- as.numeric(df_merged$Age)
  df_merged$Proportion <- df_merged[, vec_unique_celltypes[i]]/df_merged$count
  
  df_merged_without_European <- df_merged[-which(df_merged$ethnicity == "European"),]
  df_merged_without_European$ethnicity <- factor(df_merged_without_European$ethnicity, 
                                                levels = c("SG_Chinese", "SG_Indian", "SG_Malay", 
                                                           "Japanese", "Korean", "Thai"))
  
  for (str_ethnicity in unique(df_merged_without_European$ethnicity)){
    df_merged_without_European[ , str_ethnicity] <- 0
    df_merged_without_European[ , str_ethnicity][which(df_merged_without_European$ethnicity == str_ethnicity)] <- 1
  }
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = ethnicity, colour = ethnicity)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin each ethnicity out of all CD4+ T cells per donor", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion_Treg_CD4+_T_denominator_ethnicity.pdf", sep = "_"), 
         height = 4, width = 6, dpi = 300, units = "in")
  
  for (str_individual_ethnicity in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
    test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", 
                                    str_individual_ethnicity, 
                                    " + Sex:Age + Sex:", str_individual_ethnicity, 
                                    " + Age:", str_individual_ethnicity)), data = df_merged_without_European)
    
    sink(file = paste("Treg_CD4+_T_denominator_LinearModel_Figure_Log10Proportion_Individualethnicity", str_individual_ethnicity, "AllTwoWayInteractions.txt", sep = "_"))
    print(str_individual_ethnicity)
    print(summary(test_model))
    sink()
  }
  
  for (str_individual_ethnicity in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
    test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", 
                                    str_individual_ethnicity)), data = df_merged_without_European)
    
    sink(file = paste("Treg_CD4+_T_denominator_LinearModel_Figure_Log10Proportion_Individualethnicity", str_individual_ethnicity, ".txt", sep = "_"))
    print(str_individual_ethnicity)
    print(summary(test_model))
    sink()
  }
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "ethnicity.pdf", sep = "_"))
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Sex)) + 
    geom_boxplot() +
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nin each sex", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex.pdf", sep = "_"))
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = ethnicity, fill = Sex)) + 
    geom_boxplot() + #(aes(colour = ethnicity)) +
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby sex and ethnicity", collapse = " ")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "sex_ethnicity.pdf", sep = "_"))
  
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
  
  ggplot(data = df_merged_without_European, aes(y = log10(Proportion), x = Age, colour = ethnicity)) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    scale_colour_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#619CFF", "#00BA38", "#00BFC4", "#F564E3")) + 
    ggtitle(paste("Log10(Proportion) of", str_celltype, "\nby age, trend by ethnicity", collapse = " ")) + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme(legend.title=element_text(size = 12), 
          legend.text=element_text(size = 12))
  
  ggsave(filename = paste("RPlot_Figure_Log10Proportion", str_celltype, "age_by_ethnicity.pdf", sep = "_"))
  
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
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity, data = df_merged_without_European)
    
    sink(file = paste("Treg_CD4+_T_denominator_LinearModel_Figure_Log10Proportion", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Japanese + SG_Chinese + SG_Malay + SG_Indian + Thai, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualethnicityvsKorean", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + Korean + Thai + SG_Chinese + SG_Malay + SG_Indian, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion_IndividualethnicityvsJapanese", str_celltype, ".txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    for (str_individual_ethnicity in c("SG_Chinese", "SG_Indian", "Japanese", "Korean", "SG_Malay", "Thai")){
      test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", str_individual_ethnicity)), data = df_merged_without_European)
      
      sink(file = paste("LinearModel_Figure_Log10Proportion_Individualethnicity", str_celltype, str_individual_ethnicity, ".txt", sep = "_"))
      print(str_celltype)
      print(str_individual_ethnicity)
      print(summary(test_model))
      sink()
      vec_ethnicity_model <- c(vec_ethnicity_model, as.character(summary(test_model)$call)[2])
      vec_ethnicity_ethnicity <- c(vec_ethnicity_ethnicity, str_individual_ethnicity)
      vec_ethnicity_subtype <- c(vec_ethnicity_subtype, str_celltype)
      vec_ethnicity_df <- c(vec_ethnicity_df, test_model$df.residual)
      vec_ethnicity_N <- c(vec_ethnicity_N, dim(test_model$model)[1])
      vec_ethnicity_Rsquared <- c(vec_ethnicity_Rsquared, summary(test_model)$r.squared)
      vec_ethnicity_coefficient <- c(vec_ethnicity_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_ethnicity])
      vec_ethnicity_tvalue <- c(vec_ethnicity_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_ethnicity])
      vec_ethnicity_prtvalue <- c(vec_ethnicity_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_ethnicity])
    }
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity + Sex:ethnicity + Age:ethnicity, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, "Interactionsethnicity.txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
    
    test_model <- lm(log10(Proportion) ~ Age + Sex + ethnicity + Sex:Age + Sex:ethnicity + Age:ethnicity, data = df_merged_without_European)
    
    sink(file = paste("LinearModel_Figure_Log10Proportion", str_celltype, "AllTwoWayInteractions.txt", sep = "_"))
    print(str_celltype)
    print(summary(test_model))
    sink()
  }
}

##### Analysis of AIDA Singapore donors: incremental variance explained, SLAS-2 comparisons #####

##### INPUT REQUIRED: READ IN GENOTYPE PCS AS df_PCA AND SMOKING METADATA AS df_smoking #####

##### Consider only donors with at least 800 cells and Singapore ethnicities #####

df_total_cell_count <- df_AIDA_metadata[, c("DCP_ID", "count")]
df_total_cell_count <- aggregate(count ~ ., df_total_cell_count, FUN = sum)
rownames(df_total_cell_count) <- df_total_cell_count$DCP_ID
vec_donors_to_exclude <- df_total_cell_count$DCP_ID[which(df_total_cell_count$count < 800)]
df_AIDA_metadata <- df_AIDA_metadata[!(df_AIDA_metadata$DCP_ID %in% vec_donors_to_exclude), ]
df_AIDA_metadata <- df_AIDA_metadata[(df_AIDA_metadata$Country == "SG"), ]

##### Consider only non-platelet cells

vec_platelet_to_exclude <- rownames(df_AIDA_metadata)[df_AIDA_metadata$Annotation_Level1 == "Platelet"]
df_AIDA_metadata <- df_AIDA_metadata[-match(vec_platelet_to_exclude, rownames(df_AIDA_metadata)), ]

df_AIDA_metadata$Smoking <- df_smoking[df_AIDA_metadata$DCP_ID, "Smoking_Status"]
df_AIDA_metadata$PC1 <- df_PCA[df_AIDA_metadata$DCP_ID, "PC1"]
df_AIDA_metadata$PC2 <- df_PCA[df_AIDA_metadata$DCP_ID, "PC2"]
df_AIDA_metadata$PC3 <- df_PCA[df_AIDA_metadata$DCP_ID, "PC3"]

##### Level 3 against all PBMCs #####

vec_unique_celltypes <- unique(df_AIDA_metadata$Annotation_Level3)

vec_ethnicity_model <- c()
vec_ethnicity_ethnicity <- c()
vec_ethnicity_subtype <- c()
vec_ethnicity_df <- c()
vec_ethnicity_N <- c()
vec_ethnicity_Rsquared <- c()
vec_ethnicity_coefficient <- c()
vec_ethnicity_tvalue <- c()
vec_ethnicity_prtvalue <- c()

vec_incremental_Rsquaredadditional <- c()
vec_incremental_subtype <- c()
vec_incremental_covariate <- c()
vec_incremental_popn <- c()

##### Additional variance explained after base model #####

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
                                             "DCP_ID", "Age", "Sex", "Ethnicity", "BMI", 
                                             "Smoking", "PC1", "PC2", "PC3")]
  rownames(df_celltype_metadata) <- df_celltype_metadata$DCP_ID
  
  df_merged <- merge(df_celltype_metadata, df_celltype, by = "DCP_ID")
  df_merged$Age <- as.numeric(df_merged$Age)
  df_merged$Proportion <- df_merged[, vec_unique_celltypes[i]]/df_merged$count
  df_merged$ethnicity <- df_merged$Ethnicity
  df_merged$BMI <- as.numeric(df_merged$BMI)
  df_merged$Smoking <- as.character(df_merged$Smoking)
  
  df_merged <- df_merged %>% mutate(
    Smoking = ifelse(Smoking == 0, "No", Smoking),
    Smoking = ifelse(Smoking == 1, "Yes", Smoking),
    Smoking = as.factor(Smoking)
  )
  
  df_merged$PC1 <- as.numeric(df_merged$PC1)
  df_merged$PC2 <- as.numeric(df_merged$PC2)
  df_merged$PC3 <- as.numeric(df_merged$PC3)
  
  df_merged_without_European <- df_merged[-which(df_merged$ethnicity == "European"),]
  
  for (str_ethnicity in unique(df_merged_without_European$ethnicity)){
    df_merged_without_European[ , str_ethnicity] <- 0
    df_merged_without_European[ , str_ethnicity][which(df_merged_without_European$ethnicity == str_ethnicity)] <- 1
  }
  
  for (str_individual_ethnicity in c("Chinese", "Indian", "Malay")){
    test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Sex + ", str_individual_ethnicity)), data = df_merged_without_European)
    
    vec_ethnicity_model <- c(vec_ethnicity_model, as.character(summary(test_model)$call)[2])
    vec_ethnicity_ethnicity <- c(vec_ethnicity_ethnicity, str_individual_ethnicity)
    vec_ethnicity_subtype <- c(vec_ethnicity_subtype, str_celltype)
    vec_ethnicity_df <- c(vec_ethnicity_df, test_model$df.residual)
    vec_ethnicity_N <- c(vec_ethnicity_N, dim(test_model$model)[1])
    vec_ethnicity_Rsquared <- c(vec_ethnicity_Rsquared, summary(test_model)$r.squared)
    vec_ethnicity_coefficient <- c(vec_ethnicity_coefficient, summary(test_model)$coefficients[, "Estimate"][str_individual_ethnicity])
    vec_ethnicity_tvalue <- c(vec_ethnicity_tvalue, summary(test_model)$coefficients[, "t value"][str_individual_ethnicity])
    vec_ethnicity_prtvalue <- c(vec_ethnicity_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][str_individual_ethnicity])
  }
  
  ##### Examining variance explained by self-reported ethnicity, or genotype PC1 + PC2 + PC3, in separate models
  ##### Analysis of AIDA Singapore donors with all metadata fields available
  
  df_merged_without_European <- df_merged_without_European[c(which(df_merged_without_European$Smoking == "Yes"), 
                                                             which(df_merged_without_European$Smoking == "No")),]
  
  for (str_variable in c("ethnicity", "PC1 + PC2 + PC3")){
    
    test_model_all <- lm(formula(paste0("log10(Proportion) ~ Sex + Age + BMI + Smoking + ", str_variable)), 
                         data = df_merged_without_European)
    test_model_no_sex <- lm(formula(paste0("log10(Proportion) ~ Age + BMI + Smoking + ", str_variable)), 
                            data = df_merged_without_European)
    test_model_no_age <- lm(formula(paste0("log10(Proportion) ~ Sex + BMI + Smoking + ", str_variable)), 
                            data = df_merged_without_European)
    test_model_no_bmi <- lm(formula(paste0("log10(Proportion) ~ Sex + Age + Smoking + ", str_variable)), 
                            data = df_merged_without_European)
    test_model_no_smoking <- lm(formula(paste0("log10(Proportion) ~ Sex + Age + BMI + ", str_variable)), 
                                data = df_merged_without_European)
    test_model_no_popn <- lm(formula(paste0("log10(Proportion) ~ Sex + Age + BMI + Smoking")), 
                             data = df_merged_without_European)
    
    vec_incremental_Rsquaredadditional <- c(vec_incremental_Rsquaredadditional, 
                                            summary(test_model_all)$r.squared-summary(test_model_no_sex)$r.squared, 
                                            summary(test_model_all)$r.squared-summary(test_model_no_age)$r.squared, 
                                            summary(test_model_all)$r.squared-summary(test_model_no_bmi)$r.squared, 
                                            summary(test_model_all)$r.squared-summary(test_model_no_smoking)$r.squared, 
                                            summary(test_model_all)$r.squared-summary(test_model_no_popn)$r.squared)
    vec_incremental_subtype <- c(vec_incremental_subtype, rep(str_celltype, 5))
    vec_incremental_covariate <- c(vec_incremental_covariate, "Sex", "Age", "BMI", "Smoking", str_variable)
    vec_incremental_popn <- c(vec_incremental_popn, rep(str_variable, 5))
    
  }
}  

df_Singapore_singlecell_individual_ethnicities <- data.frame(Model = vec_ethnicity_model,
                                                             Covariate = vec_ethnicity_ethnicity,
                                                             Subtype = vec_ethnicity_subtype,
                                                             df = vec_ethnicity_df,
                                                             N = vec_ethnicity_N,
                                                             Rsquared = vec_ethnicity_Rsquared,
                                                             Coefficient = vec_ethnicity_coefficient,
                                                             tvalue = vec_ethnicity_tvalue,
                                                             pvalue = vec_ethnicity_prtvalue)
df_incremental <- data.frame(Covariate = vec_incremental_covariate, 
                             Subtype = vec_incremental_subtype, 
                             Additional_Rsquared = vec_incremental_Rsquaredadditional , 
                             PopnDescriptor = vec_incremental_popn)

vec_rows_to_remove <- c(which(df_incremental$Subtype == "CD4+_T"),
                        which(df_incremental$Subtype == "T"),
                        which(df_incremental$Subtype == "CD8+_T"),
                        which(df_incremental$Subtype == "NK"),
                        which(df_incremental$Subtype == "Myeloid"),
                        which(df_incremental$Subtype == "memory_B"),
                        which(df_incremental$Subtype == "ILC"),
                        which(df_incremental$Subtype == "B"),
                        which(df_incremental$Subtype == "cDC1"),
                        which(df_incremental$Subtype == "DC_SIGLEC6hi"),
                        which(df_incremental$Subtype == "dnT"),
                        which(df_incremental$Subtype == "CD34_HSPC"))
df_incremental_corecelltypes <- df_incremental[-vec_rows_to_remove, ]

ggplot(df_incremental[-vec_rows_to_remove, ][which((df_incremental[-vec_rows_to_remove, ])$PopnDescriptor == "ethnicity"), ], aes(x = Covariate, y = Additional_Rsquared)) +
  geom_boxplot() +
  ggtitle(paste("Incremental variance explained by individual covariates", collapse = " ")) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15)) + 
  scale_y_log10()

ggplot(df_incremental[-vec_rows_to_remove, ][which((df_incremental[-vec_rows_to_remove, ])$PopnDescriptor == "PC1 + PC2 + PC3"), ], aes(x = Covariate, y = Additional_Rsquared)) +
  geom_boxplot() +
  ggtitle(paste("Incremental variance explained by individual covariates", collapse = " ")) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15)) + 
  scale_y_log10()

##### SLAS2 validation of AIDA cell type proportion covariate coefficients #####

##### INPUT REQUIRED: READ IN SLAS2 DATA AS df_SLAS2 AND SUPPLEMENTARY TABLE S4 AS df_subtypes #####

df_SLAS2$Female <- 1
df_SLAS2$Female[which(df_SLAS2$Sex == 1)] <- 0
df_SLAS2$Chinese <- 0
df_SLAS2$Malay <- 0
df_SLAS2$Indian <- 0
df_SLAS2$Chinese[which(df_SLAS2$Race == 1) ] <- 1
df_SLAS2$Malay[which(df_SLAS2$Race == 2)] <- 1
df_SLAS2$Indian[which(df_SLAS2$Race == 3)] <- 1
df_SLAS2$Race <- as.character(df_SLAS2$Race)
df_SLAS2$Sex <- as.character(df_SLAS2$Sex)
df_SLAS2$Female <- as.character(df_SLAS2$Female)
df_SLAS2$Chinese <- as.character(df_SLAS2$Chinese)
df_SLAS2$Malay <- as.character(df_SLAS2$Malay)
df_SLAS2$Indian <- as.character(df_SLAS2$Indian)
df_SLAS2$Ethnicity <- df_SLAS2$Race
df_SLAS2 <- df_SLAS2[-which(df_SLAS2$Race == "4"), ]

df_SLAS2 <- df_SLAS2 %>% mutate(
  Ethnicity = ifelse(Ethnicity == 1, "Chinese", Ethnicity),
  Ethnicity = ifelse(Ethnicity == 2, "Malay", Ethnicity),
  Ethnicity = ifelse(Ethnicity == 3, "Indian", Ethnicity),
  Ethnicity = as.factor(Ethnicity)
)

### Use Column X10 for PBMCs/Single Cells/Live Cells/CD34+CD45+ cells | Counts as denominator

rownames(df_subtypes) <- df_subtypes$Subtype

vec_ethnicity_model <- c()
vec_ethnicity_ethnicity <- c()
vec_ethnicity_subtype <- c()
vec_ethnicity_df <- c()
vec_ethnicity_N <- c()
vec_ethnicity_Rsquared <- c()
vec_ethnicity_coefficient <- c()
vec_ethnicity_tvalue <- c()
vec_ethnicity_prtvalue <- c()

for (str_celltype in rownames(df_subtypes)){
  str_fieldofinterest <- paste0("X", as.character(df_subtypes[str_celltype, "RowID"])) 
  df_SLAS2$Proportion <- df_SLAS2[, str_fieldofinterest]/df_SLAS2$X10
  df_SLAS2_subset <- df_SLAS2[-which(df_SLAS2$Proportion == 0), ]
  
  for (str_individual_ethnicity in c("Chinese", "Indian", "Malay")){
    test_model <- lm(formula(paste0("log10(Proportion) ~ Age + Female + ", str_individual_ethnicity)), data = df_SLAS2_subset)
    
    vec_ethnicity_model <- c(vec_ethnicity_model, as.character(summary(test_model)$call)[2])
    vec_ethnicity_ethnicity <- c(vec_ethnicity_ethnicity, str_individual_ethnicity)
    vec_ethnicity_subtype <- c(vec_ethnicity_subtype, str_celltype)
    vec_ethnicity_df <- c(vec_ethnicity_df, test_model$df.residual)
    vec_ethnicity_N <- c(vec_ethnicity_N, dim(test_model$model)[1])
    vec_ethnicity_Rsquared <- c(vec_ethnicity_Rsquared, summary(test_model)$r.squared)
    vec_ethnicity_coefficient <- c(vec_ethnicity_coefficient, summary(test_model)$coefficients[, "Estimate"][paste0(str_individual_ethnicity, "1")])
    vec_ethnicity_tvalue <- c(vec_ethnicity_tvalue, summary(test_model)$coefficients[, "t value"][paste0(str_individual_ethnicity, "1")])
    vec_ethnicity_prtvalue <- c(vec_ethnicity_prtvalue, summary(test_model)$coefficients[, "Pr(>|t|)"][paste0(str_individual_ethnicity, "1")])
  }
}

df_individual_ethnicities <- data.frame(Model = vec_ethnicity_model, 
                                        Covariate = vec_ethnicity_ethnicity, 
                                        Subtype = vec_ethnicity_subtype, 
                                        df = vec_ethnicity_df, 
                                        Number = vec_ethnicity_N, 
                                        Adj_Rsquared = vec_ethnicity_Rsquared, 
                                        Coefficient = vec_ethnicity_coefficient, 
                                        tvalue = vec_ethnicity_tvalue, 
                                        pvalue = vec_ethnicity_prtvalue)

df_SLAS2_to_check <- df_individual_ethnicities[, c("Covariate", "Subtype", "Coefficient", "pvalue")]
colnames(df_SLAS2_to_check)[which(colnames(df_SLAS2_to_check) == "Coefficient")] <- "SLAS2_Coefficient"
colnames(df_SLAS2_to_check)[which(colnames(df_SLAS2_to_check) == "pvalue")] <- "SLAS2_pvalue"

df_merged <- merge(df_Singapore_singlecell_individual_ethnicities, df_SLAS2_to_check, 
                   by = c("Covariate", "Subtype"))

ggplot(df_merged, aes(x = Coefficient, y = SLAS2_Coefficient)) +
  geom_point() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_smooth(method = "lm") + 
  xlab("AIDA Single-cell Coefficient Estimate") + ylab("SLAS2 Flow Cytometry \nCoefficient Estimate") + 
  ggtitle("Correlation of AIDA Single-cell and \nSLAS2 flow cytometry coefficients") + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15)) +
  theme(legend.title=element_text(size = 12), 
        legend.text=element_text(size = 12))

cor.test(df_merged$Coefficient, df_merged$SLAS2_Coefficient)