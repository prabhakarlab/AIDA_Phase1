library(RCAv2)

cell.df <- data.frame()
# SG
for (i in c(1,seq(3,19,1))){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("SG_HEL_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("SG_HEL_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    PBMC_r <- readRDS(paste0("../",lib_name,
                             "/RCA_PBMC_after_doublet_removal_and_before_QC.rds"))
    PBMC_r <- dataLogNormalise(PBMC_r)
    
    ############ project to all immune panels ################
    # get projection results against global panel
    PBMC_r <- dataProject(PBMC_r, method = "GlobalPanel",
                          corMeth = "pearson", scale = TRUE)
    global.proj <- PBMC_r$projection.data
    global.proj.immune <- read.table("../rownames_of_glocal_projection_immune_cells.txt", 
                                     header = FALSE)
    global.proj <- as.data.frame(as.matrix(global.proj[global.proj.immune$V1,]))
    
    # get projection results against other two panels
    PBMC_r <- dataProjectMultiPanel(PBMC_r,method = list("NovershternPanel", 
                                                         "MonacoPanel",
                                                         "CITESeqPanel"),
                                    scale = TRUE,corMeth = "pearson")
    two.proj <- as.data.frame(as.matrix(PBMC_r$projection.data))
    
    # combine these panels
    proj.all <- rbind(global.proj,two.proj)
    proj.all <- as.matrix(proj.all)
    proj.all <- as(proj.all, "dgCMatrix")
    
    # Assign projection result to RCA object
    PBMC_r$projection.data <- proj.all
    
    #Estimate the most probable cell type label for each cell
    PBMC_r <- estimateCellTypeFromProjection(PBMC_r,confidence = NULL)
    
    saveRDS(PBMC_r,paste0("./",lib_name,".singlets.RDS"))
    
    new_add <- data.frame(lib = lib_name,
                          singlets = ncol(PBMC_r$raw.data))
    cell.df <- rbind(cell.df,
                     new_add)
  }
}




# JP
for (i in seq(1,5,1)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("JP_RIK_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("JP_RIK_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    PBMC_r <- readRDS(paste0("../",lib_name,
                             "/RCA_PBMC_after_doublet_removal_and_before_QC.rds"))
    PBMC_r <- dataLogNormalise(PBMC_r)
    
    ############ project to all immune panels ################
    # get projection results against global panel
    PBMC_r <- dataProject(PBMC_r, method = "GlobalPanel",
                          corMeth = "pearson", scale = TRUE)
    global.proj <- PBMC_r$projection.data
    global.proj.immune <- read.table("../rownames_of_glocal_projection_immune_cells.txt", 
                                     header = FALSE)
    global.proj <- as.data.frame(as.matrix(global.proj[global.proj.immune$V1,]))
    
    # get projection results against other two panels
    PBMC_r <- dataProjectMultiPanel(PBMC_r,method = list("NovershternPanel", 
                                                         "MonacoPanel",
                                                         "CITESeqPanel"),
                                    scale = TRUE,corMeth = "pearson")
    two.proj <- as.data.frame(as.matrix(PBMC_r$projection.data))
    
    # combine these panels
    proj.all <- rbind(global.proj,two.proj)
    proj.all <- as.matrix(proj.all)
    proj.all <- as(proj.all, "dgCMatrix")
    
    # Assign projection result to RCA object
    PBMC_r$projection.data <- proj.all
    
    #Estimate the most probable cell type label for each cell
    PBMC_r <- estimateCellTypeFromProjection(PBMC_r,confidence = NULL)
    
    saveRDS(PBMC_r,paste0("./",lib_name,".singlets.RDS"))
    
    new_add <- data.frame(lib = lib_name,
                          singlets = ncol(PBMC_r$raw.data))
    cell.df <- rbind(cell.df,
                     new_add)
    
  }
}



# KR
for (i in seq(2,12,1)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("KR_SGI_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("KR_SGI_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    PBMC_r <- readRDS(paste0("../",lib_name,
                             "/RCA_PBMC_after_doublet_removal_and_before_QC.rds"))
    PBMC_r <- dataLogNormalise(PBMC_r)
    
    ############ project to all immune panels ################
    # get projection results against global panel
    PBMC_r <- dataProject(PBMC_r, method = "GlobalPanel",
                          corMeth = "pearson", scale = TRUE)
    global.proj <- PBMC_r$projection.data
    global.proj.immune <- read.table("../rownames_of_glocal_projection_immune_cells.txt", 
                                     header = FALSE)
    global.proj <- as.data.frame(as.matrix(global.proj[global.proj.immune$V1,]))
    
    # get projection results against other two panels
    PBMC_r <- dataProjectMultiPanel(PBMC_r,method = list("NovershternPanel", 
                                                         "MonacoPanel",
                                                         "CITESeqPanel"),
                                    scale = TRUE,corMeth = "pearson")
    two.proj <- as.data.frame(as.matrix(PBMC_r$projection.data))
    
    # combine these panels
    proj.all <- rbind(global.proj,two.proj)
    proj.all <- as.matrix(proj.all)
    proj.all <- as(proj.all, "dgCMatrix")
    
    # Assign projection result to RCA object
    PBMC_r$projection.data <- proj.all
    
    #Estimate the most probable cell type label for each cell
    PBMC_r <- estimateCellTypeFromProjection(PBMC_r,confidence = NULL)
    
    saveRDS(PBMC_r,paste0("./",lib_name,".singlets.RDS"))
    
    new_add <- data.frame(lib = lib_name,
                          singlets = ncol(PBMC_r$raw.data))
    cell.df <- rbind(cell.df,
                     new_add)
    
  }
}

write.table(cell.df,"singlet_number_across_all_lib.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## add SG batch 20 and 21
cell.df <- read.table("singlet_number_across_all_lib.txt", sep = "\t", header = TRUE)
# SG
for (i in c(20,21)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("SG_HEL_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("SG_HEL_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    PBMC_r <- readRDS(paste0("../",lib_name,
                             "/RCA_PBMC_after_doublet_removal_and_before_QC.rds"))
    PBMC_r <- dataLogNormalise(PBMC_r)
    
    ############ project to all immune panels ################
    # get projection results against global panel
    PBMC_r <- dataProject(PBMC_r, method = "GlobalPanel",
                          corMeth = "pearson", scale = TRUE)
    global.proj <- PBMC_r$projection.data
    global.proj.immune <- read.table("../rownames_of_glocal_projection_immune_cells.txt", 
                                     header = FALSE)
    global.proj <- as.data.frame(as.matrix(global.proj[global.proj.immune$V1,]))
    
    # get projection results against other two panels
    PBMC_r <- dataProjectMultiPanel(PBMC_r,method = list("NovershternPanel", 
                                                         "MonacoPanel",
                                                         "CITESeqPanel"),
                                    scale = TRUE,corMeth = "pearson")
    two.proj <- as.data.frame(as.matrix(PBMC_r$projection.data))
    
    # combine these panels
    proj.all <- rbind(global.proj,two.proj)
    proj.all <- as.matrix(proj.all)
    proj.all <- as(proj.all, "dgCMatrix")
    
    # Assign projection result to RCA object
    PBMC_r$projection.data <- proj.all
    
    #Estimate the most probable cell type label for each cell
    PBMC_r <- estimateCellTypeFromProjection(PBMC_r,confidence = NULL)
    
    saveRDS(PBMC_r,paste0("./",lib_name,".singlets.RDS"))
    
    new_add <- data.frame(lib = lib_name,
                          singlets = ncol(PBMC_r$raw.data))
    cell.df <- rbind(cell.df,
                     new_add)
  }
}

write.table(cell.df,"singlet_number_across_all_lib.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## add JP batch 6 and 10
cell.df <- read.table("singlet_number_across_all_lib.txt", sep = "\t", header = TRUE)
for (i in seq(6,10,1)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("JP_RIK_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("JP_RIK_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    PBMC_r <- readRDS(paste0("../",lib_name,
                             "/RCA_PBMC_after_doublet_removal_and_before_QC.rds"))
    PBMC_r <- dataLogNormalise(PBMC_r)
    
    ############ project to all immune panels ################
    # get projection results against global panel
    PBMC_r <- dataProject(PBMC_r, method = "GlobalPanel",
                          corMeth = "pearson", scale = TRUE)
    global.proj <- PBMC_r$projection.data
    global.proj.immune <- read.table("../rownames_of_glocal_projection_immune_cells.txt", 
                                     header = FALSE)
    global.proj <- as.data.frame(as.matrix(global.proj[global.proj.immune$V1,]))
    
    # get projection results against other two panels
    PBMC_r <- dataProjectMultiPanel(PBMC_r,method = list("NovershternPanel", 
                                                         "MonacoPanel",
                                                         "CITESeqPanel"),
                                    scale = TRUE,corMeth = "pearson")
    two.proj <- as.data.frame(as.matrix(PBMC_r$projection.data))
    
    # combine these panels
    proj.all <- rbind(global.proj,two.proj)
    proj.all <- as.matrix(proj.all)
    proj.all <- as(proj.all, "dgCMatrix")
    
    # Assign projection result to RCA object
    PBMC_r$projection.data <- proj.all
    
    #Estimate the most probable cell type label for each cell
    PBMC_r <- estimateCellTypeFromProjection(PBMC_r,confidence = NULL)
    
    saveRDS(PBMC_r,paste0("./",lib_name,".singlets.RDS"))
    
    new_add <- data.frame(lib = lib_name,
                          singlets = ncol(PBMC_r$raw.data))
    cell.df <- rbind(cell.df,
                     new_add)
    
  }
}


write.table(cell.df,"singlet_number_across_all_lib.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
