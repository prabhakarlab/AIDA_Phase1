library(dplyr)
source("./utils/quy_BCR_utils.R")

valid_barcode_all_sg <- data.frame()

# SG
for (i in c(1,seq(3,19,1))){
  print(i)
  for (j in c(1,2)){
    if (i < 10){
      TCR_path <- paste0("/mnt/quy_1TB/AIDA_project/cellranger5_vdj_output/BCR/SG_B00",i,"_L00",
                         j,"_BCR-filtered_contig_annotations.csv")
    } else {
      TCR_path <- paste0("/mnt/quy_1TB/AIDA_project/cellranger5_vdj_output/BCR/SG_B0",i,"_L00",
                         j,"_BCR-filtered_contig_annotations.csv")
    }
    
    
    tcr.res <- read.csv(TCR_path, stringsAsFactors = FALSE)
    tcr.res.combined <- BCRprofle(tcr.res)
    tcr.res.combined$barcode <- unlist(lapply(tcr.res.combined$barcode,
                                              function(x) unlist(strsplit(x,split = "-"))[1] ))
    tcr.res.combined$barcode <- paste(tcr.res.combined$barcode, "-SG_B", i,"_L", j, sep = "")
    valid_barcode_all_sg <- rbind(valid_barcode_all_sg, tcr.res.combined)
  }
  
}


valid_barcode_all_kr <- data.frame()
# KR
for (i in seq(2,12,1)){
  print(i)
  for (j in c(1,2)){
    if (i < 10){
      TCR_path <- paste0("/mnt/quy_1TB/AIDA_project/cellranger5_vdj_output/BCR/KR_B00",i,"_L00",
                         j,"_BCR-filtered_contig_annotations.csv")
    } else {
      TCR_path <- paste0("/mnt/quy_1TB/AIDA_project/cellranger5_vdj_output/BCR/KR_B0",i,"_L00",
                         j,"_BCR-filtered_contig_annotations.csv")
    }
    
    
    tcr.res <- read.csv(TCR_path, stringsAsFactors = FALSE)
    tcr.res.combined <- BCRprofle(tcr.res)
    tcr.res.combined$barcode <- unlist(lapply(tcr.res.combined$barcode,
                                              function(x) unlist(strsplit(x,split = "-"))[1] ))
    tcr.res.combined$barcode <- paste(tcr.res.combined$barcode, "-KR_B", i,"_L", j, sep = "")
    valid_barcode_all_kr <- rbind(valid_barcode_all_kr, tcr.res.combined)
  }
  
}


valid_barcode_all_jp <- data.frame()
# JP
for (i in seq(1,5,1)){
  print(i)
  for (j in c(1,2)){
    # skip B001 L002
    if(i==1 & j==2){
      next
    }
    
    if (i < 10){
      TCR_path <- paste0("/mnt/quy_1TB/AIDA_project/cellranger5_vdj_output/BCR/JP_B00",i,"_L00",
                         j,"_BCR-filtered_contig_annotations.csv")
    } else {
      TCR_path <- paste0("/mnt/quy_1TB/AIDA_project/cellranger5_vdj_output/BCR/JP_B0",i,"_L00",
                         j,"_BCR-filtered_contig_annotations.csv")
    }
    
    
    tcr.res <- read.csv(TCR_path, stringsAsFactors = FALSE)
    tcr.res.combined <- BCRprofle(tcr.res)
    tcr.res.combined$barcode <- unlist(lapply(tcr.res.combined$barcode,
                                              function(x) unlist(strsplit(x,split = "-"))[1] ))
    tcr.res.combined$barcode <- paste(tcr.res.combined$barcode, "-JP_B", i,"_L", j, sep = "")
    valid_barcode_all_jp <- rbind(valid_barcode_all_jp, tcr.res.combined)
  }
  
}
valid_barcode_all <- rbind(valid_barcode_all_sg,valid_barcode_all_kr,valid_barcode_all_jp)

write.table(valid_barcode_all, "valid_BCR_across_all_batch.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

##### add JP batch 6-10
valid_barcode_all <- read.table("valid_BCR_across_all_batch.txt",
                                sep = "\t", header = TRUE)

valid_barcode_all_jp <- data.frame()
# JP
for (i in seq(6,10,1)){
  print(i)
  for (j in c(1,2)){
    # skip B001 L002
    if(i==1 & j==2){
      next
    }
    
    if (i < 10){
      TCR_path <- paste0("/mnt/quy_1TB/AIDA_project/cellranger5_vdj_output/BCR/JP_B00",i,"_L00",
                         j,"_BCR-filtered_contig_annotations.csv")
    } else {
      TCR_path <- paste0("/mnt/quy_1TB/AIDA_project/cellranger5_vdj_output/BCR/JP_B0",i,"_L00",
                         j,"_BCR-filtered_contig_annotations.csv")
    }
    
    
    tcr.res <- read.csv(TCR_path, stringsAsFactors = FALSE)
    tcr.res.combined <- BCRprofle(tcr.res)
    tcr.res.combined$barcode <- unlist(lapply(tcr.res.combined$barcode,
                                              function(x) unlist(strsplit(x,split = "-"))[1] ))
    tcr.res.combined$barcode <- paste(tcr.res.combined$barcode, "-JP_B", i,"_L", j, sep = "")
    valid_barcode_all_jp <- rbind(valid_barcode_all_jp, tcr.res.combined)
  }
  
}

valid_barcode_all <- rbind(valid_barcode_all, valid_barcode_all_jp)

write.table(valid_barcode_all, "valid_BCR_across_all_batch.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



