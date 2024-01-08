library(RCAv2)
library(ggplot2)

sample.df <- data.frame()
# SG
for (i in c(1,seq(3,19,1))){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("SG_HEL_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("SG_HEL_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    PBMC_r <- readRDS(paste0("./",lib_name,".singlets.RDS"))
    
    # if there is missing sample
    if (i %in% c(1,8,10,12,20)){
      demux.file <- list.files(path = paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                                             lib_name,"/demultiplexing/"),
                               pattern = "\\.missingSampleAdded.tsv$")
    } else {
      demux.file <- list.files(path = paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                                             lib_name,"/demultiplexing/"),
                               pattern = "\\.barcodeSummary.tsv$")
    }
    
    demux <- read.delim(paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                               lib_name,"/demultiplexing/",demux.file),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    demux.singlet <- demux[(demux$Barcode %in% colnames(PBMC_r$raw.data)),]
    demux.singlet$type <- unlist(lapply(demux.singlet$SampleIdentity, 
                                        function(x) unlist(strsplit(x,split = ":"))[1] ))
    demux.singlet$sample <- unlist(lapply(demux.singlet$SampleIdentity, 
                                        function(x) unlist(strsplit(x,split = ":"))[2] ))
    print(unique(demux.singlet$type))
    
    demux.singlet.sub <- demux.singlet[,c("Barcode","sample")]
    demux.singlet.sub$Barcode <- paste(demux.singlet.sub$Barcode, "-SG_B", i,"_L", j, sep = "")
    demux.singlet.sub$lib <- lib_name
    
    sample.df <- rbind(sample.df, demux.singlet.sub)

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
    PBMC_r <- readRDS(paste0("./",lib_name,".singlets.RDS"))
    
    demux.file <- list.files(path = paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                                           lib_name,"/demultiplexing/"),
                             pattern = "\\.barcodeSummary.tsv$")
    
    demux <- read.delim(paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                               lib_name,"/demultiplexing/",demux.file),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    demux.singlet <- demux[(demux$Barcode %in% colnames(PBMC_r$raw.data)),]
    demux.singlet$type <- unlist(lapply(demux.singlet$SampleIdentity, 
                                        function(x) unlist(strsplit(x,split = ":"))[1] ))
    demux.singlet$sample <- unlist(lapply(demux.singlet$SampleIdentity, 
                                          function(x) unlist(strsplit(x,split = ":"))[2] ))
    print(unique(demux.singlet$type))
    
    demux.singlet.sub <- demux.singlet[,c("Barcode","sample")]
    demux.singlet.sub$Barcode <- paste(demux.singlet.sub$Barcode, "-JP_B", i,"_L", j, sep = "")
    demux.singlet.sub$lib <- lib_name
    
    sample.df <- rbind(sample.df, demux.singlet.sub)
    
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
    PBMC_r <- readRDS(paste0("./",lib_name,".singlets.RDS"))
    
    demux.file <- list.files(path = paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                                           lib_name,"/demultiplexing/"),
                             pattern = "\\.barcodeSummary.tsv$")
    
    demux <- read.delim(paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                               lib_name,"/demultiplexing/",demux.file),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    demux.singlet <- demux[(demux$Barcode %in% colnames(PBMC_r$raw.data)),]
    demux.singlet$type <- unlist(lapply(demux.singlet$SampleIdentity, 
                                        function(x) unlist(strsplit(x,split = ":"))[1] ))
    demux.singlet$sample <- unlist(lapply(demux.singlet$SampleIdentity, 
                                          function(x) unlist(strsplit(x,split = ":"))[2] ))
    print(unique(demux.singlet$type))
    
    demux.singlet.sub <- demux.singlet[,c("Barcode","sample")]
    demux.singlet.sub$Barcode <- paste(demux.singlet.sub$Barcode, "-KR_B", i,"_L", j, sep = "")
    demux.singlet.sub$lib <- lib_name
    
    sample.df <- rbind(sample.df, demux.singlet.sub)
    
  }
}

write.table(sample.df, "all_singlets_mapped_to_their_samples.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


############ add SG batch 20, 21
sample.df <- read.table("all_singlets_mapped_to_their_samples.txt",sep = "\t", header = TRUE)

sample.df.sg <- data.frame()
# SG
for (i in c(20,21)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("SG_HEL_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("SG_HEL_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    PBMC_r <- readRDS(paste0("./",lib_name,".singlets.RDS"))
    
    # if there is missing sample
    if (i %in% c(1,8,10,12,20)){
      demux.file <- list.files(path = paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                                             lib_name,"/demultiplexing/"),
                               pattern = "\\.missingSampleAdded.tsv$")
    } else {
      demux.file <- list.files(path = paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                                             lib_name,"/demultiplexing/"),
                               pattern = "\\.barcodeSummary.tsv$")
    }
    
    demux <- read.delim(paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                               lib_name,"/demultiplexing/",demux.file),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    demux.singlet <- demux[(demux$Barcode %in% colnames(PBMC_r$raw.data)),]
    demux.singlet$type <- unlist(lapply(demux.singlet$SampleIdentity, 
                                        function(x) unlist(strsplit(x,split = ":"))[1] ))
    demux.singlet$sample <- unlist(lapply(demux.singlet$SampleIdentity, 
                                          function(x) unlist(strsplit(x,split = ":"))[2] ))
    print(unique(demux.singlet$type))
    
    demux.singlet.sub <- demux.singlet[,c("Barcode","sample")]
    demux.singlet.sub$Barcode <- paste(demux.singlet.sub$Barcode, "-SG_B", i,"_L", j, sep = "")
    demux.singlet.sub$lib <- lib_name
    
    sample.df.sg <- rbind(sample.df.sg, demux.singlet.sub)
    
  }
}
sample.df <- rbind(sample.df, sample.df.sg)
write.table(sample.df, "all_singlets_mapped_to_their_samples.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

###### plot cell viability and number
sample.df <- read.table("all_singlets_mapped_to_their_samples.txt",
                        sep = "\t", header = TRUE)
sample.df.sg <- sample.df[grep("^SG_HEL", sample.df$lib),]

cell.via <- read.table("../metadata_files/SG_samples_and_viability_info.txt", 
                       header = TRUE, sep = "\t")
sample.df.sg <- sample.df.sg[which(sample.df.sg$sample %in% cell.via$sample),]

sample.df.sg$count <- 1
sample.df.sg$Barcode <- NULL
sample.df.sg$lib <- NULL

sample.df.sg.sum <- aggregate(count ~ ., sample.df.sg, FUN = sum)

cell.via <- cell.via[(cell.via$sample %in% sample.df.sg.sum$sample),]
sample.via <- merge(x = sample.df.sg.sum, y = cell.via, by="sample")

ggplot()+
  geom_point(sample.via, mapping = aes(x = count, y=viability))+
  theme_bw(11)
ggsave("2.SG_sample_number_vs_cell_viability.pdf")




############ add JP batch 6-10
sample.df <- read.table("all_singlets_mapped_to_their_samples.txt",sep = "\t", header = TRUE)

sample.df.jp <- data.frame()
for (i in seq(6,10,1)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("JP_RIK_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("JP_RIK_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    PBMC_r <- readRDS(paste0("./",lib_name,".singlets.RDS"))
    
    demux.file <- list.files(path = paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                                           lib_name,"/demultiplexing/"),
                             pattern = "\\.barcodeSummary.tsv$")
    
    demux <- read.delim(paste0("/mnt/quy_1TB/AIDA_project/DRAGEN_pipeline/dragen_output/",
                               lib_name,"/demultiplexing/",demux.file),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    demux.singlet <- demux[(demux$Barcode %in% colnames(PBMC_r$raw.data)),]
    demux.singlet$type <- unlist(lapply(demux.singlet$SampleIdentity, 
                                        function(x) unlist(strsplit(x,split = ":"))[1] ))
    demux.singlet$sample <- unlist(lapply(demux.singlet$SampleIdentity, 
                                          function(x) unlist(strsplit(x,split = ":"))[2] ))
    print(unique(demux.singlet$type))
    
    demux.singlet.sub <- demux.singlet[,c("Barcode","sample")]
    demux.singlet.sub$Barcode <- paste(demux.singlet.sub$Barcode, "-JP_B", i,"_L", j, sep = "")
    demux.singlet.sub$lib <- lib_name
    
    sample.df.jp <- rbind(sample.df.jp, demux.singlet.sub)
    
  }
}
sample.df <- rbind(sample.df, sample.df.jp)
write.table(sample.df, "all_singlets_mapped_to_their_samples.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


