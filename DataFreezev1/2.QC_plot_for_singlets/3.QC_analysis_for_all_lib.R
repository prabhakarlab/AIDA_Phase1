library(ggplot2)
library(RCAv2)

######################## Doublets and Singlets #################
all.df <- data.frame()
# SG
for (i in c(1,seq(3,21,1))){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("SG_HEL_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("SG_HEL_B00",i,"_L00",j,"_5GEX")
    }
    
    cell.ij <- read.table(paste0("../",lib_name,"/dragen_demultiplexing_type_summary.txt"),
                          sep = "\t", header = FALSE)
    new.df <- data.frame(lib = lib_name,
                         all = sum(cell.ij$V2))
    all.df <- rbind(all.df, new.df)
  }
}
# JP
for (i in seq(1,10,1)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("JP_RIK_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("JP_RIK_B00",i,"_L00",j,"_5GEX")
    }
    
    cell.ij <- read.table(paste0("../",lib_name,"/dragen_demultiplexing_type_summary.txt"),
                          sep = "\t", header = FALSE)
    new.df <- data.frame(lib = lib_name,
                         all = sum(cell.ij$V2))
    all.df <- rbind(all.df, new.df)
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
    
    cell.ij <- read.table(paste0("../",lib_name,"/dragen_demultiplexing_type_summary.txt"),
                          sep = "\t", header = FALSE)
    new.df <- data.frame(lib = lib_name,
                         all = sum(cell.ij$V2))
    all.df <- rbind(all.df, new.df)
  }
}

singlet.df <- read.table("./singlet_number_across_all_lib.txt", 
                         sep = "\t", header = TRUE)
all.df.merge <- merge(x = all.df,
                y = singlet.df,
                by = "lib")
all.df.merge$doublet <- all.df.merge$all - all.df.merge$singlets

singlet.1 <- data.frame(lib = all.df.merge$lib,
                        num = all.df.merge$singlets,
                        type = "Singlet")
doublet.1 <- data.frame(lib = all.df.merge$lib,
                        num = all.df.merge$doublet,
                        type = "Doublet")
all.1 <- rbind(singlet.1, doublet.1)
all.1$lib <- factor(all.1$lib,
                    levels = all.df$lib)

all.1$type <- factor(all.1$type,levels = c("Doublet","Singlet"))


ggplot(all.1, aes(fill=type, y=num, x=lib)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("3.Singlet_Doublet_number_across_all_libraries.pdf",
         width = 15, height = 6)



ggplot(all.1, aes(fill=type, y=num, x=lib)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("3.Singlet_Doublet_percentage_across_all_libraries.pdf",
         width = 15, height = 6)


####################### NODG and pMito of all singlets #########################
all.df.sg <- data.frame()
# SG
for (i in c(1,seq(3,21,1))){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("SG_HEL_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("SG_HEL_B00",i,"_L00",j,"_5GEX")
    }
    
    print(lib_name)
    cell.ij <- readRDS(paste0("../",lib_name,"/RCA_PBMC_after_doublet_removal_and_before_QC.rds"))
    
    rawdata <- cell.ij$raw.data
    nGeneVec <- Matrix::colSums(rawdata>0)
    # Select mito genes
    mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
    # Compute percent.mito vector
    pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)
    
    pMitoVec <- pMitoVec[names(nGeneVec)]
    qc_matrix <- data.frame(barcode = paste(colnames(rawdata), "-SG_B", i,"_L", j, sep = ""),
                            NODG = nGeneVec,
                            pMt = pMitoVec,
                            lib = lib_name)
    
    all.df.sg <- rbind(all.df.sg, qc_matrix)
  }
}


# JP 
all.df.jp <- data.frame()
for (i in seq(1,10,1)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("JP_RIK_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("JP_RIK_B00",i,"_L00",j,"_5GEX")
    }
    
    print(lib_name)
    cell.ij <- readRDS(paste0("../",lib_name,"/RCA_PBMC_after_doublet_removal_and_before_QC.rds"))
    
    rawdata <- cell.ij$raw.data
    nGeneVec <- Matrix::colSums(rawdata>0)
    # Select mito genes
    mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
    # Compute percent.mito vector
    pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)
    
    pMitoVec <- pMitoVec[names(nGeneVec)]
    qc_matrix <- data.frame(barcode = paste(colnames(rawdata), "-JP_B", i,"_L", j, sep = ""),
                            NODG = nGeneVec,
                            pMt = pMitoVec,
                            lib = lib_name)
    
    all.df.jp <- rbind(all.df.jp, qc_matrix)
  }
}


# KR
all.df.kr <- data.frame()
for (i in seq(2,12,1)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("KR_SGI_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("KR_SGI_B00",i,"_L00",j,"_5GEX")
    }
    print(lib_name)
    cell.ij <- readRDS(paste0("../",lib_name,"/RCA_PBMC_after_doublet_removal_and_before_QC.rds"))
    
    rawdata <- cell.ij$raw.data
    nGeneVec <- Matrix::colSums(rawdata>0)
    # Select mito genes
    mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
    # Compute percent.mito vector
    pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)
    
    pMitoVec <- pMitoVec[names(nGeneVec)]
    qc_matrix <- data.frame(barcode = paste(colnames(rawdata), "-KR_B", i,"_L", j, sep = ""),
                            NODG = nGeneVec,
                            pMt = pMitoVec,
                            lib = lib_name)
    
    all.df.kr <- rbind(all.df.kr, qc_matrix)
  }
}

all.df.sg$country <- "SG"
all.df.jp$country <- "JP"
all.df.kr$country <- "KR"

all.df <- rbind(all.df.sg, all.df.jp,all.df.kr )

all.df$lib <- factor(all.df$lib,
                         levels = unique(all.df$lib))

ggplot(all.df, aes(x=lib, y=NODG)) +
  geom_boxplot(aes(fill=country),outlier.shape = NA) + 
  scale_fill_manual(values = c("red","lightblue","lightyellow"))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  xlab("")+ylab("NODG")+ylim(c(0,4000))
ggsave("3.NODG_for_all_samples_across_library.pdf", width = 20, height = 5)


ggplot(all.df, aes(x=lib, y=pMt)) +
  geom_boxplot(aes(fill=country),outlier.shape = NA) + 
  scale_fill_manual(values = c("red","lightblue","lightyellow"))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  xlab("")+ylab("pMt")+ylim(c(0,0.1))
ggsave("3.pMt_for_all_samples_across_library.pdf", width = 20, height = 5)



######################## Donor distriubtion #################
sample.all <- read.delim("./all_singlets_mapped_to_their_samples.txt",
                         header = TRUE, sep = "\t")
donor.df <- sample.all[,c("sample","lib")]
donor.df$count <- 1

donor.df.a <- aggregate(count ~ ., donor.df, FUN = sum)
donor.df.a$country <- unlist(lapply(donor.df.a$lib,
                                    function(x) unlist(strsplit(x,split = "_"))[1] ))
donor.df.a$lib <- factor(donor.df.a$lib,
                         levels = unique(sample.all$lib))

# remove JP_RIK_B001_L002_5GEX
donor.df.a <- donor.df.a[(donor.df.a$lib != "JP_RIK_B001_L002_5GEX"),]
donor.df.a$lib <- factor(donor.df.a$lib,
                         levels = c("SG_HEL_B001_L001_5GEX","SG_HEL_B001_L002_5GEX",
                                    "SG_HEL_B003_L001_5GEX","SG_HEL_B003_L002_5GEX",
                                    "SG_HEL_B004_L001_5GEX","SG_HEL_B004_L002_5GEX",
                                    "SG_HEL_B005_L001_5GEX","SG_HEL_B005_L002_5GEX",
                                    "SG_HEL_B006_L001_5GEX","SG_HEL_B006_L002_5GEX",
                                    "SG_HEL_B007_L001_5GEX","SG_HEL_B007_L002_5GEX",
                                    "SG_HEL_B008_L001_5GEX","SG_HEL_B008_L002_5GEX",
                                    "SG_HEL_B009_L001_5GEX","SG_HEL_B009_L002_5GEX",
                                    "SG_HEL_B010_L001_5GEX","SG_HEL_B010_L002_5GEX",
                                    "SG_HEL_B011_L001_5GEX","SG_HEL_B011_L002_5GEX",
                                    "SG_HEL_B012_L001_5GEX","SG_HEL_B012_L002_5GEX",
                                    "SG_HEL_B013_L001_5GEX","SG_HEL_B013_L002_5GEX",
                                    "SG_HEL_B014_L001_5GEX","SG_HEL_B014_L002_5GEX",
                                    "SG_HEL_B015_L001_5GEX","SG_HEL_B015_L002_5GEX",
                                    "SG_HEL_B016_L001_5GEX","SG_HEL_B016_L002_5GEX",
                                    "SG_HEL_B017_L001_5GEX","SG_HEL_B017_L002_5GEX",
                                    "SG_HEL_B018_L001_5GEX","SG_HEL_B018_L002_5GEX",
                                    "SG_HEL_B019_L001_5GEX","SG_HEL_B019_L002_5GEX",
                                    "SG_HEL_B020_L001_5GEX","SG_HEL_B020_L002_5GEX",
                                    "SG_HEL_B021_L001_5GEX","SG_HEL_B021_L002_5GEX",
                                    "JP_RIK_B001_L001_5GEX",
                                    "JP_RIK_B002_L001_5GEX","JP_RIK_B002_L002_5GEX",
                                    "JP_RIK_B003_L001_5GEX","JP_RIK_B003_L002_5GEX",
                                    "JP_RIK_B004_L001_5GEX","JP_RIK_B004_L002_5GEX",
                                    "JP_RIK_B005_L001_5GEX","JP_RIK_B005_L002_5GEX",
                                    "JP_RIK_B006_L001_5GEX","JP_RIK_B006_L002_5GEX",
                                    "JP_RIK_B007_L001_5GEX","JP_RIK_B007_L002_5GEX",
                                    "JP_RIK_B008_L001_5GEX","JP_RIK_B008_L002_5GEX",
                                    "JP_RIK_B009_L001_5GEX","JP_RIK_B009_L002_5GEX",
                                    "JP_RIK_B010_L001_5GEX","JP_RIK_B010_L002_5GEX",
                                    "KR_SGI_B001_L001_5GEX","KR_SGI_B001_L002_5GEX",
                                    "KR_SGI_B002_L001_5GEX","KR_SGI_B002_L002_5GEX",
                                    "KR_SGI_B003_L001_5GEX","KR_SGI_B003_L002_5GEX",
                                    "KR_SGI_B004_L001_5GEX","KR_SGI_B004_L002_5GEX",
                                    "KR_SGI_B005_L001_5GEX","KR_SGI_B005_L002_5GEX",
                                    "KR_SGI_B006_L001_5GEX","KR_SGI_B006_L002_5GEX",
                                    "KR_SGI_B007_L001_5GEX","KR_SGI_B007_L002_5GEX",
                                    "KR_SGI_B008_L001_5GEX","KR_SGI_B008_L002_5GEX",
                                    "KR_SGI_B009_L001_5GEX","KR_SGI_B009_L002_5GEX",
                                    "KR_SGI_B010_L001_5GEX","KR_SGI_B010_L002_5GEX",
                                    "KR_SGI_B011_L001_5GEX","KR_SGI_B011_L002_5GEX",
                                    "KR_SGI_B012_L001_5GEX","KR_SGI_B012_L002_5GEX"))

ggplot(donor.df.a, aes(x=lib, y=count)) +
  geom_boxplot(aes(fill=country),outlier.shape = NA) + 
  geom_jitter(height=0, stroke=0, size=1, width=0.25)+
  scale_fill_manual(values = c("red","lightblue","lightyellow"))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  xlab("")+ylab("number")+ylim(c(0,7500))+
  geom_hline(yintercept = 1000, color="red")
ggsave("3.cell_number_for_all_samples_across_library.pdf", width = 40, height = 5)



ggplot(donor.df.a, aes(x=lib, y=count)) +
  geom_boxplot(aes(fill=country),outlier.shape = NA) + 
  geom_jitter(height=0, stroke=0, size=1, width=0.25)+
  scale_fill_manual(values = c("red","lightblue","lightyellow"))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  xlab("")+ylab("number")+ylim(c(0,2000))+
  geom_hline(yintercept = 1000, color="red")
ggsave("3.cell_number_for_all_samples_across_library_zoom_in.pdf", width = 40, height = 5)


