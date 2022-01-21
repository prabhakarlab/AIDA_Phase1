library(RCAv2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rstatix)


all.cell.df <- read.table("../combined_RCA_QC/4.cell_annotation_for_good_cells_after_QC_healthy_only.txt",
                          header = TRUE, sep = "\t")
all.cell <- all.cell.df$barcode

PBMC.list <- list()

# SG
for (i in c(1,seq(3,21,1))){
  if (i >= 10){
    lib_name1 <- paste0("SG_HEL_B0",i,"_L001_5GEX")
    lib_name2 <- paste0("SG_HEL_B0",i,"_L002_5GEX")
  } else {
    lib_name1 <- paste0("SG_HEL_B00",i,"_L001_5GEX")
    lib_name2 <- paste0("SG_HEL_B00",i,"_L002_5GEX")
  }
  
  print(lib_name1)
  PBMC_r1 <- readRDS(paste0("../all_singlet_RDS/",lib_name1,
                            ".singlets.RDS"))
  PBMC_r2 <- readRDS(paste0("../all_singlet_RDS/",lib_name2,
                            ".singlets.RDS"))
  
  raw.data1 <- PBMC_r1$raw.data
  colnames(raw.data1) <- paste(colnames(raw.data1), "-SG_B", i,"_L1", sep = "")
  cell1 <- colnames(raw.data1)
  cell1.good <- cell1[(cell1 %in% all.cell)]
  raw.data1 <- raw.data1[,cell1.good]
  
  
  raw.data2 <- PBMC_r2$raw.data
  colnames(raw.data2) <- paste(colnames(raw.data2), "-SG_B", i,"_L2", sep = "")
  cell2 <- colnames(raw.data2)
  cell2.good <- cell2[(cell2 %in% all.cell)]
  raw.data2 <- raw.data2[,cell2.good]
  
  raw.data.all <- Seurat::RowMergeSparseMatrices(raw.data1, raw.data2)
  
  #Generate a Seurat object
  seu.i <- CreateSeuratObject(counts = raw.data.all, 
                              min.cells = 0, 
                              min.features = 0)
  name.i <- paste0("SG_",i)
  PBMC.list[[name.i]] <- seu.i
}


# JP
## only keep library 1 for batch 1
PBMC_r1 <- readRDS("../all_singlet_RDS/JP_RIK_B001_L001_5GEX.singlets.RDS")
raw.data1 <- PBMC_r1$raw.data
colnames(raw.data1) <- paste(colnames(raw.data1), "-JP_B1_L1", sep = "")
cell1 <- colnames(raw.data1)
cell1.good <- cell1[(cell1 %in% all.cell)]
raw.data1 <- raw.data1[,cell1.good]
#Generate a Seurat object
seu.i <- CreateSeuratObject(counts = raw.data1, 
                            min.cells = 0, 
                            min.features = 0)
name.i <- "JP_1"
PBMC.list[[name.i]] <- seu.i


for (i in seq(2,10,1)){
  if (i >= 10){
    lib_name1 <- paste0("JP_RIK_B0",i,"_L001_5GEX")
    lib_name2 <- paste0("JP_RIK_B0",i,"_L002_5GEX")
  } else {
    lib_name1 <- paste0("JP_RIK_B00",i,"_L001_5GEX")
    lib_name2 <- paste0("JP_RIK_B00",i,"_L002_5GEX")
  }
  
  print(lib_name1)
  PBMC_r1 <- readRDS(paste0("../all_singlet_RDS/",lib_name1,
                            ".singlets.RDS"))
  PBMC_r2 <- readRDS(paste0("../all_singlet_RDS/",lib_name2,
                            ".singlets.RDS"))
  
  raw.data1 <- PBMC_r1$raw.data
  colnames(raw.data1) <- paste(colnames(raw.data1), "-JP_B", i,"_L1", sep = "")
  cell1 <- colnames(raw.data1)
  cell1.good <- cell1[(cell1 %in% all.cell)]
  raw.data1 <- raw.data1[,cell1.good]
  
  raw.data2 <- PBMC_r2$raw.data
  colnames(raw.data2) <- paste(colnames(raw.data2), "-JP_B", i,"_L2", sep = "")
  cell2 <- colnames(raw.data2)
  cell2.good <- cell2[(cell2 %in% all.cell)]
  raw.data2 <- raw.data2[,cell2.good]
  
  raw.data.all <- Seurat::RowMergeSparseMatrices(raw.data1, raw.data2)
  
  #Generate a Seurat object
  seu.i <- CreateSeuratObject(counts = raw.data.all, 
                              min.cells = 0, 
                              min.features = 0)
  name.i <- paste0("JP_",i)
  PBMC.list[[name.i]] <- seu.i
}


# KR
for (i in seq(2,12,1)){
  if (i >= 10){
    lib_name1 <- paste0("KR_SGI_B0",i,"_L001_5GEX")
    lib_name2 <- paste0("KR_SGI_B0",i,"_L002_5GEX")
  } else {
    lib_name1 <- paste0("KR_SGI_B00",i,"_L001_5GEX")
    lib_name2 <- paste0("KR_SGI_B00",i,"_L002_5GEX")
  }
  
  print(lib_name1)
  PBMC_r1 <- readRDS(paste0("../all_singlet_RDS/",lib_name1,
                            ".singlets.RDS"))
  PBMC_r2 <- readRDS(paste0("../all_singlet_RDS/",lib_name2,
                            ".singlets.RDS"))
  
  raw.data1 <- PBMC_r1$raw.data
  colnames(raw.data1) <- paste(colnames(raw.data1), "-KR_B", i,"_L1", sep = "")
  cell1 <- colnames(raw.data1)
  cell1.good <- cell1[(cell1 %in% all.cell)]
  raw.data1 <- raw.data1[,cell1.good]
  
  raw.data2 <- PBMC_r2$raw.data
  colnames(raw.data2) <- paste(colnames(raw.data2), "-KR_B", i,"_L2", sep = "")
  cell2 <- colnames(raw.data2)
  cell2.good <- cell2[(cell2 %in% all.cell)]
  raw.data2 <- raw.data2[,cell2.good]
  
  raw.data.all <- Seurat::RowMergeSparseMatrices(raw.data1, raw.data2)
  
  #Generate a Seurat object
  seu.i <- CreateSeuratObject(counts = raw.data.all, 
                              min.cells = 0, 
                              min.features = 0)
  name.i <- paste0("KR_",i)
  PBMC.list[[name.i]] <- seu.i
}

saveRDS(PBMC.list, "PBMC.list.RDS")


for (i in seq(1,41,1)){
  PBMC.list.i <- PBMC.list[[i]]
  message( paste0(names(PBMC.list)[i], ": ",ncol(PBMC.list.i)) )
}

#normalize and select features
PBMC.list <- lapply(X = PBMC.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})


# select features for downstream integration, and run PCA on each object in the list
features <- SelectIntegrationFeatures(object.list = PBMC.list)
PBMC.list <- lapply(X = PBMC.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


# using KR_7 as reference (most abundant)
reference_dataset <- which(names(PBMC.list) == "KR_7")
anchors <- FindIntegrationAnchors(object.list = PBMC.list, 
                                  reference = reference_dataset, 
                                  reduction = "rpca", 
                                  dims = 1:30)
saveRDS(anchors,"anchors.for.reciprocal.pca.rds")

# integrate data
PBMC.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

saveRDS(PBMC.integrated,"PBMC.integrated.after.QC.rds")


DefaultAssay(PBMC.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
PBMC.integrated <- ScaleData(PBMC.integrated)
PBMC.integrated <- RunPCA(PBMC.integrated, npcs = 30)

pdf("1.Elbowplot.pdf")
ElbowPlot(PBMC.integrated, ndims =30)
dev.off()


#UMAP
PBMC.integrated <- RunUMAP(PBMC.integrated, reduction = "pca", dims = 1:17)

# clustering
PBMC.integrated <- FindNeighbors(PBMC.integrated, dims = 1:17)
PBMC.integrated <- FindClusters(PBMC.integrated, resolution = 2)

umap.df <- FetchData(PBMC.integrated, vars = c("UMAP_1","UMAP_2"))
umap.df$barcode <- rownames(umap.df)
write.table(umap.df,"umap.df.txt",col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

PBMC.integrated <- FindNeighbors(object = PBMC.integrated, 
                                 reduction = "pca", dims = 1:17, 
                                 k.param = 500, verbose = TRUE,compute.SNN = FALSE)
saveRDS(PBMC.integrated,"PBMC.integrated.after.QC.rds")


graph_nnMat <- Matrix::Matrix(PBMC.integrated@meta.data[,"integrated_nn"], sparse = TRUE)
saveRDS(graph_nnMat,"graph_nnMat.rds")

########### plot UMAP
umap.df <- read.table("umap.df.txt", sep = "\t", header = TRUE)
all.cell.df <- read.table("../combined_RCA_QC/4.cell_annotation_for_good_cells_after_QC_healthy_only.txt",
                          header = TRUE, sep = "\t")

all.cell.df.sub <- all.cell.df[,c("barcode","cell","sample")]
umap.df <- merge(x=umap.df,y=all.cell.df.sub, by="barcode")

write.table(umap.df,"umap.df.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

ggplot(data = umap.df, 
       mapping = aes(x = UMAP_1, y = UMAP_2, color=cell)) + 
  geom_point(size = .2, stroke=0) + 
  scale_color_manual(values=c("blue","brown","darkgreen",
                              "darkorange","red","pink",
                              "turquoise","steelblue")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("1.integrated_UMAP_with_annotation_after_QC_and_healthy_only.png",width=9,height=7)


########## plot TCR and BCR
umap.df <- read.table("umap.df.txt", sep = "\t", header = TRUE)
valid_barcode_tcr <- read.table("../TCR_BCR/valid_TCR_across_all_batch.txt",
                                header = TRUE, sep = "\t")
valid_barcode_bcr <- read.table("../TCR_BCR/valid_BCR_across_all_batch.txt",
                                header = TRUE, sep = "\t")
umap.df$TCR <- 0
umap.df$BCR <- 0
umap.df[(umap.df$barcode %in% valid_barcode_tcr$barcode),"TCR"] <- 1
umap.df[(umap.df$barcode %in% valid_barcode_bcr$barcode),"BCR"] <- 1
write.table(umap.df,"umap.df.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

ggplot() + 
  geom_point(data = umap.df, 
             mapping = aes(x = UMAP_1, y = UMAP_2),
             size = .2, stroke=0,col="grey") + 
  geom_point(data = umap.df[which(umap.df$TCR==1),], 
             mapping = aes(x = UMAP_1, y = UMAP_2),
             size = .2, stroke=0,col="red") +
  geom_point(data = umap.df[which(umap.df$BCR==1),], 
             mapping = aes(x = UMAP_1, y = UMAP_2),
             size = .2, stroke=0,col="blue") +
  theme_bw(10) 
ggsave("1.integrated_UMAP_highlighting_TCR_and_BCR.png",width=7,height=7)


##### plot different country
umap.df <- read.table("umap.df.txt", sep = "\t", header = TRUE)
umap.df$lib <- unlist(lapply(umap.df$barcode,
                             function(x) unlist(strsplit(x,split = "-"))[2] ))
umap.df$country <- unlist(lapply(umap.df$lib,
                                 function(x) unlist(strsplit(x,split = "_"))[1] ))
write.table(umap.df,"umap.df.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

unique.country <- c("JP","KR","SG")
color.country <- c("steelblue","darkgreen","red")
for (i in seq_along(unique.country)){
  unique.country.i <- unique.country[i]
  
  ggplot()+
    geom_point(umap.df, mapping = aes(x=UMAP_1,y=UMAP_2), color="grey",size=.15,stroke=0)+
    geom_point(umap.df[(umap.df$country==unique.country.i),], 
               mapping = aes(x=UMAP_1,y=UMAP_2), color=color.country[i],size=.15,stroke=0)+
    theme_bw(10)+
    ggtitle(paste0("n = ",nrow(umap.df[(umap.df$country==unique.country.i),])))
  ggsave(paste0("2.UMAP_highlighting_country_",unique.country.i,".png"),
           device = "png", width = 5, height = 5)
  
}
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
########## remove sample CD_21_11570_XN_HS_ARR in JP batch 10 who has withdrawn from AIDA
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

umap.df <- read.table("umap.df.txt", sep = "\t", header = TRUE)
umap.df.consent <- umap.df[(umap.df$sample != "CD_21_11570_XN_HS_ARR"),]
write.table(umap.df.consent,"umap.df.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
##########################################################################################
########## discard one Japan sample and recluster and visulization
##########################################################################################
umap.df <- read.table("umap.df.txt", sep="\t", header=TRUE)
PBMC.integrated.final <- subset(PBMC.integrated, cells=umap.df$barcode)

#UMAP
PBMC.integrated.final <- RunUMAP(PBMC.integrated.final, reduction = "pca", dims = 1:17)

# clustering
PBMC.integrated.final <- FindNeighbors(PBMC.integrated.final, dims = 1:17)
PBMC.integrated.final <- FindClusters(PBMC.integrated.final, resolution = 2)

umap.df.final <- FetchData(PBMC.integrated.final, vars = c("UMAP_1","UMAP_2"))
umap.df.final$barcode <- rownames(umap.df.final)
write.table(umap.df.final,"umap.df.final.txt",
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

saveRDS(PBMC.integrated.final,"PBMC.integrated.after.QC.and.removing.Japan.withdrawn.sample.rds")

### add meta data 
umap.df.final <- read.table("umap.df.final.txt", header=TRUE, sep="\t")

rownames(umap.df.final) <- umap.df.final$barcode
umap.df.final <- umap.df.final[(colnames(PBMC.integrated.final)),]
PBMC.integrated.final$cell.annotation <- umap.df.final$cell
PBMC.integrated.final$lib <- umap.df.final$lib
PBMC.integrated.final$country <- umap.df.final$country
PBMC.integrated.final$DCP_ID <- umap.df.final$DCP_ID
PBMC.integrated.final$age  <- umap.df.final$age
PBMC.integrated.final$sex  <- umap.df.final$sex
PBMC.integrated.final$ethnicity  <- umap.df.final$ethnicity
saveRDS(PBMC.integrated.final,"PBMC.integrated.after.QC.and.removing.Japan.withdrawn.sample.withMetadata.rds")



########### plot UMAP
umap.df.final <- read.table("umap.df.final.txt", sep = "\t", header = TRUE)
all.cell.df <- read.table("../combined_RCA_QC/4.cell_annotation_for_good_cells_after_QC_healthy_only.txt",
                          header = TRUE, sep = "\t")

all.cell.df.sub <- all.cell.df[(all.cell.df$barcode %in% umap.df.final$barcode),
                               c("barcode","cell","sample")]
umap.df.final <- merge(x=umap.df.final,y=all.cell.df.sub, by="barcode")

write.table(umap.df.final,"umap.df.final.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

ggplot(data = umap.df.final, 
       mapping = aes(x = UMAP_1, y = UMAP_2, color=cell)) + 
  geom_point(size = .2, stroke=0) + 
  scale_color_manual(values=c("blue","brown","darkgreen",
                              "darkorange","red","pink",
                              "turquoise","steelblue")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("1.integrated_UMAP_with_annotation_after_QC_and_healthy_only_discarded_one_JP_sample.png",width=9,height=7)


########## plot TCR and BCR
umap.df.final <- read.table("umap.df.final.txt", sep = "\t", header = TRUE)
valid_barcode_tcr <- read.table("../TCR_BCR/valid_TCR_across_all_batch.txt",
                                header = TRUE, sep = "\t")
valid_barcode_bcr <- read.table("../TCR_BCR/valid_BCR_across_all_batch.txt",
                                header = TRUE, sep = "\t")
umap.df.final$TCR <- 0
umap.df.final$BCR <- 0
umap.df.final[(umap.df.final$barcode %in% valid_barcode_tcr$barcode),"TCR"] <- 1
umap.df.final[(umap.df.final$barcode %in% valid_barcode_bcr$barcode),"BCR"] <- 1
write.table(umap.df.final,"umap.df.final.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

ggplot() + 
  geom_point(data = umap.df.final, 
             mapping = aes(x = UMAP_1, y = UMAP_2),
             size = .2, stroke=0,col="grey") + 
  geom_point(data = umap.df.final[which(umap.df.final$TCR==1),], 
             mapping = aes(x = UMAP_1, y = UMAP_2),
             size = .2, stroke=0,col="red") +
  geom_point(data = umap.df.final[which(umap.df.final$BCR==1),], 
             mapping = aes(x = UMAP_1, y = UMAP_2),
             size = .2, stroke=0,col="blue") +
  theme_bw(10) 
ggsave("1.integrated_UMAP_highlighting_TCR_and_BCR_discarded_one_JP_sample.png",width=7,height=7)


##### plot different country
umap.df.final <- read.table("umap.df.final.txt", sep = "\t", header = TRUE)
umap.df.final$lib <- unlist(lapply(umap.df.final$barcode,
                             function(x) unlist(strsplit(x,split = "-"))[2] ))
umap.df.final$country <- unlist(lapply(umap.df.final$lib,
                                 function(x) unlist(strsplit(x,split = "_"))[1] ))
write.table(umap.df.final,"umap.df.final.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

unique.country <- c("JP","KR","SG")
color.country <- c("steelblue","darkgreen","red")
for (i in seq_along(unique.country)){
  unique.country.i <- unique.country[i]
  
  ggplot()+
    geom_point(umap.df.final, mapping = aes(x=UMAP_1,y=UMAP_2), color="grey",size=.15,stroke=0)+
    geom_point(umap.df.final[(umap.df.final$country==unique.country.i),], 
               mapping = aes(x=UMAP_1,y=UMAP_2), color=color.country[i],size=.15,stroke=0)+
    theme_bw(10)+
    ggtitle(paste0("n = ",nrow(umap.df.final[(umap.df.final$country==unique.country.i),])))
  ggsave(paste0("2.UMAP_highlighting_country_",unique.country.i,"_discarded_one_JP_sample.png"),
         device = "png", width = 5, height = 5)
  
}



################### add sample
umap.df.final <- read.table("umap.df.final.txt", sep = "\t", header = TRUE)
umap.df.final[(umap.df.final$sample == "SG_HEL_H007"), "sample"] <- "SG_HEL_H07a"

sample.df <- read.table("../metadata_files/sample_age_sex_ethnicity_JP1-10_SG1-21_KR2-12.txt", header = TRUE,
                        sep = "\t")
sample.df[(sample.df$sex == "F") , "sex"] <- "Female"
sample.df[(sample.df$sex == "M") , "sex"] <- "Male"

unique_sample <- unique(umap.df.final$sample)
unique_sample <- unique_sample[(unique_sample %in% sample.df$sample)]

sample.df.sub <- sample.df[(sample.df$sample %in% unique_sample),]
x <- data.frame(table(sample.df.sub$sample))
duplicated_lonza <- c("CD_21_09271_XN_HS_ARR","CD_21_09240_XN_HS_ARR",
                      "CD_21_09193_XN_HS_ARR","CD_21_09209_XN_HS_ARR")
sample.df.sub <- sample.df.sub[!(sample.df.sub$sample %in% duplicated_lonza),]
unique_lonza <- data.frame(sample = duplicated_lonza,
                           DCP_ID = c("LONZA3038099","LONZA3038097",
                                      "LONZA3038016","LONZA3038306"),
                           age = c(23,41,21,30),
                           sex = c("Male","Female","Male","Male"),
                           ethnicity = rep("Caucasian",length(duplicated_lonza)))
sample.df.sub <- rbind(sample.df.sub,unique_lonza)

### change lonza ID 
sample.df.sub[(grep("_LONZA", sample.df.sub$DCP_ID)), "DCP_ID"] <- unlist(lapply( sample.df.sub[(grep("_LONZA", sample.df.sub$DCP_ID)), "DCP_ID"],
                                                                                  function(x) unlist(strsplit(x,split = "_"))[3]  ) )

umap.df.final.m <- merge(x = umap.df.final, y = sample.df.sub, by = "sample")

write.table(umap.df.final.m,"umap.df.final.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

####################################################################################
########################### profile ethnicity enrichment ###########################
####################################################################################

PBMC.integrated.final <- FindNeighbors(object = PBMC.integrated.final, 
                                       reduction = "pca", dims = 1:17, 
                                       k.param = 500, verbose = TRUE,compute.SNN = FALSE)

graph_nnMat_final <- Matrix::Matrix(PBMC.integrated.final@meta.data[,"integrated_nn"], sparse = TRUE)
saveRDS(graph_nnMat_final,"graph_nnMat_final.rds")
graph_nnMat_final <- readRDS("graph_nnMat_final.rds")

umap.df.final <- read.table("./umap.df.final.txt", header = TRUE, sep = "\t")


comp <- c("Japanese", "Korean","Malay","Indian","Caucasian")

ref <- "Chinese"
ref_num <- nrow(umap.df.final[(umap.df.final$ethnicity==ref),,drop=FALSE])
# calculate total number of ref cells in each cell's 500 neighbours
graph_nnMat_final.ref <- graph_nnMat_final[,umap.df.final[(umap.df.final$ethnicity==ref),"barcode"]]
graph_nnMat_final.ref.sum <- Matrix::rowSums(graph_nnMat_final.ref)


for (i in seq_along(comp)){
  comp.i <- comp[i]
  
  
  # calculate number of ref and comp cells in each cell's 300 neighbours
  graph_nnMat_final.comp <- graph_nnMat_final[, umap.df.final[(umap.df.final$ethnicity==comp.i),"barcode"] ]
  graph_nnMat_final.comp.sum <- Matrix::rowSums(graph_nnMat_final.comp)
  
  comp.ref.df <- data.frame(comp = graph_nnMat_final.comp.sum,
                            ref = graph_nnMat_final.ref.sum)
  
  comp_FC <- (nrow(umap.df.final[(umap.df.final$ethnicity==comp.i),,drop=FALSE])) / (ref_num)
  
  comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
  comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)
  
  comp.ref.df[(comp.ref.df$comp == 0 & comp.ref.df$ref == 0),"log2Ratio"] <- 0
  
  comp.ref.df <- comp.ref.df[umap.df.final$barcode,]
  
  umap.df.final$log2Ratio <- comp.ref.df$log2Ratio
  umap.df.final[(umap.df.final$log2Ratio > 3),"log2Ratio"] <- 3
  umap.df.final[(umap.df.final$log2Ratio < -3),"log2Ratio"] <- -3
  
  ggplot()+
    geom_point(umap.df.final, 
               mapping = aes(x=UMAP_1, y=UMAP_2, 
                             color=log2Ratio), 
               size=0.2, stroke=0)+
    scale_color_gradientn(colors = c("darkorange","#e8ded3","#E8E8E8","#d5e0e8","darkblue" ),
                          values=c(1, .54, .5, .46, 0),
                          limits=c(-3,3))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
    ggsave(paste0("3.Comparative_enrichment_plot_of_",
                  comp.i,"_vs_",ref,"_discarded_one_JP_sample.png"),device = "png",
           width = 7, height = 6)
  
}


##### profile commercial controls in each country
umap.df.final <- read.table("./umap.df.final.txt", header = TRUE, sep = "\t")
umap.df.lonza <- umap.df.final[(umap.df.final$ethnicity == "Caucasian"), ]

umap.df.lonza$country <- factor(umap.df.lonza$country,
                                levels =  c("KR","SG","JP"))

ggplot()+
  geom_point(umap.df.final, mapping = aes(x=UMAP_1,y=UMAP_2), color="grey",size=.15,stroke=0)+
  geom_point(umap.df.lonza, 
             mapping = aes(x=UMAP_1,y=UMAP_2, color = country),size=.15,stroke=0)+
  theme_bw(10)+
  scale_color_manual(values = c("darkgreen","red","steelblue"))+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave(paste0("4.UMAP_highlighting_commercial_control_in_lonza_all_discarded_one_JP_sample.png"),
       device = "png", width = 6, height = 5)

unique_lonzaID <- unique(umap.df.lonza$DCP_ID)

for(i in seq_along(unique_lonzaID)){
  unique_lonzaID.i <- unique_lonzaID[i]
  
  ggplot()+
    geom_point(umap.df.final, mapping = aes(x=UMAP_1,y=UMAP_2), color="grey",size=.15,stroke=0)+
    geom_point(umap.df.lonza[(umap.df.lonza$DCP_ID==unique_lonzaID.i),], 
               mapping = aes(x=UMAP_1,y=UMAP_2, color = country),size=.3,stroke=0)+
    theme_bw(10)+
    scale_color_manual(values = c("darkgreen","red","steelblue"))+
    guides(colour = guide_legend(override.aes = list(size=5)))
  ggsave(paste0("4.UMAP_highlighting_commercial_control_in_lonza",unique_lonzaID.i,
                "_discarded_one_JP_sample.png"),
         device = "png", width = 6, height = 5)
}

country.c <- c("KR","SG","JP")
color.c <- c("darkgreen","red","steelblue")

for(i in seq_along(country.c)){
  country.i <- country.c[i]
  
  ggplot()+
    geom_point(umap.df.final, mapping = aes(x=UMAP_1,y=UMAP_2), 
               color="grey",size=.15,stroke=0)+
    geom_point(umap.df.lonza[(umap.df.lonza$country==country.i),], 
               mapping = aes(x=UMAP_1,y=UMAP_2),
               size=.2,stroke=0, color=color.c[i])+
    theme_bw(10)+
    guides(colour = guide_legend(override.aes = list(size=5)))+
    ggtitle(paste0("n = ", nrow(umap.df.lonza[(umap.df.lonza$country==country.i),])))
  ggsave(paste0("4.UMAP_highlighting_commercial_control_in_",country.i,
                "_discarded_one_JP_sample.png"),
         device = "png", width = 6, height = 5)
}






######### profile commercial control enrichments across different country (combine all lots)
graph_nnMat_final <- readRDS("graph_nnMat_final.rds")

umap.df.final <- read.table("./umap.df.final.txt", header = TRUE, sep = "\t")
umap.df.lonza <- umap.df.final[(umap.df.final$ethnicity == "Caucasian"), ]

######### using SG as reference 
comp <- c("JP", "KR")

ref <- "SG"
ref_num <- nrow(umap.df.lonza[(umap.df.lonza$country==ref),,drop=FALSE])
# calculate total number of ref cells in each cell's 500 neighbours
graph_nnMat_final.ref <- graph_nnMat_final[,umap.df.lonza[(umap.df.lonza$country==ref),"barcode"]]
graph_nnMat_final.ref.sum <- Matrix::rowSums(graph_nnMat_final.ref)


for (i in seq_along(comp)){
  comp.i <- comp[i]
  
  
  # calculate number of ref and comp cells in each cell's 300 neighbours
  graph_nnMat_final.comp <- graph_nnMat_final[, umap.df.lonza[(umap.df.lonza$country==comp.i),"barcode"] ]
  graph_nnMat_final.comp.sum <- Matrix::rowSums(graph_nnMat_final.comp)
  
  comp.ref.df <- data.frame(comp = graph_nnMat_final.comp.sum,
                            ref = graph_nnMat_final.ref.sum)
  
  comp_FC <- (nrow(umap.df.lonza[(umap.df.lonza$country==comp.i),,drop=FALSE])) / (ref_num)
  
  comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
  comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)
  
  comp.ref.df[(comp.ref.df$comp <= 5 & comp.ref.df$ref <= 5),"log2Ratio"] <- 0
  
  comp.ref.df <- comp.ref.df[umap.df.final$barcode,]
  
  umap.df.final$log2Ratio <- comp.ref.df$log2Ratio
  umap.df.final[(umap.df.final$log2Ratio > 3),"log2Ratio"] <- 3
  umap.df.final[(umap.df.final$log2Ratio < -3),"log2Ratio"] <- -3
  
  ggplot()+
    geom_point(umap.df.final, 
               mapping = aes(x=UMAP_1, y=UMAP_2, 
                             color=log2Ratio), 
               size=0.2, stroke=0)+
    scale_color_gradientn(colors = c("darkorange","#e8ded3","#E8E8E8","#d5e0e8","darkblue" ),
                          values=c(1, .54, .5, .46, 0),
                          limits=c(-3,3))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  ggsave(paste0("4.2.Comparative_enrichment_plot_in_commercial_controls_of_",
                comp.i,"_vs_",ref,"_discarded_one_JP_sample.png"),device = "png",
         width = 7, height = 6)
  
}

### profile lonza relative enriched ratio in each lot and each country
unique_lonza <- unique(umap.df.lonza$DCP_ID)
unique_country <- unique(umap.df.lonza$country)
for (i in seq_along(unique_lonza)){
  for (j in seq_along(unique_country)){
    unique_lonza.i <- unique_lonza[i]
    unique_country.j <- unique_country[j]
    umap.df.lonza.ij <- umap.df.lonza[which(umap.df.lonza$DCP_ID == unique_lonza.i &
                                              umap.df.lonza$country ==  unique_country.j),]
    if (nrow(umap.df.lonza.ij)==0){
      next
    }
    graph_nnMat.ij <- graph_nnMat_final[,umap.df.lonza.ij$barcode]
    graph_nnMat.ij.sum <- Matrix::rowSums(graph_nnMat.ij)
    graph_nnMat.ij.ratio <- graph_nnMat.ij.sum/500
    
    expected_ratio <- nrow(umap.df.lonza.ij)/nrow(umap.df.final)
    
    graph_nnMat.ij.ratio <- graph_nnMat.ij.ratio[umap.df.final$barcode]
    umap.df.final$ratio <- graph_nnMat.ij.ratio
    umap.df.final$adjusted <- log2(umap.df.final$ratio/expected_ratio)
    umap.df.final[(umap.df.final$adjusted == -Inf),"adjusted"] <- -4
    
    umap.df.final[(umap.df.final$adjusted > 4),"adjusted"] <- 4
    umap.df.final[(umap.df.final$adjusted < -4),"adjusted"] <- -4
    
    ggplot()+
      geom_point(umap.df.final, 
                 mapping = aes(x=UMAP_1, y=UMAP_2, 
                               color=adjusted), 
                 size=0.2, stroke=0)+
      scale_color_gradientn(colors = c("darkorange","#e8ded3","#E8E8E8","#d5e0e8","darkblue" ),
                            values=c(1, .58, .5, .42, 0),
                            limits=c(-4,4))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    ggsave(paste0("4.3.integrated_UMAP_with_",unique_lonza.i,"_in_",unique_country.j,".png"),
           device = "png",width = 7, height = 6)
  }
}



###################### cell composition in each ethnicity
# get celltype proportion for different country
umap.df.final <- read.table("./umap.df.final.txt", header = TRUE, sep = "\t")

comp <- c("Japanese", "Korean","Malay","Indian","Chinese")

for (i in seq_along(comp)){
  comp.i <- comp[i]
  umap.df.i <- umap.df.final[which(umap.df.final$ethnicity==comp.i),]
  umap.df.i.c <- as.data.frame(table(umap.df.i$cell))
  umap.df.i.c <- umap.df.i.c %>% mutate(Var1 = factor(Var1, levels = umap.df.i.c$Var1),
                                                label = paste0(Var1, " ", round(Freq / sum(Freq) * 100, 1), "%"))
  pdf(paste0("5.cell_type_composition_in_",comp.i,".pdf"))
  pie(x = umap.df.i.c$Freq, labels = umap.df.i.c$label, 
      col =c("blue","brown","darkgreen",
             "darkorange","red","pink","turquoise","steelblue"),
      main=paste0("n = ",nrow(umap.df.i)))
  dev.off()
}

## profile cell type proportion in lonza for different country
umap.df.lonza <- umap.df.final[(umap.df.final$ethnicity == "Caucasian"), ]
country.c <- c("JP","SG","KR")
for (i in seq_along(country.c)) {
  country.i <- country.c[i]
  
  umap.df.i <- umap.df.lonza[which(umap.df.lonza$country==country.i),]
  umap.df.i.c <- as.data.frame(table(umap.df.i$cell))
  umap.df.i.c <- umap.df.i.c %>% mutate(Var1 = factor(Var1, levels = umap.df.i.c$Var1),
                                        label = paste0(Var1, " ", round(Freq / sum(Freq) * 100, 1), "%"))
  pdf(paste0("5.2.cell_type_composition_in_",country.i,"_commercial_controls.pdf"))
  pie(x = umap.df.i.c$Freq, labels = umap.df.i.c$label, 
      col =c("blue","brown","darkgreen",
             "darkorange","red","pink","turquoise","steelblue"),
      main=paste0("n = ",nrow(umap.df.i)))
  dev.off()
}


##### test1: compare JP lonza in one library with other samples
umap.df <- read.table("./umap.df.final.txt", header = TRUE, sep = "\t")
umap.df.jp <- umap.df[(umap.df$country == "JP"),]

# example LONZA3038306

colors37 = c("#466791","#60BF37","#953ADA","#4FBE6C","#CE49D3","#A7B43D","#5A51DC",
             "#D49F36","#552095","#507F2D","#DB37AA","#84B67C","#A06FDA","#DF462A",
             "#5B83DB","#C76C2D","#4F49A3","#82702D","#DD6BBB","#334C22","#D83979",
             "#55BAAD","#DC4555","#62AAD3","#8C3025","#417D61","#862977","#BBA672",
             "#403367","#DA8A6D","#A79CD4","#71482C","#C689D0","#6B2940","#D593A7",
             "#895C8B","#BD5975")

umap.df.jp.sub <- umap.df.jp[(umap.df.jp$lib %in% c("JP_B7_L1","JP_B7_L2")),]
unique_sample <- unique(umap.df.jp.sub$DCP_ID)
for(i in seq_along(unique_sample)){
  umap.df.jp.sub.i <- umap.df.jp.sub[(umap.df.jp.sub$DCP_ID == unique_sample[i]),]
  
  umap.df.jp.sub.i$lib <- factor(umap.df.jp.sub.i$lib)
  ggplot()+
    geom_point(umap.df, mapping = aes(x=UMAP_1, y=UMAP_2),
               color="lightgrey", stroke=0, size=0.2)+
    geom_point(umap.df.jp.sub.i, mapping = aes(x=UMAP_1, y=UMAP_2, col=lib),
               stroke=0, size=0.7)+
    scale_color_manual(values = colors37[(i*2-1):(i*2)])+
    theme_bw(10)+
    guides(colour = guide_legend(override.aes = list(size=5)))+
    ggtitle(unique_sample[i])
  ggsave(paste0("6.test_integrated_UMAP_highlighting_JP_sample_",
                unique_sample[i],".png"), device = "png",
         width = 6, height = 5)
    
    
    
}



umap.df.sg.sub <- umap.df[(umap.df$lib %in% c("SG_B3_L1","SG_B3_L2")),]
unique_sample <- unique(umap.df.sg.sub$DCP_ID)
for(i in seq_along(unique_sample)){
  umap.df.sg.sub.i <- umap.df.sg.sub[(umap.df.sg.sub$DCP_ID == unique_sample[i]),]
  
  umap.df.sg.sub.i$lib <- factor(umap.df.sg.sub.i$lib)
  ggplot()+
    geom_point(umap.df, mapping = aes(x=UMAP_1, y=UMAP_2),
               color="lightgrey", stroke=0, size=0.2)+
    geom_point(umap.df.sg.sub.i, mapping = aes(x=UMAP_1, y=UMAP_2, col=lib),
               stroke=0, size=0.7)+
    scale_color_manual(values = colors37[(i*2-1):(i*2)])+
    theme_bw(10)+
    guides(colour = guide_legend(override.aes = list(size=5)))+
    ggtitle(unique_sample[i])
  ggsave(paste0("6.test_integrated_UMAP_highlighting_SG_sample_",
                unique_sample[i],".png"), device = "png",
         width = 6, height = 5)
  
  
  
}




umap.df.kr.sub <- umap.df[(umap.df$lib %in% c("KR_B5_L1","KR_B5_L2")),]
unique_sample <- unique(umap.df.kr.sub$DCP_ID)
for(i in seq_along(unique_sample)){
  umap.df.kr.sub.i <- umap.df.kr.sub[(umap.df.kr.sub$DCP_ID == unique_sample[i]),]
  
  umap.df.kr.sub.i$lib <- factor(umap.df.kr.sub.i$lib)
  ggplot()+
    geom_point(umap.df, mapping = aes(x=UMAP_1, y=UMAP_2),
               color="lightgrey", stroke=0, size=0.2)+
    geom_point(umap.df.kr.sub.i, mapping = aes(x=UMAP_1, y=UMAP_2, col=lib),
               stroke=0, size=0.7)+
    scale_color_manual(values = colors37[(i*2-1):(i*2)])+
    theme_bw(10)+
    guides(colour = guide_legend(override.aes = list(size=5)))+
    ggtitle(unique_sample[i])
  ggsave(paste0("6.test_integrated_UMAP_highlighting_KR_sample_",
                unique_sample[i],".png"), device = "png",
         width = 6, height = 5)
  
  
  
}

########### Major cell type abundance analysis ###########################
#### remove sample with cell number less than 250 #######################
## ethnicity
umap.df.final <- read.table("umap.df.final.txt", header = TRUE, sep = "\t")
umap.df.final$batch <- paste(umap.df.final$country,
                             unlist(lapply(umap.df.final$lib, 
                                           function(x) unlist(strsplit(x, split = "_"))[2] )), 
                             sep = "_")
umap.df.final$sampleBatch <- paste(umap.df.final$DCP_ID,
                                   umap.df.final$batch,
                                   sep = "-")

umap.df.sub <- umap.df.final[,c("cell","sampleBatch","ethnicity")]
umap.df.sub$count <- 1
abundance.df.final <- aggregate(count ~ ., umap.df.sub, FUN = sum)

abundance.df.final[(abundance.df.final$cell == "mDC"),"cell"] <- "cDC"


unique_sampleBatch <- unique(abundance.df.final$sampleBatch)
unique_celltype <- unique(abundance.df.final$cell)



abundance.df.final_new <- data.frame()
# calculate the proportion for each donor with clinical conditions
for (i in seq_along(unique_sampleBatch)){
  sampleBatch.i <- unique_sampleBatch[i]
  abundance.df.final.i <- abundance.df.final[which(abundance.df.final$sampleBatch==sampleBatch.i),]
  
  ### minimum cell number is 250
  if(sum(abundance.df.final.i$count) < 250){
    next
  }
  
  # add row without corresponding cell types
  unique_celltype.no <- unique_celltype[!(unique_celltype %in% abundance.df.final.i$cell)]
  if(length(unique_celltype.no) > 0){
    no.df <- data.frame(cell = unique_celltype.no,
                        sampleBatch = rep(sampleBatch.i, length(unique_celltype.no)),
                        ethnicity= rep(abundance.df.final.i$ethnicity[1], 
                                      length(unique_celltype.no)),
                        count = rep(0,length(unique_celltype.no)))
    abundance.df.final.i <- rbind(abundance.df.final.i, no.df)
  }
  
  
  
  sum.i <- sum(abundance.df.final.i$count)
  abundance.df.final.i$proportion <- 100*abundance.df.final.i$count/sum.i
  abundance.df.final_new <- rbind(abundance.df.final_new,abundance.df.final.i)
}



unique_celltype <- c("B cells","cDC","Monocytes","NK cells",
                     "pDC","Plasma B","Platelets","T cells")
#ylim.c <- c(30, 3, 50, 40, 2, 3,1)
## adjusted p-value
for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new1 <- abundance.df.final_new[which(abundance.df.final_new$cell==unique_celltype.i),]
  
  abundance.df.final_new1$ethnicity <- factor(abundance.df.final_new1$ethnicity,
                                             levels = c("Caucasian","Chinese","Indian",
                                                        "Japanese","Korean","Malay"))
  
  
  stat.test <- wilcox_test(data=abundance.df.final_new1,
                           proportion ~ ethnicity,
                           p.adjust.method = "BH")
  
  
  
  stat.test <- stat.test[(stat.test$p.adj<0.05),]
  
  stat.test <- stat.test %>% add_xy_position(x = "ethnicity")
  
  
  kruskal.i <- kruskal.test(proportion~ethnicity, data=abundance.df.final_new1)
  kruskal.i.res <- kruskal.i$p.value
  
  ggboxplot(abundance.df.final_new1, 
            x = "ethnicity", y = "proportion", fill = "ethnicity",
            outlier.size = 0, alpha=0.7,
            outlier.shape = NA) + 
    geom_jitter(size=1.5, stroke=0, width = .2 )+
    stat_pvalue_manual(stat.test, label = "p.adj")+
    ylab("proportion (%)")+
    xlab("")+
    scale_fill_manual(values = c("#E1BC00","#A5D700","#00997F",
                                 "#1313A9","#B20086","#F20000"))+
    ggtitle(paste0("Kruskal p=",kruskal.i.res))
  ggsave(paste0("7.differential_abundance_across_samples_in_",
                gsub(" ", "_", unique_celltype.i),"_for_ethnicity.pdf"),width = 6, height = 4)
  
  ggboxplot(abundance.df.final_new1, 
            x = "ethnicity", y = "proportion", fill = "ethnicity",
            outlier.size = 0, alpha=0.7,
            outlier.shape = NA) + 
    geom_jitter(size=1.5, stroke=0, width = .2 )+
    ylab("proportion (%)")+
    xlab("")+
    scale_fill_manual(values = c("#E1BC00","#A5D700","#00997F",
                                 "#1313A9","#B20086","#F20000"))+
    ggtitle(paste0("Kruskal p=",kruskal.i.res))
  ggsave(paste0("7.differential_abundance_across_samples_in_",
                gsub(" ", "_", unique_celltype.i),"_for_ethnicity3.pdf"),width = 6, height = 4)
  
}

## age
age.df <- unique(umap.df.final[,c("sampleBatch","age")])
age.df <- age.df[(age.df$sampleBatch %in% abundance.df.final_new$sampleBatch),]
abundance.df.age <- merge(x = abundance.df.final_new,
                          y = age.df,
                          by = "sampleBatch")
abundance.df.age <- abundance.df.age[(abundance.df.age$age != "N.A."),]
# remove lonza samples
abundance.df.age <- abundance.df.age[(abundance.df.age$ethnicity != "Caucasian"),]

abundance.df.age$age <- as.numeric(as.character(abundance.df.age$age))
unique_celltype <- c("B cells","cDC","Monocytes","NK cells",
                     "pDC","Plasma B","Platelets","T cells")

for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.age.i <- abundance.df.age[(abundance.df.age$cell == unique_celltype.i),]
  
  
  
  cor.res <- cor.test(abundance.df.age.i$age, abundance.df.age.i$proportion,
                      method = "pearson")
  ggplot()+
    geom_point(abundance.df.age.i, 
               mapping = aes(x=age, y=proportion))+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          panel.background = element_blank())+
    ggtitle(paste0("cell: ",unique_celltype.i,"\nPearson cor = ", cor.res$estimate,
                   "\np-value = ", cor.res$p.value))
  ggsave(paste0("8.correlation_between_age_and_cell_abundance_in_",
                gsub(" ", "_", unique_celltype.i),".pdf" ))
}

#### different gender
sex.df <- unique(umap.df.final[,c("sampleBatch","sex")])
sex.df <- sex.df[(sex.df$sampleBatch %in% abundance.df.final_new$sampleBatch),]
abundance.df.sex <- merge(x = abundance.df.final_new,
                          y = sex.df,
                          by = "sampleBatch")
abundance.df.sex <- abundance.df.sex[(abundance.df.sex$sex != "N.A."),]
abundance.df.sex <- abundance.df.sex[(abundance.df.sex$ethnicity != "Caucasian"),]
unique_celltype <- c("B cells","cDC","Monocytes","NK cells",
                     "pDC","Plasma B","Platelets","T cells")

for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.sex.i <- abundance.df.sex[(abundance.df.sex$cell == unique_celltype.i),]
  
  abundance.df.sex.i$ethnicity <- factor(abundance.df.sex.i$ethnicity,
                                              levels = c("Chinese","Indian","Malay",
                                                         "Japanese","Korean"))
  abundance.df.sex.i$sex <- factor(abundance.df.sex.i$sex,
                                   levels = c("Male","Female"))
  
  
  ggboxplot(abundance.df.sex.i, x = "ethnicity", y = "proportion",
            color = "sex", add = "jitter", add.params = list(size=.5))+
    stat_compare_means(aes(group = sex),method = "wilcox.test", label = "p.signif")+
    labs(caption = "ns: p > 0.05; *: p <= 0.05; \n**: p <= 0.01; ***: p <= 0.001; \n****: p <= 0.0001")
  ggsave(paste0("9.cell_abundance_between_sex_among_different_ethnicity_in_",
                  gsub(" ", "_", unique_celltype.i),".pdf"),height = 4)
}


############### multiple regression across ethnicity, age and sex for cell abundance analysis
#### remove sample with cell number less than 250 #######################
umap.df.final <- read.table("umap.df.final.txt", header = TRUE, sep = "\t")
umap.df.final$batch <- paste(umap.df.final$country,
                             unlist(lapply(umap.df.final$lib, 
                                           function(x) unlist(strsplit(x, split = "_"))[2] )), 
                             sep = "_")
umap.df.final$sampleBatch <- paste(umap.df.final$DCP_ID,
                                   umap.df.final$batch,
                                   sep = "-")

umap.df.sub <- umap.df.final[,c("cell","sampleBatch","ethnicity")]
umap.df.sub$count <- 1
abundance.df.final <- aggregate(count ~ ., umap.df.sub, FUN = sum)

abundance.df.final[(abundance.df.final$cell == "mDC"),"cell"] <- "cDC"


unique_sampleBatch <- unique(abundance.df.final$sampleBatch)
unique_celltype <- unique(abundance.df.final$cell)

abundance.df.final_new <- data.frame()
# calculate the proportion for each donor with clinical conditions
for (i in seq_along(unique_sampleBatch)){
  sampleBatch.i <- unique_sampleBatch[i]
  abundance.df.final.i <- abundance.df.final[which(abundance.df.final$sampleBatch==sampleBatch.i),]
  
  ### minimum cell number is 250
  if(sum(abundance.df.final.i$count) < 250){
    next
  }
  
  # add row without corresponding cell types
  unique_celltype.no <- unique_celltype[!(unique_celltype %in% abundance.df.final.i$cell)]
  if(length(unique_celltype.no) > 0){
    no.df <- data.frame(cell = unique_celltype.no,
                        sampleBatch = rep(sampleBatch.i, length(unique_celltype.no)),
                        ethnicity= rep(abundance.df.final.i$ethnicity[1], 
                                       length(unique_celltype.no)),
                        count = rep(0,length(unique_celltype.no)))
    abundance.df.final.i <- rbind(abundance.df.final.i, no.df)
  }
  
  
  
  sum.i <- sum(abundance.df.final.i$count)
  abundance.df.final.i$proportion <- 100*abundance.df.final.i$count/sum.i
  abundance.df.final_new <- rbind(abundance.df.final_new,abundance.df.final.i)
}

# remove lonza
abundance.df.final_new <- abundance.df.final_new[(abundance.df.final_new$ethnicity != "Caucasian"),]
# add age
age.df <- unique(umap.df.final[,c("sampleBatch","age")])
age.df <- age.df[(age.df$sampleBatch %in% abundance.df.final_new$sampleBatch),]
abundance.df.final_new2 <- merge(x = abundance.df.final_new,
                          y = age.df,
                          by = "sampleBatch")
abundance.df.final_new2 <- abundance.df.final_new2[(abundance.df.final_new2$age != "N.A."),]

# add gender
sex.df <- unique(umap.df.final[,c("sampleBatch","sex")])
sex.df <- sex.df[(sex.df$sampleBatch %in% abundance.df.final_new2$sampleBatch),]
abundance.df.final_new3 <- merge(x = abundance.df.final_new2,
                                 y = sex.df,
                                 by = "sampleBatch")
abundance.df.final_new3 <- abundance.df.final_new3[(abundance.df.final_new3$age != "N.A."),]
abundance.df.final_new3$sex_code <- 2
abundance.df.final_new3[(abundance.df.final_new3$sex == "Female"),"sex_code"] <- 1
abundance.df.final_new3$age <- as.numeric(as.character(abundance.df.final_new3$age))
abundance.df.final_new3$ethnicity_code <- 1
abundance.df.final_new3[(abundance.df.final_new3$ethnicity == "Japanese"), "ethnicity_code"] <- 2
abundance.df.final_new3[(abundance.df.final_new3$ethnicity == "Korean"), "ethnicity_code"] <- 3
abundance.df.final_new3[(abundance.df.final_new3$ethnicity == "Malay"), "ethnicity_code"] <- 4
abundance.df.final_new3[(abundance.df.final_new3$ethnicity == "Indian"), "ethnicity_code"] <- 5

unique_celltype <- c("B cells","cDC","Monocytes","NK cells",
                     "pDC","Plasma B","Platelets","T cells")
model.all <- data.frame()
for(i in seq_along(unique_celltype)) {
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new3.i <- abundance.df.final_new3[(abundance.df.final_new3$cell == unique_celltype.i),]
  
  model <- lm(proportion ~ age + sex_code + ethnicity_code,
              data = abundance.df.final_new3.i)
  
  x <- as.data.frame(summary(model)$coefficient)
  new_data <- data.frame(cell = rep(unique_celltype.i,nrow(x)-1),
                         type = c("age","sex","ethnicity"),
                         estimates = x[c("age","sex_code","ethnicity_code"),
                                       "Estimate"],
                         se = x[c("age","sex_code","ethnicity_code"),
                                "Std. Error"],
                         pvalue = x[c("age","sex_code","ethnicity_code"),
                                    "Pr(>|t|)"])
  model.all <- rbind(model.all, new_data)
}

model.all$p_sign <- "ns"
model.all[(model.all$pvalue <= 0.05 & model.all$pvalue > 0.01),"p_sign"] <- "*"
model.all[(model.all$pvalue <= 0.01 & model.all$pvalue > 0.001),"p_sign"] <- "**"
model.all[(model.all$pvalue <= 0.001),"p_sign"] <- "***"

model.all$type <- factor(model.all$type,
                         levels = c("age","sex","ethnicity"))

ggplot(model.all, aes(x=cell, y=estimates, colour=type)) + 
  geom_point(position=position_dodge(.5)) + 
  geom_errorbar(aes(ymin=estimates-se, ymax=estimates+se), width=.1, position=position_dodge(.5)) + 
  geom_text(aes(label = p_sign, y = estimates+se), 
            vjust = -1.5,position=position_dodge(.5)) +
  ylim(c(-5,5))+
  scale_color_manual(values = c("orange","royalblue","darkgreen")) +
  labs(caption = "ns: p > 0.05; *: p <= 0.05; \n**: p <= 0.01; ***: p <= 0.001
       Sex: Female - 1; Male - 2
       Ethnicity: Chinese - 1; Japanese - 2;
                  Korean - 3; Malay - 4; Indian - 5")+
  xlab("")+
  ylab("Coefficient estimate")+
  theme_bw()
ggsave("10.coefficient_estimates_of_multiple_linear_regression_in_cell_abundance_across_sex_age_and_ethnicity.pdf",
       width = 10, height = 5)

##################### CHECK: If you calculate the average (“pseudobulk”) NODG 
##################### of each sample, do AIDA B cells have more sample to 
##################### sample NODG variation than, say, T cells or monocytes?

qc_matrix1 <- read.table("../combined_RCA_QC/2.QC_metrics_for_part1.txt", header = TRUE, sep = "\t")
qc_matrix2 <- read.table("../combined_RCA_QC/2.QC_metrics_for_part2.txt", header = TRUE, sep = "\t")
qc_matrix3 <- read.table("../combined_RCA_QC/2.QC_metrics_for_part3.txt", header = TRUE, sep = "\t")
qc_matrix4 <- read.table("../combined_RCA_QC/2.QC_metrics_for_part4.txt", header = TRUE, sep = "\t")
qc_matrix.all <- rbind(qc_matrix1, qc_matrix2, qc_matrix3,qc_matrix4)

final.df <- read.table("../combined_seurat_integration/umap.df.final.txt", header = TRUE,
                       sep = "\t")
final.df <- final.df[(final.df$ethnicity != "Caucasian"),]

sample.table <- as.data.frame(table(final.df$DCP_ID))
sample.table <- sample.table[(sample.table$Freq > 250),]
unique_sample <- as.character(sample.table$Var1)
unique_cell <- unique(final.df$cell)


NODG_sample_cell <- data.frame()
for(i in seq(unique_sample)){
  unique_sample.i <- unique_sample[i]
  for(j in seq_along(unique_cell)){
    unique_cell.j <- unique_cell[j]
    
    barcode.ij <- final.df[(final.df$DCP_ID == unique_sample.i & final.df$cell == unique_cell.j),"barcode"]
    if (length(barcode.ij) > 0){
      qc_matrix.ij <- qc_matrix.all[(qc_matrix.all$barcode %in% barcode.ij),]
      mean.nodg <- mean(qc_matrix.ij$NODG)
      
      new.ij <- data.frame(sample = unique_sample.i,
                           cell = unique_cell.j,
                           NODG = mean.nodg)
      
      NODG_sample_cell <- rbind(NODG_sample_cell, new.ij)
    } 
    
  }
}

ggboxplot(NODG_sample_cell, x = "cell", y = "NODG",
          color = "cell", palette = "jco",
          add = "jitter", add.params = list(size=0.1))
ggsave("test.NODG_across_different_cell_types.pdf", width = 7, height = 5)




########### profile marker genes for major cell types
PBMC <- readRDS("PBMC.integrated.after.QC.and.removing.Japan.withdrawn.sample.withMetadata.rds")
umap.df <- FetchData(PBMC, vars = c("UMAP_1","UMAP_2", "cell.annotation"))
# marker genes
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B",
                  "CD27","CD38","TNFRSF17","ITGAM","FCER1A","CD1C","LILRA4","IL3RA")


umap.df[(umap.df$cell.annotation == "mDC"), "cell.annotation"] <- "cDC"
unique.cluster <- c("T cells","NK cells","B cells","Monocytes","cDC","pDC","Plasma B","Platelets")

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(PBMC@assays$RNA@data)) {
    marker.i <- PBMC@assays$RNA@data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- rownames(umap.df[which(umap.df$cell.annotation==cluster.j),])
      marker.i.j <- marker.i[,barcode.j]
      exp_per.i.j <- sum(marker.i.j>0)*100/length(marker.i.j)
      exp_mean.i.j <- mean(marker.i.j)
      exp_per.i <- c(exp_per.i,exp_per.i.j)
      exp_mean.i <- c(exp_mean.i,exp_mean.i.j)
    }
    names(exp_per.i) <- unique.cluster
    
    exp_mean.i<- as.vector(scale(exp_mean.i,center = TRUE, scale = TRUE))
    names(exp_mean.i) <- unique.cluster
    marker.df.tmp <- data.frame(gene=rep(i,length(exp_mean.i)),
                                cluster = seq(1,length(exp_mean.i),1),
                                exp_per=exp_per.i,
                                exp_mean = exp_mean.i)
    
    marker.df <- rbind(marker.df,marker.df.tmp)
  }
  
}

marker.df$exp_per <- marker.df$exp_per/20
colnames(marker.df) <- c("gene","cluster","expresson_per","scaled_expression")

pdf("test.Marker_genes_across_different_major_cell_types.pdf", width = 5, height = 5)
ggplot()+geom_point(data = marker.df, 
                    mapping = aes(x=cluster,y=gene,
                                  size=expresson_per, color=scaled_expression))+
  scale_colour_gradient(low = "white",high = "red")+
  scale_x_continuous(breaks=seq(1,length(unique.cluster),1),labels=unique.cluster)+
  scale_y_continuous(breaks=seq(1,length(marker_genes),1),labels=marker_genes)+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
dev.off()


