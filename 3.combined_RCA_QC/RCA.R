library(RCAv2)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rstatix)

################################################################################
########################### read singlets rds #################################
################### for SG 1,3-19, JP 1-5 and KR 2-12 ##########################
################################################################################
celltype.all <- c()
# SG
for (i in c(1,seq(3,19,1))){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("SG_HEL_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("SG_HEL_B00",i,"_L00",j,"_5GEX")
    }
    
    print(lib_name)
    PBMC_r <- readRDS(paste0("../all_singlet_RDS/",lib_name,
                             ".singlets.RDS"))
    raw.ij <- PBMC_r$raw.data
    nom.ij <- PBMC_r$data
    colnames(raw.ij) <- paste(colnames(raw.ij), "-SG_B", i,"_L", j, sep = "")
    colnames(nom.ij) <- paste(colnames(nom.ij), "-SG_B", i,"_L", j, sep = "")
    
    estimate.cell.ij <- unlist(PBMC_r$cell.Type.Estimate)
    names(estimate.cell.ij) <- paste(colnames(PBMC_r$projection.data), 
                                     "-SG_B", i,"_L", j, sep = "")
    if (i==1 & j==1){
      raw.data.all  <- raw.ij
      norm.data.all <- nom.ij
    } else {
      raw.data.all <- Seurat::RowMergeSparseMatrices(raw.data.all, raw.ij)
      norm.data.all <- Seurat::RowMergeSparseMatrices(norm.data.all, nom.ij)
    }
    
    
    celltype.all <- c(celltype.all, estimate.cell.ij)
  }
}

# JP
for (i in seq(1,5,1)){
  for (j in c(1,2)){
    # skip JP B001 L002
    if (i==1 & j==2){
      next
    }
    
    if (i >= 10){
      lib_name <- paste0("JP_RIK_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("JP_RIK_B00",i,"_L00",j,"_5GEX")
    }
    
    print(lib_name)
    PBMC_r <- readRDS(paste0("../all_singlet_RDS/",lib_name,
                             ".singlets.RDS"))
    raw.ij <- PBMC_r$raw.data
    nom.ij <- PBMC_r$data
    colnames(raw.ij) <- paste(colnames(raw.ij), "-JP_B", i,"_L", j, sep = "")
    colnames(nom.ij) <- paste(colnames(nom.ij), "-JP_B", i,"_L", j, sep = "")
    
    estimate.cell.ij <- unlist(PBMC_r$cell.Type.Estimate)
    names(estimate.cell.ij) <- paste(colnames(PBMC_r$projection.data), 
                                     "-JP_B", i,"_L", j, sep = "")
      
    raw.data.all <- Seurat::RowMergeSparseMatrices(raw.data.all, raw.ij)
    norm.data.all <- Seurat::RowMergeSparseMatrices(norm.data.all, nom.ij)
    
    celltype.all <- c(celltype.all, estimate.cell.ij)
  }
}

# KR
for (i in seq(2,10,1)){
  for (j in c(1,2)){
    
    if (i >= 10){
      lib_name <- paste0("KR_SGI_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("KR_SGI_B00",i,"_L00",j,"_5GEX")
    }
    
    print(lib_name)
    PBMC_r <- readRDS(paste0("../all_singlet_RDS/",lib_name,
                             ".singlets.RDS"))
    raw.ij <- PBMC_r$raw.data
    nom.ij <- PBMC_r$data
    colnames(raw.ij) <- paste(colnames(raw.ij), "-KR_B", i,"_L", j, sep = "")
    colnames(nom.ij) <- paste(colnames(nom.ij), "-KR_B", i,"_L", j, sep = "")
    
    estimate.cell.ij <- unlist(PBMC_r$cell.Type.Estimate)
    names(estimate.cell.ij) <- paste(colnames(PBMC_r$projection.data), 
                                     "-KR_B", i,"_L", j, sep = "")
    
    raw.data.all <- Seurat::RowMergeSparseMatrices(raw.data.all, raw.ij)
    norm.data.all <- Seurat::RowMergeSparseMatrices(norm.data.all, nom.ij)
    
    celltype.all <- c(celltype.all, estimate.cell.ij)
  }
}



#### the matrix is too huge, need to separate the process
# KR
celltype.kr <- c()
for (i in seq(11,12,1)){
  for (j in c(1,2)){
    
    if (i >= 10){
      lib_name <- paste0("KR_SGI_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("KR_SGI_B00",i,"_L00",j,"_5GEX")
    }
    
    print(lib_name)
    PBMC_r <- readRDS(paste0("../all_singlet_RDS/",lib_name,
                             ".singlets.RDS"))
    raw.ij <- PBMC_r$raw.data
    nom.ij <- PBMC_r$data
    colnames(raw.ij) <- paste(colnames(raw.ij), "-KR_B", i,"_L", j, sep = "")
    colnames(nom.ij) <- paste(colnames(nom.ij), "-KR_B", i,"_L", j, sep = "")
    
    estimate.cell.ij <- unlist(PBMC_r$cell.Type.Estimate)
    names(estimate.cell.ij) <- paste(colnames(PBMC_r$projection.data), 
                                     "-KR_B", i,"_L", j, sep = "")
    
    if (i==11 & j==1){
      raw.data.kr  <- raw.ij
      norm.data.kr <- nom.ij
    } else {
      raw.data.kr <- Seurat::RowMergeSparseMatrices(raw.data.kr, raw.ij)
      norm.data.kr <- Seurat::RowMergeSparseMatrices(norm.data.kr, nom.ij)
    }
    
    
    celltype.kr <- c(celltype.kr, estimate.cell.ij)
  }
}

################################################################################
########################### create RCA object #################################
################### for SG 1,3-19, JP 1-5 and KR 2-12 ##########################
# due to huge cell number, split data into three clusters and annotate separately#
################################################################################
## split into 3 parts
all.cell <- data.frame(cell = colnames(RCA.all$raw.data))
all.cell$country <- unlist(lapply(all.cell$cell, function(x) unlist(strsplit(x,split="-"))[2] ))
all.cell.1 <- all.cell[(all.cell$country %in% c("SG_B1_L1","SG_B1_L2",
  "SG_B3_L1","SG_B3_L2",
  "SG_B4_L1","SG_B4_L2",
  "SG_B5_L1","SG_B5_L2",
  "SG_B6_L1","SG_B6_L2",
  "SG_B7_L1","SG_B7_L2",

  "KR_B2_L1","KR_B2_L2",
  "KR_B3_L1","KR_B3_L2",
  "KR_B4_L1","KR_B4_L2",
  "KR_B5_L1","KR_B5_L2",
  
  "JP_B1_L1",
  "JP_B2_L1","JP_B2_L2")),]

all.cell.2 <- all.cell[(all.cell$country %in% c("SG_B8_L1","SG_B8_L2",
  "SG_B9_L1","SG_B9_L2",
  "SG_B10_L1","SG_B10_L2",
  "SG_B11_L1","SG_B11_L2",
  "SG_B12_L1","SG_B12_L2",
  "SG_B13_L1","SG_B13_L2",

  "KR_B6_L1","KR_B6_L2",
  "KR_B7_L1","KR_B7_L2",
  "KR_B8_L1","KR_B8_L2",
  "KR_B9_L1","KR_B9_L2",

  "JP_B3_L1","JP_B3_L2")),]

all.cell.3 <- all.cell[(all.cell$country %in% c("SG_B14_L1","SG_B14_L2",
  "SG_B15_L1","SG_B15_L2",
  "SG_B16_L1","SG_B16_L2",
  "SG_B17_L1","SG_B17_L2",
  "SG_B18_L1","SG_B18_L2",
  "SG_B19_L1","SG_B19_L2",

  "KR_B10_L1","KR_B10_L2",

  "JP_B4_L1","JP_B4_L2",
  "JP_B5_L1","JP_B5_L2")),]

raw.all <- RCA.all$raw.data
norm.all <- RCA.all$data
estimate.cells <- unlist(RCA.all$cell.Type.Estimate)

## part 1 create RCA object
raw.all.1 <- raw.all[,all.cell.1$cell]
norm.all.1 <- norm.all[,all.cell.1$cell]
estimate.cells.1 <- estimate.cells[all.cell.1$cell]
raw.all.1 <- raw.all.1[rownames(norm.all.1),]

RCA.all.1 <- createRCAObject(rawData = raw.all.1, 
                           normData = norm.all.1)
estimate.cells.1 <- estimate.cells.1[colnames(raw.all.1)]
RCA.all.1$cell.Type.Estimate <- as.list(estimate.cells.1)
saveRDS(RCA.all.1,"RCA.combined.part1.rds")

## part 2 create RCA object
raw.all.2 <- raw.all[,all.cell.2$cell]
norm.all.2 <- norm.all[,all.cell.2$cell]
estimate.cells.2 <- estimate.cells[all.cell.2$cell]
raw.all.2 <- raw.all.2[rownames(norm.all.2),]

RCA.all.2 <- createRCAObject(rawData = raw.all.2, 
                           normData = norm.all.2)
estimate.cells.2 <- estimate.cells.2[colnames(raw.all.2)]
RCA.all.2$cell.Type.Estimate <- as.list(estimate.cells.2)
saveRDS(RCA.all.2,"RCA.combined.part2.rds")

## part 3 create RCA object
raw.all.3 <- raw.all[,all.cell.3$cell]
norm.all.3 <- norm.all[,all.cell.3$cell]
estimate.cells.3 <- estimate.cells[all.cell.3$cell]
raw.data.3.kr <- Seurat::RowMergeSparseMatrices(raw.all.3, raw.data.kr)
norm.data.3.kr <- Seurat::RowMergeSparseMatrices(norm.all.3, norm.data.kr)
estimate.cells.3.kr <- c(estimate.cells.3,celltype.kr)

raw.data.3.kr <- raw.data.3.kr[rownames(norm.data.3.kr),]

RCA.all.3 <- createRCAObject(rawData = raw.data.3.kr, 
                           normData = norm.data.3.kr)
estimate.cells.3.kr <- estimate.cells.3.kr[colnames(raw.data.3.kr)]
RCA.all.3$cell.Type.Estimate <- as.list(estimate.cells.3.kr)
saveRDS(RCA.all.3,"RCA.combined.part3.rds")



#RCA.all <- createRCAObject(rawData = raw.data.all, 
#                           normData = norm.data.all)
#celltype.all <- celltype.all[colnames(raw.data.all)]
#RCA.all$cell.Type.Estimate <- as.list(celltype.all)

#saveRDS(RCA.all,"RCA.combined.all.rds")

############ project to all immune panels ################


########################################################################################## 
############################## part 1 ##############################
########################################################################################## 
# remove gene "^MT-|^ERCC|^RPS|^RPL|^HSP"
gene.row <- rownames(RCA.all.1$raw.data)
gene.row.good <- gene.row[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                               x = gene.row, invert = TRUE)]

RCA.all.1$data <- RCA.all.1$data[gene.row.good,]


# get projection results against global panel
RCA.all.1 <- dataProject(RCA.all.1, method = "GlobalPanel",
                       corMeth = "pearson", scale = TRUE)
global.proj <- as.data.frame(as.matrix(RCA.all.1$projection.data))
global.proj.immune <- read.table("./rownames_of_glocal_projection_immune_cells.txt", 
                                 header = FALSE)
global.proj <- global.proj[global.proj.immune$V1,]

# get projection results against other two panels
RCA.all.1 <- dataProjectMultiPanel(RCA.all.1,method = list("NovershternPanel", 
                                                       "MonacoPanel", "CITESeqPanel"),
                                 scale = TRUE,
                                 corMeth = "pearson")
two.proj <- as.data.frame(as.matrix(RCA.all.1$projection.data))

# combine these panels
proj.all <- rbind(global.proj,two.proj)
proj.all <- as.matrix(proj.all)
proj.all <- as(proj.all, "dgCMatrix")

# Assign projection result to RCA object
RCA.all.1$projection.data <- proj.all

png("RCA.all.1.elbowplot.png")
RCAv2:::elbowPlot(RCA.all.1, nPCs = 50, approx = TRUE)
dev.off()
######################## dataSClust ########################
tempS<-Seurat::CreateSeuratObject(RCA.all.1$raw.data)
#From Seurat RunPCA
nPCs <- 13
res <- 2
npcs <- min(nPCs, nrow(RCA.all.1$projection.data) - 1)
pca.results <- irlba::irlba(A =t(RCA.all.1$projection.data), nv = npcs)
feature.loadings <- pca.results$v
sdev <- pca.results$d/sqrt(max(1, ncol(RCA.all.1$projection.data) - 1))
projection <- pca.results$u %*% diag(pca.results$d)

rownames(x = feature.loadings) <- rownames(x = RCA.all.1$projection.data)
colnames(x = feature.loadings) <- paste0("PC_", 1:npcs)
rownames(x = projection) <- colnames(x = RCA.all.1$projection.data)
colnames(x = projection) <- colnames(x = feature.loadings)
total.variance <- sum(apply(X=RCA.all.1$projection.data,MARGIN=2,FUN=var))
tempS@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = projection,
  loadings = feature.loadings,
  assay = "RNA",
  stdev = sdev,
  key = "PC_",
  misc = list(total.variance = total.variance))
tempS<-Seurat::FindNeighbors(object = tempS)
tempS<-Seurat::FindClusters(tempS,resolution = res)
# Convert labels to colours for each tree cut
if (length(unique(tempS$seurat_clusters))<41){
  dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
} else {
  if (require(randomcolorR) & require(plotrix)){
    clusterColors<-randomcoloR::distinctColorPalette(length(unique(tempS$seurat_clusters)))
    clusterColors<-sapply(clusterColors,plotrix::color.id)
    clusterColors<-sapply(clusterColors,function(x){return(x[1])})
    names(clusterColors)<-unique(clusteringResult$cluster)
    dynamicColorsList<-list(Colors=clusterColors[as.character(clusteringResult$cluster)])
  } else{
    dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
  }
}
names(dynamicColorsList)<-c("Clusters")
# Assign clustering result to RCA object
RCA.all.1$clustering.out <- list(
  "cellTree" = tempS$seurat_clusters,
  "dynamicColorsList" = dynamicColorsList
)

######################## dataSClust ########################

saveRDS(RCA.all.1, "RCA.combined.part1.rds")

### annotate each cluster
color.cluster <- data.frame(cell=unlist(RCA.all.1$cell.Type.Estimate),
                            color=RCA.all.1$clustering.out$dynamicColorsList[[1]])
color.cluster$count <- 1
color.cluster <- aggregate(count ~ ., color.cluster, FUN = sum)
cluster.final <- color.cluster %>% group_by(color) %>% top_n(n = 5, wt = count)
write.table(cluster.final,"cluster.final.part1.txt", sep="\t",
  col.names=FALSE, row.names=FALSE, quote=FALSE)

umap.df <- data.frame(barcode = colnames(RCA.all.1$raw.data),
                      cluster = RCA.all.1$clustering.out$dynamicColorsList[[1]])


umap.df$cell <- "unknown"
umap.df[which(umap.df$cluster=="bisque4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="black"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="brown4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="darkmagenta"),"cell"] <- "T cells" #check
umap.df[which(umap.df$cluster=="darkolivegreen"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkorange"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkorange2"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkslateblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="floralwhite"),"cell"] <- "T cells"

umap.df[which(umap.df$cluster=="green"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="ivory"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="lightcyan1"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="lightsteelblue1"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "Monocytes"

umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="mediumpurple3"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="orange"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="orangered4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="paleturquoise"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="plum1"),"cell"] <- "pDC"
umap.df[which(umap.df$cluster=="plum2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="red"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="saddlebrown"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="sienna3"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="skyblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="skyblue3"),"cell"] <- "mDC"
umap.df[which(umap.df$cluster=="steelblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="thistle1"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="thistle2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="violet"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="white"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="yellowgreen"),"cell"] <- "T cells"

# marker genes
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CXCR2","CD27","CD38","TNFRSF17","ITGAM","LILRA4","FCER1A","IL3RA","CD1C",
                  "CD74","HLA-DQB1","HLA-DRA")
# CD1C mDC

unique.cluster <- unique(umap.df$cluster)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(RCA.all.1$data)) {
    marker.i <- RCA.all.1$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cluster==cluster.j),"barcode"]
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

pdf("1.Marker_genes_across_different_clusters_part1.pdf", width = 12, height = 6)
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



umap.df$type <- "part1"
write.table(umap.df,"1.cell_annotation_for_part1.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

### get NODG and pMito
# plot QC
rawdata <- RCA.all.1$raw.data
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)

pMitoVec <- pMitoVec[names(nGeneVec)]
qc_matrix <- data.frame(barcode = names(pMitoVec),
                        NODG = nGeneVec,
                        pMt = pMitoVec,
                        type = "part1")
write.table(qc_matrix,"2.QC_metrics_for_part1.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)



########################################################################################## 
############################## part 2 ##############################
########################################################################################## 
# remove gene "^MT-|^ERCC|^RPS|^RPL|^HSP"
gene.row <- rownames(RCA.all.2$raw.data)
gene.row.good <- gene.row[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                               x = gene.row, invert = TRUE)]

RCA.all.2$data <- RCA.all.2$data[gene.row.good,]


# get projection results against global panel
RCA.all.2 <- dataProject(RCA.all.2, method = "GlobalPanel",
                       corMeth = "pearson", scale = TRUE)
global.proj <- as.data.frame(as.matrix(RCA.all.2$projection.data))
global.proj.immune <- read.table("./rownames_of_glocal_projection_immune_cells.txt", 
                                 header = FALSE)
global.proj <- global.proj[global.proj.immune$V1,]

# get projection results against other two panels
RCA.all.2 <- dataProjectMultiPanel(RCA.all.2,method = list("NovershternPanel", 
                                                       "MonacoPanel", "CITESeqPanel"),
                                 scale = TRUE,
                                 corMeth = "pearson")
two.proj <- as.data.frame(as.matrix(RCA.all.2$projection.data))

# combine these panels
proj.all <- rbind(global.proj,two.proj)
proj.all <- as.matrix(proj.all)
proj.all <- as(proj.all, "dgCMatrix")

# Assign projection result to RCA object
RCA.all.2$projection.data <- proj.all

RCAv2:::elbowPlot(RCA.all.2, nPCs = 50, approx = TRUE)
######################## dataSClust ########################
tempS<-Seurat::CreateSeuratObject(RCA.all.2$raw.data)
#From Seurat RunPCA
nPCs <- 13
res <- 2
npcs <- min(nPCs, nrow(RCA.all.2$projection.data) - 1)
pca.results <- irlba::irlba(A =t(RCA.all.2$projection.data), nv = npcs)
feature.loadings <- pca.results$v
sdev <- pca.results$d/sqrt(max(1, ncol(RCA.all.2$projection.data) - 1))
projection <- pca.results$u %*% diag(pca.results$d)

rownames(x = feature.loadings) <- rownames(x = RCA.all.2$projection.data)
colnames(x = feature.loadings) <- paste0("PC_", 1:npcs)
rownames(x = projection) <- colnames(x = RCA.all.2$projection.data)
colnames(x = projection) <- colnames(x = feature.loadings)
total.variance <- sum(apply(X=RCA.all.2$projection.data,MARGIN=2,FUN=var))
tempS@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = projection,
  loadings = feature.loadings,
  assay = "RNA",
  stdev = sdev,
  key = "PC_",
  misc = list(total.variance = total.variance))
tempS<-Seurat::FindNeighbors(object = tempS)
tempS<-Seurat::FindClusters(tempS,resolution = res)
# Convert labels to colours for each tree cut
if (length(unique(tempS$seurat_clusters))<41){
  dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
} else {
  if (require(randomcolorR) & require(plotrix)){
    clusterColors<-randomcoloR::distinctColorPalette(length(unique(tempS$seurat_clusters)))
    clusterColors<-sapply(clusterColors,plotrix::color.id)
    clusterColors<-sapply(clusterColors,function(x){return(x[1])})
    names(clusterColors)<-unique(clusteringResult$cluster)
    dynamicColorsList<-list(Colors=clusterColors[as.character(clusteringResult$cluster)])
  } else{
    dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
  }
}
names(dynamicColorsList)<-c("Clusters")
# Assign clustering result to RCA object
RCA.all.2$clustering.out <- list(
  "cellTree" = tempS$seurat_clusters,
  "dynamicColorsList" = dynamicColorsList
)

######################## dataSClust ########################
saveRDS(RCA.all.2, "RCA.combined.part2.rds")

### annotate each cluster
color.cluster <- data.frame(cell=unlist(RCA.all.2$cell.Type.Estimate),
                            color=RCA.all.2$clustering.out$dynamicColorsList[[1]])
color.cluster$count <- 1
color.cluster <- aggregate(count ~ ., color.cluster, FUN = sum)
cluster.final <- color.cluster %>% group_by(color) %>% top_n(n = 5, wt = count)
write.table(cluster.final,"cluster.final.part2.txt", sep="\t",
  col.names=FALSE, row.names=FALSE, quote=FALSE)

umap.df <- data.frame(barcode = colnames(RCA.all.2$raw.data),
                      cluster = RCA.all.2$clustering.out$dynamicColorsList[[1]])


umap.df$cell <- "unknown"
umap.df[which(umap.df$cluster=="bisque4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="black"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkmagenta"),"cell"] <- "mDC" #check
umap.df[which(umap.df$cluster=="darkolivegreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkorange"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="darkorange2"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="darkslateblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="floralwhite"),"cell"] <- "pDC"

umap.df[which(umap.df$cluster=="green"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="ivory"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="lightcyan1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightsteelblue1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "B cells"

umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="mediumpurple3"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="orange"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="orangered4"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="paleturquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="plum1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="plum2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="red"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="saddlebrown"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="sienna3"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="skyblue"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="skyblue3"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="steelblue"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="thistle1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="thistle2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="violet"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="white"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="yellowgreen"),"cell"] <- "Monocytes"



# marker genes
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CXCR2","CD27","CD38","TNFRSF17","ITGAM","LILRA4","FCER1A","IL3RA","CD1C",
                  "CD74","HLA-DQB1","HLA-DRA")
# CD1C mDC

unique.cluster <- unique(umap.df$cluster)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(RCA.all.2$data)) {
    marker.i <- RCA.all.2$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cluster==cluster.j),"barcode"]
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

pdf("1.Marker_genes_across_different_clusters_part2.pdf", width = 12, height = 6)
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


umap.df$type <- "part2"
write.table(umap.df,"1.cell_annotation_for_part2.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

### get NODG and pMito
# plot QC
rawdata <- RCA.all.2$raw.data
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)

pMitoVec <- pMitoVec[names(nGeneVec)]
qc_matrix <- data.frame(barcode = names(pMitoVec),
                        NODG = nGeneVec,
                        pMt = pMitoVec,
                        type = "part2")
write.table(qc_matrix,"2.QC_metrics_for_part2.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)




########################################################################################## 
############################## part 3 ##############################
########################################################################################## 
# remove gene "^MT-|^ERCC|^RPS|^RPL|^HSP"
gene.row <- rownames(RCA.all.3$raw.data)
gene.row.good <- gene.row[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                               x = gene.row, invert = TRUE)]

RCA.all.3$data <- RCA.all.3$data[gene.row.good,]


# get projection results against global panel
RCA.all.3 <- dataProject(RCA.all.3, method = "GlobalPanel",
                       corMeth = "pearson", scale = TRUE)
global.proj <- as.data.frame(as.matrix(RCA.all.3$projection.data))
global.proj.immune <- read.table("./rownames_of_glocal_projection_immune_cells.txt", 
                                 header = FALSE)
global.proj <- global.proj[global.proj.immune$V1,]

# get projection results against other two panels
RCA.all.3 <- dataProjectMultiPanel(RCA.all.3,method = list("NovershternPanel", 
                                                       "MonacoPanel", "CITESeqPanel"),
                                 scale = TRUE,
                                 corMeth = "pearson")
two.proj <- as.data.frame(as.matrix(RCA.all.3$projection.data))

# combine these panels
proj.all <- rbind(global.proj,two.proj)
proj.all <- as.matrix(proj.all)
proj.all <- as(proj.all, "dgCMatrix")

# Assign projection result to RCA object
RCA.all.3$projection.data <- proj.all

png("RCA.all.3.elbowplot.png")
RCAv2:::elbowPlot(RCA.all.3, nPCs = 50, approx = TRUE)
dev.off()
######################## dataSClust ########################
tempS<-Seurat::CreateSeuratObject(RCA.all.3$raw.data)
#From Seurat RunPCA
nPCs <- 13
res <- 2
npcs <- min(nPCs, nrow(RCA.all.3$projection.data) - 1)
pca.results <- irlba::irlba(A =t(RCA.all.3$projection.data), nv = npcs)
feature.loadings <- pca.results$v
sdev <- pca.results$d/sqrt(max(1, ncol(RCA.all.3$projection.data) - 1))
projection <- pca.results$u %*% diag(pca.results$d)

rownames(x = feature.loadings) <- rownames(x = RCA.all.3$projection.data)
colnames(x = feature.loadings) <- paste0("PC_", 1:npcs)
rownames(x = projection) <- colnames(x = RCA.all.3$projection.data)
colnames(x = projection) <- colnames(x = feature.loadings)
total.variance <- sum(apply(X=RCA.all.3$projection.data,MARGIN=2,FUN=var))
tempS@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = projection,
  loadings = feature.loadings,
  assay = "RNA",
  stdev = sdev,
  key = "PC_",
  misc = list(total.variance = total.variance))
tempS<-Seurat::FindNeighbors(object = tempS)
tempS<-Seurat::FindClusters(tempS,resolution = res)
# Convert labels to colours for each tree cut
if (length(unique(tempS$seurat_clusters))<41){
  dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
} else {
  if (require(randomcolorR) & require(plotrix)){
    clusterColors<-randomcoloR::distinctColorPalette(length(unique(tempS$seurat_clusters)))
    clusterColors<-sapply(clusterColors,plotrix::color.id)
    clusterColors<-sapply(clusterColors,function(x){return(x[1])})
    names(clusterColors)<-unique(clusteringResult$cluster)
    dynamicColorsList<-list(Colors=clusterColors[as.character(clusteringResult$cluster)])
  } else{
    dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
  }
}
names(dynamicColorsList)<-c("Clusters")
# Assign clustering result to RCA object
RCA.all.3$clustering.out <- list(
  "cellTree" = tempS$seurat_clusters,
  "dynamicColorsList" = dynamicColorsList
)

######################## dataSClust ########################

saveRDS(RCA.all.3, "RCA.combined.part3.rds")

### annotate each cluster
color.cluster <- data.frame(cell=unlist(RCA.all.3$cell.Type.Estimate),
                            color=RCA.all.3$clustering.out$dynamicColorsList[[1]])
color.cluster$count <- 1
color.cluster <- aggregate(count ~ ., color.cluster, FUN = sum)
cluster.final <- color.cluster %>% group_by(color) %>% top_n(n = 5, wt = count)
write.table(cluster.final,"cluster.final.part3.txt", sep="\t",
  col.names=FALSE, row.names=FALSE, quote=FALSE)

umap.df <- data.frame(barcode = colnames(RCA.all.3$raw.data),
                      cluster = RCA.all.3$clustering.out$dynamicColorsList[[1]])


umap.df$cell <- "unknown"

umap.df[which(umap.df$cluster=="bisque4"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="black"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown4"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkmagenta"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="darkolivegreen"),"cell"] <- "mDC"
umap.df[which(umap.df$cluster=="darkorange"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkorange2"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkslateblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="floralwhite"),"cell"] <- "NK cells"

umap.df[which(umap.df$cluster=="green"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="ivory"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightcyan1"),"cell"] <- "pDC"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="lightsteelblue1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "T cells"

umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="mediumpurple3"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="orange"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="orangered4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="paleturquoise"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="plum1"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="plum2"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="red"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="saddlebrown"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="sienna3"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="skyblue"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="skyblue3"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="steelblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="thistle1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="thistle2"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="violet"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="white"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="yellowgreen"),"cell"] <- "Monocytes"


# marker genes
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CXCR2","CD27","CD38","TNFRSF17","ITGAM","LILRA4","FCER1A","IL3RA","CD1C",
                  "CD74","HLA-DQB1","HLA-DRA")
# CD1C mDC

unique.cluster <- unique(umap.df$cluster)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(RCA.all.3$data)) {
    marker.i <- RCA.all.3$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cluster==cluster.j),"barcode"]
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

pdf("1.Marker_genes_across_different_clusters_part3.pdf", width = 12, height = 6)
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


umap.df$type <- "part3"
write.table(umap.df,"1.cell_annotation_for_part3.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

### get NODG and pMito
# plot QC
rawdata <- RCA.all.3$raw.data
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)

pMitoVec <- pMitoVec[names(nGeneVec)]
qc_matrix <- data.frame(barcode = names(pMitoVec),
                        NODG = nGeneVec,
                        pMt = pMitoVec,
                        type = "part3")
write.table(qc_matrix,"2.QC_metrics_for_part3.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)


################################################################################
######################### below for SG 20-21, JP 6-10 #########################
############################### as part 4 ######################################
################################################################################


################################################################################
########################### read singlets rds #################################
#########################  for SG 20-21, JP 6-10 ##############################
############################### as part 4 ######################################
################################################################################

celltype.all <- c()
# SG
for (i in seq(20,21,1)){
  for (j in c(1,2)){
    if (i >= 10){
      lib_name <- paste0("SG_HEL_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("SG_HEL_B00",i,"_L00",j,"_5GEX")
    }
    
    print(lib_name)
    PBMC_r <- readRDS(paste0("../all_singlet_RDS/",lib_name,
                             ".singlets.RDS"))
    raw.ij <- PBMC_r$raw.data
    nom.ij <- PBMC_r$data
    colnames(raw.ij) <- paste(colnames(raw.ij), "-SG_B", i,"_L", j, sep = "")
    colnames(nom.ij) <- paste(colnames(nom.ij), "-SG_B", i,"_L", j, sep = "")
    
    estimate.cell.ij <- unlist(PBMC_r$cell.Type.Estimate)
    names(estimate.cell.ij) <- paste(colnames(PBMC_r$projection.data), 
                                     "-SG_B", i,"_L", j, sep = "")
    if (i==20 & j==1){
      raw.data.all  <- raw.ij
      norm.data.all <- nom.ij
    } else {
      raw.data.all <- Seurat::RowMergeSparseMatrices(raw.data.all, raw.ij)
      norm.data.all <- Seurat::RowMergeSparseMatrices(norm.data.all, nom.ij)
    }
    
    
    celltype.all <- c(celltype.all, estimate.cell.ij)
  }
}

# JP
for (i in seq(6,10,1)){
  for (j in c(1,2)){
    # skip JP B001 L002
    if (i==1 & j==2){
      next
    }
    
    if (i >= 10){
      lib_name <- paste0("JP_RIK_B0",i,"_L00",j,"_5GEX")
    } else {
      lib_name <- paste0("JP_RIK_B00",i,"_L00",j,"_5GEX")
    }
    
    print(lib_name)
    PBMC_r <- readRDS(paste0("../all_singlet_RDS/",lib_name,
                             ".singlets.RDS"))
    raw.ij <- PBMC_r$raw.data
    nom.ij <- PBMC_r$data
    colnames(raw.ij) <- paste(colnames(raw.ij), "-JP_B", i,"_L", j, sep = "")
    colnames(nom.ij) <- paste(colnames(nom.ij), "-JP_B", i,"_L", j, sep = "")
    
    estimate.cell.ij <- unlist(PBMC_r$cell.Type.Estimate)
    names(estimate.cell.ij) <- paste(colnames(PBMC_r$projection.data), 
                                     "-JP_B", i,"_L", j, sep = "")
    
    raw.data.all <- Seurat::RowMergeSparseMatrices(raw.data.all, raw.ij)
    norm.data.all <- Seurat::RowMergeSparseMatrices(norm.data.all, nom.ij)
    
    celltype.all <- c(celltype.all, estimate.cell.ij)
  }
}

raw.data.all <- raw.data.all[rownames(norm.data.all),]

RCA.all.4 <- createRCAObject(rawData = raw.data.all, 
                             normData = norm.data.all)
celltype.all <- celltype.all[colnames(raw.data.all)]
RCA.all.4$cell.Type.Estimate <- as.list(celltype.all)
saveRDS(RCA.all.4,"RCA.combined.part4.rds")

# remove gene "^MT-|^ERCC|^RPS|^RPL|^HSP"
gene.row <- rownames(RCA.all.4$raw.data)
gene.row.good <- gene.row[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                               x = gene.row, invert = TRUE)]

RCA.all.4$data <- RCA.all.4$data[gene.row.good,]


# get projection results against global panel
RCA.all.4 <- dataProject(RCA.all.4, method = "GlobalPanel",
                         corMeth = "pearson", scale = TRUE)
global.proj <- as.data.frame(as.matrix(RCA.all.4$projection.data))
global.proj.immune <- read.table("../rownames_of_glocal_projection_immune_cells.txt", 
                                 header = FALSE)
global.proj <- global.proj[global.proj.immune$V1,]

# get projection results against other two panels
RCA.all.4 <- dataProjectMultiPanel(RCA.all.4,method = list("NovershternPanel", 
                                                           "MonacoPanel", "CITESeqPanel"),
                                   scale = TRUE,
                                   corMeth = "pearson")
two.proj <- as.data.frame(as.matrix(RCA.all.4$projection.data))

# combine these panels
proj.all <- rbind(global.proj,two.proj)
proj.all <- as.matrix(proj.all)
proj.all <- as(proj.all, "dgCMatrix")

# Assign projection result to RCA object
RCA.all.4$projection.data <- proj.all

png("RCA.all.4.elbowplot.png")
RCAv2:::elbowPlot(RCA.all.4, nPCs = 50, approx = TRUE)
dev.off()
######################## dataSClust ########################
tempS<-Seurat::CreateSeuratObject(RCA.all.4$raw.data)
#From Seurat RunPCA
nPCs <- 13
res <- 2
npcs <- min(nPCs, nrow(RCA.all.4$projection.data) - 1)
pca.results <- irlba::irlba(A =t(RCA.all.4$projection.data), nv = npcs)
feature.loadings <- pca.results$v
sdev <- pca.results$d/sqrt(max(1, ncol(RCA.all.4$projection.data) - 1))
projection <- pca.results$u %*% diag(pca.results$d)

rownames(x = feature.loadings) <- rownames(x = RCA.all.4$projection.data)
colnames(x = feature.loadings) <- paste0("PC_", 1:npcs)
rownames(x = projection) <- colnames(x = RCA.all.4$projection.data)
colnames(x = projection) <- colnames(x = feature.loadings)
total.variance <- sum(apply(X=RCA.all.4$projection.data,MARGIN=2,FUN=var))
tempS@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = projection,
  loadings = feature.loadings,
  assay = "RNA",
  stdev = sdev,
  key = "PC_",
  misc = list(total.variance = total.variance))
tempS<-Seurat::FindNeighbors(object = tempS)
tempS<-Seurat::FindClusters(tempS,resolution = res)
# Convert labels to colours for each tree cut
if (length(unique(tempS$seurat_clusters))<41){
  dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
} else {
  if (require(randomcolorR) & require(plotrix)){
    clusterColors<-randomcoloR::distinctColorPalette(length(unique(tempS$seurat_clusters)))
    clusterColors<-sapply(clusterColors,plotrix::color.id)
    clusterColors<-sapply(clusterColors,function(x){return(x[1])})
    names(clusterColors)<-unique(clusteringResult$cluster)
    dynamicColorsList<-list(Colors=clusterColors[as.character(clusteringResult$cluster)])
  } else{
    dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
  }
}
names(dynamicColorsList)<-c("Clusters")
# Assign clustering result to RCA object
RCA.all.4$clustering.out <- list(
  "cellTree" = tempS$seurat_clusters,
  "dynamicColorsList" = dynamicColorsList
)

######################## dataSClust ########################

saveRDS(RCA.all.4, "RCA.combined.part4.rds")

### annotate each cluster
color.cluster <- data.frame(cell=unlist(RCA.all.4$cell.Type.Estimate),
                            color=RCA.all.4$clustering.out$dynamicColorsList[[1]])
color.cluster$count <- 1
color.cluster <- aggregate(count ~ ., color.cluster, FUN = sum)
cluster.final <- color.cluster %>% group_by(color) %>% top_n(n = 5, wt = count)


umap.df <- data.frame(barcode = colnames(RCA.all.4$raw.data),
                      cluster = RCA.all.4$clustering.out$dynamicColorsList[[1]])


umap.df$cell <- "unknown"

umap.df[which(umap.df$cluster=="bisque4"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="black"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkmagenta"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="darkolivegreen"),"cell"] <- "mDC"
umap.df[which(umap.df$cluster=="darkorange"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkorange2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="darkslateblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="floralwhite"),"cell"] <- "B cells"

umap.df[which(umap.df$cluster=="green"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="ivory"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="lightcyan1"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="lightsteelblue1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "Monocytes"

umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="mediumpurple3"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="orange"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="orangered4"),"cell"] <- "pDC"
umap.df[which(umap.df$cluster=="paleturquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="plum1"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="plum2"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="red"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="saddlebrown"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="sienna3"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="skyblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="skyblue3"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="steelblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="thistle2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="violet"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="white"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="yellowgreen"),"cell"] <- "Monocytes"


# marker genes
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CXCR2","CD27","CD38","TNFRSF17","ITGAM","LILRA4","FCER1A","IL3RA","CD1C",
                  "CD74","HLA-DQB1","HLA-DRA")
# CD1C mDC

unique.cluster <- unique(umap.df$cluster)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(RCA.all.4$data)) {
    marker.i <- RCA.all.4$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cluster==cluster.j),"barcode"]
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

pdf("1.Marker_genes_across_different_clusters_part4.pdf", width = 12, height = 6)
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


umap.df$type <- "part4"
write.table(umap.df,"1.cell_annotation_for_part4.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

### get NODG and pMito
# plot QC
rawdata <- RCA.all.4$raw.data
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)

pMitoVec <- pMitoVec[names(nGeneVec)]
qc_matrix <- data.frame(barcode = names(pMitoVec),
                        NODG = nGeneVec,
                        pMt = pMitoVec,
                        type = "part4")
write.table(qc_matrix,"2.QC_metrics_for_part4.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

################################################################################
########################### combined all parts #################################
########################### QC on each cluster #################################
################################################################################
qc_matrix1 <- read.table("2.QC_metrics_for_part1.txt", header = TRUE, sep = "\t")
qc_matrix2 <- read.table("2.QC_metrics_for_part2.txt", header = TRUE, sep = "\t")
qc_matrix3 <- read.table("2.QC_metrics_for_part3.txt", header = TRUE, sep = "\t")
qc_matrix4 <- read.table("2.QC_metrics_for_part4.txt", header = TRUE, sep = "\t")
qc_matrix.all <- rbind(qc_matrix1, qc_matrix2, qc_matrix3,qc_matrix4)

cell_anno1 <- read.table("1.cell_annotation_for_part1.txt", header = TRUE, sep = "\t")
cell_anno2 <- read.table("1.cell_annotation_for_part2.txt", header = TRUE, sep = "\t")
cell_anno3 <- read.table("1.cell_annotation_for_part3.txt", header = TRUE, sep = "\t")
cell_anno4 <- read.table("1.cell_annotation_for_part4.txt", header = TRUE, sep = "\t")
cell_anno.all <- rbind(cell_anno1, cell_anno2, cell_anno3, cell_anno4)

##### QC on each cluster
unique_celltype <- unique(cell_anno.all$cell)
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- cell_anno.all[which(cell_anno.all$cell == cell.i), "barcode"]
  qc.i <- qc_matrix.all[(qc_matrix.all$barcode %in% barcode.i),]

  ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = 0.1, stroke=0) +
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  ggsave(paste0("2.QC_plot_of_each_cell_type_",gsub(" ","_",cell.i),".png"),
         width = 3, height = 3,device = "png")
}

# QC on each cluster

##### QC filtering
unique_celltype <- c("Monocytes",
                     "B cells",
                     "T cells",
                     "NK cells",
                     "pDC",
                     "mDC",
                     "Plasma B",
                     "Platelets")

nodg.low <- c(500,1000, 1100, 1100,1700,1500,1000,200)
nodg.up <- c(4000,3000,2900, 2900,4500,5500,6500,1200)
pmt.low <- c(0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001,0.0001)
pmt.up <- c(0.08, 0.08, 0.08, 0.08,0.08,0.08, 0.125, 0.125)
good.barcode <- c()

for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- cell_anno.all[which(cell_anno.all$cell == cell.i), "barcode"]
  qc.i <- qc_matrix.all[(qc_matrix.all$barcode %in% barcode.i),]
  rownames(qc.i) <- qc.i$barcode
  nodg.low.i <- nodg.low[i]
  nodg.up.i <- nodg.up[i]
  pmt.low.i <- pmt.low[i]
  pmt.up.i <- pmt.up[i]
  qc.good.i <- qc.i[which(qc.i$NODG >= nodg.low.i & qc.i$NODG <= nodg.up.i &
                            qc.i$pMt >= pmt.low.i & qc.i$pMt <= pmt.up.i),]
  good.barcode <- c(good.barcode, rownames(qc.good.i))
  qc.i$quality <- "bad"
  qc.i[rownames(qc.good.i),"quality"] <- "good"
  qc.i$quality <- factor(qc.i$quality, levels = c("good","bad"))
  ggplot2::ggplot() +
    ggplot2::geom_point(data = qc.i, ggplot2::aes(x = NODG, y = pMt,color=quality),size = .1,stroke=0) + 
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d(data = qc.i, mapping=aes(x = NODG, y = pMt)) +
    ggplot2::scale_color_manual(values = c("black","grey"))
  ggsave(paste0("3.QC_plot_of_each_cell_type_",gsub(" ","_",cell.i),"_highlighting_good.png"),
         width = 4, height = 3,device = "png")
}



# after QC

cell_anno.all.good <- cell_anno.all[(cell_anno.all$barcode %in% good.barcode),]
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- cell_anno.all.good[which(cell_anno.all.good$cell == cell.i), "barcode"]
  qc.i <- qc_matrix.all[(qc_matrix.all$barcode %in% barcode.i),]
  
  ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = 0.2, stroke=0) + 
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  ggsave(paste0("4.QC_plot_of_each_cell_type_",gsub(" ","_",cell.i),"_after_QC.png"),
         width = 3, height = 3,device = "png")
}

write.table(cell_anno.all.good, "3.cell_annotation_for_good_cells_after_QC.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#######################################################################################
##################################### add samples info and get the healthy only for sg
#######################################################################################

samples.df <- read.table("../all_singlet_RDS/all_singlets_mapped_to_their_samples.txt", 
                         header = TRUE, sep = "\t")
cell_anno.all.good <- read.table("./3.cell_annotation_for_good_cells_after_QC.txt", sep = "\t",
                                 header = TRUE)
samples.df.good <- samples.df[(samples.df$Barcode %in% cell_anno.all.good$barcode),]

cell_anno.all.good.sample <- merge(x = cell_anno.all.good, y = samples.df.good,
                                   by.x = "barcode", by.y = "Barcode")


no.healthy.list <- read.table("../all_singlet_RDS/nonhealthy_list_sg.txt", header = FALSE,
                              sep = "\t")
cell_anno.all.good.sample.healthy <- cell_anno.all.good.sample[!(cell_anno.all.good.sample$sample %in% no.healthy.list$V1),]

write.table(cell_anno.all.good.sample.healthy,
            "4.cell_annotation_for_good_cells_after_QC_healthy_only.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


### cell number per sample
sample.count <- as.data.frame(table(cell_anno.all.good.sample.healthy$sample))
sample.count <- sample.count[order(sample.count$Freq),]
sample.count$Var1 <- as.character(sample.count$Var1)
sample.count$Var1 <- factor(sample.count$Var1,
                            levels = sample.count$Var1)

ggplot(sample.count, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity")+
  geom_hline(yintercept = 300, col="red")


############################ QC after discarding one JP samples due to withdrawing
qc_matrix1 <- read.table("2.QC_metrics_for_part1.txt", header = TRUE, sep = "\t")
qc_matrix2 <- read.table("2.QC_metrics_for_part2.txt", header = TRUE, sep = "\t")
qc_matrix3 <- read.table("2.QC_metrics_for_part3.txt", header = TRUE, sep = "\t")
qc_matrix4 <- read.table("2.QC_metrics_for_part4.txt", header = TRUE, sep = "\t")
qc_matrix.all <- rbind(qc_matrix1, qc_matrix2, qc_matrix3,qc_matrix4)

final.df <- read.table("../combined_seurat_integration/umap.df.final.txt", header = TRUE,
                       sep = "\t")
qc_matrix.all.sub <- qc_matrix.all[(qc_matrix.all$barcode %in% final.df$barcode),]

unique_celltype <- c("Monocytes",
                     "B cells",
                     "T cells",
                     "NK cells",
                     "pDC",
                     "mDC",
                     "Plasma B",
                     "Platelets")

for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- final.df[which(final.df$cell == cell.i), "barcode"]
  qc.i <- qc_matrix.all.sub[(qc_matrix.all.sub$barcode %in% barcode.i),]
  
  ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = 0.2, stroke=0) + 
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  ggsave(paste0("5.QC_plot_of_each_cell_type_",gsub(" ","_",cell.i),"_after_QC_and_removing_1_JP_sample.png"),
         width = 3, height = 3,device = "png")
}

## QC metrcis between Chinese, Mala, and Indian
qc_matrix1 <- read.table("2.QC_metrics_for_part1.txt", header = TRUE, sep = "\t")
qc_matrix2 <- read.table("2.QC_metrics_for_part2.txt", header = TRUE, sep = "\t")
qc_matrix3 <- read.table("2.QC_metrics_for_part3.txt", header = TRUE, sep = "\t")
qc_matrix4 <- read.table("2.QC_metrics_for_part4.txt", header = TRUE, sep = "\t")
qc_matrix.all <- rbind(qc_matrix1, qc_matrix2, qc_matrix3,qc_matrix4)

final.df <- read.table("../combined_seurat_integration/umap.df.final.txt", header = TRUE,
                       sep = "\t")
qc_matrix.all.sub <- qc_matrix.all[(qc_matrix.all$barcode %in% final.df$barcode),]

eth.df <- final.df[,c("barcode","ethnicity")]

qc_matrix.eth.df <- merge(x = qc_matrix.all.sub,
                          y = eth.df, by="barcode")

qc_matrix.eth.df <- qc_matrix.eth.df[(qc_matrix.eth.df$ethnicity %in% c("Chinese","Malay","Indian")),]

qc_matrix.eth.df$ethnicity <- factor(qc_matrix.eth.df$ethnicity,
                                     levels = c("Chinese","Malay","Indian"))


ggplot(qc_matrix.eth.df, aes(x=ethnicity, y=NODG)) +
  geom_boxplot(aes(fill=ethnicity),outlier.shape = NA) + 
  scale_fill_manual(values = c("red","lightblue","lightyellow"))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  xlab("")+ylab("NODG")+ylim(c(0,4000))
ggsave("6.NODG_across_Chinese_Indian_and_Malay.pdf")
