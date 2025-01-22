###0. overview ####
#AIMS: 
## this code is for 1st (B, pDC + Myeloid, T + NK) and 2nd (CD4+T and nonC4+T cells) round of sub-clustering analysis 

#NOTES:
#1. should change the output directories (plotdir)
#2. change the directory for the seurat object input (seurat_obj)
#3. change the title (title_change), cell type name (cell_type), and resolution (res_parameter) accordingly
  
###1. set up ####
#load the library
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)

#set the output
plotdir = ""
dir.create(plotdir)
setwd(plotdir)

#parameter set up
res_parameter = 
cell_type = ""      #for files' name, please insert: "B_" , "T_NK_" , "pDC_Myeloid_", "T_NK_nonCD4_",or "T_NK_CD4_"
title_change = paste0("Subclustering_RPCA_",cell_type, "_res",  res_parameter)   #for plot titles

#load marker list
load(file='/mnt/sod2/csb6/eliora/AIDA_data_freeze_v2/script_template/marker_list.rda')

### 2. load the input####
#import the rds 
print("loading datset...")
seurat_obj <- readRDS("")
colnames(seurat_obj[[]])
print("seurat obejct is loaded")

### 3. clustering ####
print("start clustering...")
DefaultAssay(seurat_obj) <- 'integrated'
vec_features <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = vec_features)
seurat_obj <- RunPCA(seurat_obj)
png(paste0(cell_type,"elbow_plot_PC.png"), type = "cairo", height=500, width= 500)
ElbowPlot(seurat_obj,ndims =30)
dev.off()
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = res_parameter)
print("clustering is done")
saveRDS(seurat_obj,paste0(cell_type,"1.subclustering_reintegrated.rds"))
print("seurat object is saved")

### 4. clusters  and marker gene visualization ####
DefaultAssay(seurat_obj) <- 'RNA'
print("plotting UMAPs, dot plots, and violin plot...")

#### plot seurat cluster #### 
png(paste0(cell_type,"3.UMAP_seurat_cluster.png"), type = "cairo", height=800, width= 800)
DimPlot(seurat_obj,pt.size = 0.5) + labs(title = title_change)
dev.off()
png(paste0(cell_type,"3.UMAP_split_by_cluster.png"), type = "cairo", height=1000, width= 1000)
DimPlot(seurat_obj,split.by='seurat_clusters',pt.size = 0.5,ncol=6)+ labs(title = paste0(title_change,"_split_by_cluster"))+ NoLegend()
dev.off()

#### plot gene expression (dot plots, feature plots, and violin plots) ####
##### plotting basic marker #### 
png(file = file.path(plotdir,paste0(cell_type,"4.pan_marker_dot_plot.png")),width=2000, height=900, type = "cairo")
DotPlot(seurat_obj, features = basic.marker , assay = 'RNA') +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  RotatedAxis() + labs(title = paste0(title_change,"_pan_marker"))
dev.off()

##### plotting platelet' genes #### 
#checking any platelet's gene contamination 
png(file = file.path(plotdir,paste0(cell_type,"3.feature_plot_platelet.png")),width=900, height=1100, type = "cairo")
FeaturePlot(seurat_obj, slot = 'data', features = platelet_gene, ncol = 3) 
dev.off()
png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_platelet.png")),width=1500, height=900, type = "cairo")
VlnPlot(seurat_obj,features = platelet_gene, pt.size = 0, ncol = 3) 
dev.off()

##### plotting cell type specific-marker genes ####
##different lineage will have different markers to be plot

if(cell_type == "B"){
  
  ##overall
  png(file = file.path(plotdir,paste0(cell_type,"3.feature_plot_marker_genes.png")),width=1000, height=1200, type = "cairo")
  print(FeaturePlot(seurat_obj, slot = 'data', features = unlist(marker_genes_B_plasma), ncol = 4))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_overall_marker.png")),width=500, height=700, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.basic , assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_basic_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_overall_marker_splitted.png")),width=800, height=700, type = "cairo")
  print(DotPlot(seurat_obj, features = marker_genes_B_plasma , assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_overall_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_all.png")),width=1500, height=1000, type = "cairo")
  print(VlnPlot(seurat_obj,features =unlist(marker_genes_B_plasma) , pt.size = 0, ncol = 4))
  dev.off()
  
  ##naive
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_naive_transional.png")),width=600, height=700, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.naive, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_naiveB_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_naive_transional.png")),width=1700, height=900, type = "cairo")
  print(VlnPlot(seurat_obj,features =markers.naive , pt.size = 0, ncol = 4))
  dev.off()
  
  ##memory
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_memory.png")),width=600, height=700, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.memory, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_memoryB_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_memory.png")),width=2000, height=600, type = "cairo")
  print(VlnPlot(seurat_obj,features =markers.memory , pt.size = 0, ncol = 6))
  dev.off()
  
  ##atypical
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_atypical.png")),width=600, height=700, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.atypical, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_atypicalB_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_atypical.png")),width=1600, height=900, type = "cairo")
  print(VlnPlot(seurat_obj,features =markers.atypical , pt.size = 0, ncol = 4))
  dev.off()
  
  ##plasma -- checking if there is plasma cells' contamination
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_plasma.png")),width=400, height=700, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.plasma, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_plasma_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_plasma.png")),width=1300, height=900, type = "cairo")
  print(VlnPlot(seurat_obj,features =markers.plasma , pt.size = 0, ncol = 3))
  dev.off()
  
}else if (cell_type == "pDC_Myeloid") {
  
  ##overall
  png(file = file.path(plotdir,paste0(cell_type,"3.feature_plot_marker_genes.png")),width=2100, height=1200, type = "cairo")
  print(FeaturePlot(seurat_obj, slot = 'data', features = unlist(marker_genes_pDC_M), ncol = 7))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.all_dot_plot.png")),width=1200, height=800, type = "cairo")
  print(DotPlot(seurat_obj, features = marker_genes_pDC_M , assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis()+ labs(title = paste0(title_change,"_overall_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_all.png")),width=2300, height=700, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = unlist(marker_genes_pDC_M), ncol = 6, pt.size = 0))
  dev.off()
  ## DC
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_DC.png")),width=1200, height=900, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.DC, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_DC_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_DC.png")),width=3000, height=2000, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = unlist(markers.DC), ncol = 6, pt.size = 0))
  dev.off()
  ##monocytes
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_monocytes.png")),width=900, height=900, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.mono, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis()+ labs(title = paste0(title_change,"_monocytes_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_monocytes.png")),width=2200, height=1200, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = unlist(markers.mono), ncol = 5, pt.size = 0))
  dev.off()
  
}else if (cell_type %in% c("T_NK","T_NK_nonCD4_","T_NK_CD4_")) {
  
  ##overall
  png(file = file.path(plotdir,paste0(cell_type,"3.feature_plot_marker_genes.png")),width=2500, height=1500, type = "cairo")
  print(FeaturePlot(seurat_obj, slot = 'data', features = unlist(marker_genes_T_NK), ncol = 10))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.all_dot_plot.png")),width=1500, height=1200, type = "cairo")
  print(DotPlot(seurat_obj, features = marker_genes_T_NK , assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_overall_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_basic_marker.png")),width=1500, height=600, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4"), ncol = 3, pt.size = 0))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_all.png")),width=4500, height=1200, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = unlist(marker_genes_T_NK), ncol = 11, pt.size = 0))
  dev.off()
  
  ## Tregs
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_Tregs.png")),width=1000, height=1200, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.Tregs, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_Tregs_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_Tregs.png")),width=2700, height=800, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = markers.Tregs, ncol = 6, pt.size = 0))
  dev.off()
  
  ## CD8
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_CD8T.png")),width=1000, height=1200, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.CD8, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_CD8+T_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_CD8T.png")),width=2700, height=1200, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = markers.CD8, ncol = 6, pt.size = 0))
  dev.off()
  
  ##CD4
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_CD4.png")),width=1000, height=1200, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.CD4, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_CD4+T_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_CD4.png")),width=3000, height=1500, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = markers.CD4, ncol = 7, pt.size = 0))
  dev.off()
  
  ##NK + ILC
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_NK.png")),width=1000, height=1200, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.NK, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_NK_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_NK.png")),width=2700, height=800, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = unlist(markers.NK), ncol = 6, pt.size = 0))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.dot_plot_ILC.png")),width=1000, height=1200, type = "cairo")
  print(DotPlot(seurat_obj, features = markers.ILC, assay = 'RNA') +
          geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
          guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
          RotatedAxis() + labs(title = paste0(title_change,"_ILC_gene")))
  dev.off()
  png(file = file.path(plotdir,paste0(cell_type,"4.violin_plot_ILC.png")),width=2700, height=800, type = "cairo")
  print(VlnPlot(seurat_obj, slot = 'data', features = unlist(markers.ILC), ncol = 6, pt.size = 0))
  dev.off()
}

### 5. proportion per cluster ####
#create and save metadata
metadata.df <- seurat_obj[[]]
write.table(metadata.df,file.path(plotdir,paste0(cell_type,"5.metadata_allcells_information.txt")), sep = "\t")

#### check the proportion of country, ethnicity, and sex each cluster ####
column_list = c("Country", "Ethnicity", "Sex")
for(h in 1:length(column_list)){
  k = column_list[h]
  lib_list = unique(seurat_obj[[k]])
  data <- data.frame()
  for(i in as.numeric(levels(seurat_obj))){
    denominator.cluster = nrow(filter(metadata.df,seurat_clusters ==i))
    for(j in rownames(lib_list)){
      denominator.lib = nrow(filter(metadata.df,metadata.df[,k] == lib_list[j,1]))
      frequency = nrow(filter(metadata.df,seurat_clusters ==i,metadata.df[,k] ==lib_list[j,1]))
      data.temp <- data.frame(cluster = i, column = lib_list[j,1], freq = frequency, 
                              prop_lib = round(frequency/denominator.lib,4), prop_cluster = round(frequency/denominator.cluster,4))
      data = rbind(data,data.temp)
    }
  }
  write.csv(data,file.path(plotdir,paste0(cell_type,"5.Number_of_Cells_per", k ,"_percluster.csv")),row.names= FALSE)
  png(paste0(cell_type,"5.barplot_stacked_prop_",k,"_percluster_.png"), type = "cairo", height=500, width= 800)
  print(ggplot(data,aes(fill = column, y=prop_lib, x=cluster))+ 
          geom_bar(position='fill',stat='identity') + labs(title = paste0(title_change,"_stacked_barplot_",k,"_proportion_percluster")))
  dev.off()
  png(paste0(cell_type,"5.barplot_prop_",k,"_percluster.png"), type = "cairo", height=500, width= 800)
  print(ggplot(data,aes(fill = column, y=prop_lib, x=cluster))+
          geom_bar(stat='identity') + labs(title = paste0(title_change,"_barplot_",k,"_proportion_percluster")))
  dev.off()
  png(paste0(cell_type,"5.UMAP_grouped_by_",k,".png"), type = "cairo", height=500, width= 500)
  print(DimPlot(seurat_obj,group.by=k,pt.size = 0.3)+ labs(title = paste0(title_change,"_UMAP_grouped_by",k)))
  dev.off()
  png(paste0(cell_type,"5.UMAP_split_by_",k,".png"), type = "cairo", height=400, width=300*nrow(unique(seurat_obj[[k]])))
  print(DimPlot(seurat_obj,group.by=k,pt.size = 0.3, split.by = k)+ labs(title = paste0(title_change,"_UMAP_split_by",k)))
  dev.off()
  png(paste0(cell_type,"5.barplot_stacked_freq_",k,"_percluster.png"), type = "cairo", height=500, width= 800)
  print(ggplot(data,aes(fill = column, y=freq, x=cluster))+
          geom_bar(position='fill',stat='identity')+ labs(title = paste0(title_change,"_stacked_barplot_",k,"_frequency_percluster")))
  dev.off()
}

#### check the NODG and number of cells of each cluster ####
info.cluster = data.frame(table(seurat_obj$seurat_clusters))
for(i in 1:nrow(info.cluster)){
  metadata.temp = filter(metadata.df,seurat_clusters == (i-1))
  info.cluster$Mean_NODG[i]= round(mean(metadata.temp$nFeature_RNA),0)
  info.cluster$Mean_count[i]= round(mean(metadata.temp$nCount_RNA),0)
  info.cluster$Median_NODG[i]= round(median(metadata.temp$nFeature_RNA),0)
  info.cluster$Median_count[i]= round(median(metadata.temp$nCount_RNA),0)
}
colnames(info.cluster) = c("Cluster","Num_cells","Mean_NODG","Mean_count","Median_NODG","Median_count")
write.csv(info.cluster,paste0(cell_type,"5.Number_of_cell_NOGD_percluster.csv"))
png(file = file.path(plotdir,paste0(cell_type,"5.Violin_plot_NODG_percluster.png")),width=1000, height=500, type = "cairo")
VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) + NoLegend()
dev.off()
png(file = file.path(plotdir,paste0(cell_type,"5.Violin_plot_count_percluster.png")),width=1000, height=500, type = "cairo")
VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) + NoLegend() 
dev.off()

### 6. highlight cluster ######
print("plotting UMAP higlighted per cluster...")
dir = paste0(plotdir,"/",cell_type, "6.UMAP_highlight_percluster")
dir.create(dir)

for(i in 0:(length(levels(seurat_obj))-1)){
  cells = WhichCells(seurat_obj, idents = i)
  png(file.path(dir,paste0(cell_type,"6.UMAP_flagging_cluster_",i,".png")), type = "cairo",height=500, width= 500)
  print(DimPlot(seurat_obj,group.by = NULL,pt.size = 0.5,cells.highlight = cells,cols.highlight = "#DE2D26") + 
          NoLegend() + labs(title = paste0(title_change,"_UMAP_highligting_cluster_", i)) +
          theme(plot.title = element_text(size=12)))
  dev.off()
}

### 7. library proportion per cluster ####
dir.prop = paste0(plotdir,"/", cell_type,"7.Proportion_perlibrary")
dir.create(dir.prop)
dir.freq = paste0(plotdir,"/", cell_type,"7.Frequency_perlibrary")
dir.create(dir.freq)

num.cluster = c(0:(length(levels(seurat_obj))-1))
lib_list = unique(seurat_obj$Library)
freq.data <- data.frame()
prop.data <- data.frame()

for(i in 0:(length(num.cluster)-1)){
  freq.combined = c()
  prop.combined = c()
  filter.df = filter(metadata.df,seurat_clusters == i)
  for(j in 1:length(lib_list)){
    if(lib_list[j] %in% unique(filter.df$Library)){
      ncell = filter(metadata.df,Library == lib_list[j])
      freq = nrow(subset(filter.df, Library == lib_list[j]))
      prop = round(100*nrow(subset(filter.df, Library == lib_list[j]))/nrow(ncell),4)
    }else{
      prop = 0
      freq = 0}
    prop.combined = append(prop.combined,prop)
    freq.combined = append(freq.combined,freq)
  }
  if(i ==0){
    freq.data = freq.combined
    prop.data = prop.combined
  }else{
    freq.data = rbind(freq.data,freq.combined)
    prop.data = rbind(prop.data,prop.combined)
  }
}
rownames(freq.data) = paste0("cluster_",num.cluster)
rownames(prop.data) = paste0("cluster_",num.cluster)
colnames(freq.data) = lib_list
colnames(prop.data) = lib_list

freq.data = data.frame(t(freq.data))
prop.data = data.frame(t(prop.data))
freq.data$Library = rownames(freq.data)
prop.data$Library = rownames(prop.data)
freq.data$country = sapply(freq.data$Library,function(x) strsplit(x,"_")[[1]][1])
prop.data$country = sapply(prop.data$Library,function(x) strsplit(x,"_")[[1]][1])
write.csv(freq.data,file.path(plotdir,paste0(cell_type,"7.Number_of_Cells_perlib_percluster.csv")),row.names= FALSE)
write.csv(prop.data,file.path(plotdir,paste0(cell_type,"7.Prop_perlib_percluster.csv")),row.names= FALSE)

for(i in 1:length(num.cluster)){
  png(file.path(dir.prop,paste0(cell_type,"7.barplot_prop_cluster",num.cluster[i],"_perlib.png")), type = "cairo", height=400, width= 1000)
  print(ggplot(prop.data,aes(y=prop.data[,i], x=Library, fill= country))+
          geom_bar(stat='identity')+ ylab('Proportion (%)') + 
          theme(axis.text.x = element_text(hjust = 1, size = 7, angle=45),axis.text.y = element_text(size = 7),
                legend.text=element_text(size = 7),legend.title=element_text(size = 7)) + 
          labs(title = paste0(title_change,"_barplot_proportion_cells_from_cluster_", i-1, "_in_each_library")))
  dev.off()
}

for(i in 1:length(num.cluster)){
  png(file.path(dir.freq,paste0(cell_type,"7.barplot_cluster",num.cluster[i],"_perlib.png")), type = "cairo", height=400, width= 1000)
  print(ggplot(freq.data,aes(y=freq.data[,i], x=Library, fill= country))+
          geom_bar(stat='identity')+
          theme(axis.text.x = element_text(hjust = 1, size = 7, angle=45),axis.text.y = element_text(size = 7),
                legend.text=element_text(size = 7),legend.title=element_text(size = 7)) + 
          labs(title = paste0(title_change,"_barplot_frequency_cells_from_cluster_", i-1, "_in_each_library")))
  dev.off()
}
print("calculating library proportion is done")

### 8. UMAP per cluster split by metadata ####
# print UMAP per cluster split by ethnicity, country, and sex
column_list = c("Ethnicity","Country", "Sex")

for(lib in column_list){
  dir = paste0(plotdir,paste0("/",cell_type,"8.UMAP_per",lib,"_percluster"))
  dir.create(dir)
  
  for(i in 0:(length(levels(seurat_obj))-1)){
    seurat_obj_subset = subset(seurat_obj, idents = i)
    png(file.path(dir,paste0(cell_type,"8.UMAP_split_",lib, "_cluster_",i,".png")), 
        type = "cairo",height=300, width= 200*(nrow(unique(seurat_obj[[lib]]))))
    print(DimPlot(seurat_obj_subset,group.by = lib,split.by = lib,pt.size = 0.3) 
          + NoLegend() + labs(title = paste0(title_change," UMAP split by ",lib ," for cluster ", i))
          + theme(plot.title = element_text(size=12)))
    dev.off()
  }
}

### 9. DEG analysis ####
DefaultAssay(seurat_obj) <- 'RNA'
pbmc.markers <- FindAllMarkers(seurat_obj,min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers,file.path(plotdir,paste0(cell_type,"2.DEG_result_min.pct0.25_threshold0.25.csv")))
pbmc.markers.filtered = data.frame(pbmc.markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC))
write.csv(pbmc.markers.filtered,file.path(plotdir,paste0(cell_type,"2.DEG_result_filtered_top20_min.pct0.25_threshold0.25.csv")))
print("DEG is done")
gene.features = unique(data.frame(pbmc.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC))[,7])
png(file = file.path(plotdir,paste0(cell_type,"2.dot_plot_DEG_top5.png")),width=1800, height=600, type = "cairo")
DotPlot(seurat_obj, features = gene.features, assay = 'RNA') +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  RotatedAxis() + labs(title = paste0(title_change,"_top5DEG"))
dev.off()
print("DEG is done")


### 10. Checking expression of selected amrker genes (percentage , mean) ####
#marker gene
if(cell_type == "B"){
  marker_gene = marker.per.mean.B
}else if (cell_type == "pDC_Myeloid"){
  marker_gene = marker.per.mean.pDC.M
}else if (cell_type == "T_NK"){
  marker_gene = marker.per.mean.TNK
}else if (cell_type == "T_NK_CD4_"){
  marker_gene = marker.per.mean.TNK.CD4
}else if (cell_type == "T_NK_nonCD4_"){
  marker_gene = marker.per.mean.TNK.nonCD4
}
  
##check the mean expression
data <- data.frame()
num = grep("TRUE",rownames(seurat_obj) %in% marker_gene)
for(i in as.numeric(levels(seurat_obj))){
  barcode = WhichCells(seurat_obj,idents = i)
  ncell = length(barcode)
  exp.df = GetAssayData(object = seurat_obj, assay = "RNA", slot = "data")[marker_gene,barcode]
  #exp.df = exp.df[num,]
  exp.mean = apply(exp.df,1, FUN=mean)
  data.temp = round(c(i, ncell, exp.mean),4)
  if(i ==0){data = data.temp}else{data = rbind(data,data.temp)}
}
colnames(data) = c("cluster","num.cell",names(exp.mean))
rownames(data) = NULL
write.csv(data,file.path(plotdir,paste0(cell_type,"9.mean_expression_marker_Genes_percluster.csv")))

##check the mean expression - non 0
check_mean_exp <- function(vector){
  return(mean(vector[vector>0]))
}
data <- data.frame()
num = grep("TRUE",rownames(seurat_obj) %in% marker_gene)
for(i in as.numeric(levels(seurat_obj))){
  barcode = WhichCells(seurat_obj,idents = i)
  ncell = length(barcode)
  exp.df = GetAssayData(object = seurat_obj, assay = "RNA", slot = "data")[marker_gene,barcode]
  exp.mean = apply(exp.df,1, FUN=check_mean_exp)
  data.temp = round(c(i, ncell, exp.mean),4)
  if(i ==0){data = data.temp}else{data = rbind(data,data.temp)}
}
colnames(data) = c("cluster","num.cell",names(exp.mean))
rownames(data) = NULL
write.csv(data,file.path(plotdir,paste0(cell_type,"9.mean_expression_non0_marker_Genes_percluster.csv")))

## check the percentage
check_exp <- function(vector){
  return(sum(vector>0))
}
data <- data.frame()
num = grep("TRUE",rownames(seurat_obj) %in% marker_gene)
for(i in as.numeric(levels(seurat_obj))){
  barcode = WhichCells(seurat_obj,idents = i)
  ncell = length(barcode)
  exp.df = GetAssayData(object = seurat_obj, assay = "RNA", slot = "data")[marker_gene,barcode]
  num.cell = apply(exp.df,1, FUN=check_exp)
  data.temp = round(c(i, ncell, round(num.cell/ncell,2)),4)
  if(i ==0){data = data.temp}else{data = rbind(data,data.temp)}
}
colnames(data) = c("cluster","num.cell",names(num.cell))
rownames(data) = NULL
write.csv(data,file.path(plotdir,paste0(cell_type,"9.percentage_marker_Genes_percluster.csv")))
