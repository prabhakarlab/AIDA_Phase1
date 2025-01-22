## A. marker gene list for ploting gene expression #####
### 1. BASIC ####
basic.marker <- list('monocytes' = c("LYZ",'CD14','CSF3R','S100A8','S100A9','FCGR3A','CDKN1C','CX3CR1'),
                     'DC' = c("CD74", "HLA-DPA1", "HLA-DPB1","CLEC9A","THBD","CLEC10A","CD1C","AXL","SIGLEC6",'LILRA4','MZB1'),
                     'T' = c("CD3D", "CD3E", "CD3G","CD4","CD8A","CD8B","FOXP3","TRGV9","TRAV1-2"),
                     'NK' = c("NKG7","GNLY","FGFBP2","SPON2","FCER1G","NCAM1","KLRC1","XCL2"), 
                     'B' = c("MS4A1","CD27","IGHD","IGHM","TCL1A","FCRL5"), 
                     'plasma' = c("JCHAIN","TNFRSF17","ITM2C"),
                     'platelet' = c("ITGA2B","PF4","TUBB1","PPBP","GP1BB","NRGN"))
platelet_gene = c("GP1BB","ITGA2B","PF4","TUBB1","NRGN","CLU","CAVIN2","PRKAR2B","SPARC","CCL5","PPBP")


### 2. B + plasma cells ####
marker_genes_B_plasma <- list( c("MS4A1","CD79A"),
                               'plasma' = c("MZB1","TNFRSF17"),
                               'naive' = c("TCL1A","IGHD"),
                               'memory' = c("CD27","TNFRSF13B","IGHM"),
                               'atypical' = c("FGR","FCRL5","FCRL3","ITGAX"),
                               'other_genes' = c("PPBP","GP1BB","IFI6","IFI44L","MX1","IFIT3"))

markers.basic = c("MS4A1","CD27","TNFRSF13B","IGHD","IGHM","TCL1A","FGR","FCRL5","ITGAX","PPBP")
markers.naive = c("MS4A1","CD27","TNFRSF13B","TCL1A","CD38","PPP1R14A","PLD4","SOX4","VPREB3","CD72","IL4R","PLPP5","SELL","FCER2")
markers.memory = c("MS4A1","TCL1A","CD27","TNFRSF13B","IGHD","IGHM","LINC01857","COBLL1","CD1C","AP3B1","CRIP1","CRIP2","TAGLN2","ITGB1","S100A4","S100A10","AHNAK")
markers.atypical <- c("MS4A1","CD27","TNFRSF13B","IGHD","IGHM","TCL1A",
                      "FGR","FCRL5","FCRL3","TNFRSF1B","ITGAX","PPP1R14A","MPP6","RGS2","HCK","RHOB","LTB","CD24")
markers.plasma <- c("MZB1","JCHAIN","TNFRSF17","ITM2C","DERL3","POU2AF1","IGHA1","TXNDC11","CD79A")

### 3. pDC + myeloid ####
marker_genes_pDC_M <- list( c("CD68","LYZ"), 
                            'pDC' = c('ITM2C','LILRA4','MZB1'),
                            'ASDC' = c('AXL','SIGLEC6','CD22'),
                            'cDC1' = c('CLEC9A','THBD','BATF3'),
                            'cDC2' = c('CLEC10A','FCER1A','CD1C'),
                            'CD14+mono' = c('CD14','CSF3R','S100A8'),
                            'CD16+mono' = c('FCGR3A','CDKN1C','CX3CR1'),
                            'other gene' = c("C1QA","PPBP","GP1BB","IFI6","IFI44L","MX1","IFIT3"))

markers.DC<- list('ASDC' = c('AXL','SIGLEC6', 'CD22', 'DAB2', 'GAS6', 'PPP1R14A'),
                  'pDC' = c('ITM2C','LILRA4','IL3RA','MZB1','IRF4'),
                  'cDC' = c('IFI30','ITGAX','FGR','LYZ','LY86','GLIPR2','ENTPD1'),
                  'cDC1' = c('CLEC9A','THBD','BATF3','C1orf54','CADM1'),
                  'cDC2' = c('CLEC10A','FCER1A','CD1C','HLA-DQA1','CD14',"S100A9","S100A8"))


markers.mono <- list( 'CD14+mono' = c("CD14","CSF3R","S100A9","S100A8","S100A12","FCN1"),
                      'CD16+mono' = c("FCGR3A","CX3CR1","CDKN1C","CSF1R","ITGAL"),
                      'other gene' = c("C1QA","C1QB","C1QC","PPBP","GP1BB","IFI6","IFI44L","MX1","IFIT3"))

### 4. T and NK cells + ILC + dnT ####
marker_genes_T_NK <- list( c("CD3D", "CD3E", "CD3G","CD4","CD8A","CD8B"),
                           'Tregs' = c("FOXP3","CTLA4","TCF7","HLA-DRB1"),
                           'gdT' = c("TRDC","TRGV9","TRDV2","TRDV1","TRDV3"),
                           'MAIT' = c("TRAV1-2","NCR3","SLC4A10"),
                           'ILC' = c("SOX4",'TNFRSF18','LINC01229','TTLL10','KIT'),
                           'dnT' = c('FXYD2','AC004585.1','NUCB2','PTPN3','MIR4422HG'),
                           'common T' = c("CCR7","SELL","LEF1","IL7R","ITGB1","KLRB1","CCL5","GZMK",
                                          "GZMB","GZMH","GZMA","CST7"),
                           'NK' = c("NKG7","GNLY","TYROBP","FCGR3A","FGFBP2","SPON2","FCER1G","NCAM1","KLRC1","XCL2"),
                           'other gene' = c("MKI67","STMN1","RRM2","PPBP","GP1BB","IFI6","IFI44L","MX1","IFIT3"))

markers.Tregs<- c("CD3D", "CD3E", "CD3G","CD4","CD8A","CD8B",
                  "FOXP3","RTKN2","IL2RA","CTLA4",
                  "CD74","HLA-DRB1","HLA-DRB5","HLA-DRA",
                  "TCF7","CCR7","EEF1B2","PLAC8",
                  "LIMS1","KLRB1","AQP3","GPR25")

markers.CD8 <- c("CD3D", "CD3E", "CD3G", "CD8A","CD8B","CD4",
                 "MKI67","KLRB1","TRAV1-2", "NCR3", "SLC4A10",
                 "CCR7","SELL","LEF1","TCF7",
                 "CD27","IL7R","GZMK","CCL5",
                 "CST7","NKG7","GZMA","PRF1","GZMH","GZMB","GNLY","KLRD1","TBX21","CX3CR1","FGFBP2","LAG3","IFNG",
                 "SOX4","IFI6","IFI44L","PPBP")

markers.CD4<- c("CD3D", "CD3E", "CD3G", "CD4","CD8A","CD8B",
                "MKI67","FOXP3","RTKN2","IL2RA",
                "CCR7","SELL","LEF1","AIF1","BACH2",
                "MAL", "CD27",
                "LTB","AQP3","KLF6","ITGB1","PASK","S100A4","KLRB1","AHNAK",
                "CCL5","LYAR","ITGB7","GZMA","GZMK","GZMB","GZMH","CST7","GNLY","NKG7","FGFBP2",
                "PPBP","CXCR5","LTA","TNF","CXCR3","GATA3","CCR4","RORA","IFI6","IFI44L","SOX4")

markers.NK<- list("CD3D", "CD3E", "CD3G","CD4","CD8A","CD8B","IFI6","IFI44L","PPBP","GP1BB","MKI67","STMN1","RRM2",
                  'NK' = c("TYROBP","CD247","NKG7","GNLY","KLRD1","KLRF1",'GZMB'),
                  'NK CD16' = c( "FCGR3A","FGFBP2","SPON2","FCER1G","KLRC2"),
                  'NK CD56' = c("NCAM1","KLRC1","XCL1","XCL2"))

markers.ILC<- list('ILC' = c("SOX4",'TNFRSF18','LINC01229','TTLL10','KIT'),
                   'ILC1' = c("CD3D","CD27","GZMK","TNFRSF8","TNFRSF1B","TNFRSF10A","CCL5","LAG3"),
                   'ILC2' = c("KLRG1","GATA3","MAF","HPGD","IL10RA","IL32"),
                   'ILC3' = c("HLA-DRA","HLA-DRB1","XCL1","XCL2","KLRC1","FCER1G"))

## B. marker gene list for checking the eprcetnage and mean expression #####
marker.per.mean.B = c("MS4A1","CD27","TNFRSF13B","IGHD","IGHM","TCL1A","FGR","FCRL5","FCRL3")
marker.per.mean.pDC.M = c("CD14","FCGR3A","CD1C","CLEC10A","FCER1A","CLEC9A","THBD","MZB1","LILRA4","AXL","SIGLEC6")
marker.per.mean.TNK = c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","TRGV9","TRAV1-2","FOXP3","FCGR3A","NCAM1","FCER1G")
marker.per.mean.TNK.CD4 = c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CCR7","ITGB1","KLRB1","CCL5","GZMK","GZMB","FOXP3","CTLA4",'FXYD2','HLA-DRB1',"TCF7","LIMS1")
marker.per.mean.TNK.nonCD4 = c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","TRGV9","TRAV1-2","CCR7","CCL5","GZMK","GZMB","FCGR3A","NCAM1","FCER1G","KLRC2")

## C. save all markers in 1 rda file ####
setwd("/mnt/sod2/csb6/eliora/AIDA_data_freeze_v2/script_template")
save(basic.marker,platelet_gene,marker_genes_B_plasma,markers.basic,markers.naive,markers.memory,markers.atypical,markers.plasma,
     marker_genes_pDC_M,markers.DC,markers.mono,
     marker_genes_T_NK,markers.Tregs,markers.CD8,markers.CD4,markers.NK,markers.ILC, 
     marker.per.mean.B, marker.per.mean.pDC.M, marker.per.mean.TNK, marker.per.mean.TNK.CD4, marker.per.mean.TNK.nonCD4,
     file = "marker_list.rda")
