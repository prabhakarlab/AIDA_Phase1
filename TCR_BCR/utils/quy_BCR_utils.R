tcr1_lines <- c("TCR1", "cdr3_aa1", "cdr3_nt1")
tcr2_lines <- c("TCR2", "cdr3_aa2", "cdr3_nt2")
heavy_lines <- c("IGH", "cdr3_aa1", "cdr3_nt1", "vgene1")
light_lines <- c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")

BCRprofle <- function(bcr.res){
    chain1 <- "IGH"
    chain2 <- "IGK"
    chain3 <- "IGL"
    cellType <- "B"
    
    bcr.res <- subset(bcr.res, chain != "Multi")
    bcr.res <- subset(bcr.res, chain %in% c(chain1, chain2, chain3))
    bcr.res <- subset(bcr.res, productive %in% c(TRUE, "TRUE", "True","true"))
    
    bcr.res <- organizeGenes(cellType, bcr.res, chain1, chain2, chain3)
    # keeps cells at least one heavy chain and one light chain
    bcr.res <- completeBCR(bcr.res)

    unique_df <- unique(bcr.res$barcode)
    Con.df <- data.frame(matrix(NA, length(unique_df), 9))
    colnames(Con.df) <- c("barcode",heavy_lines, light_lines)
    
    Con.df$barcode <- unique_df
    Con.df <- organizeBCR(Con.df, bcr.res)
    Con.df$BCRgene <- paste(Con.df$IGH, Con.df$IGLC, sep = "_")
    return(Con.df)
}



#Sorting the V/D/J/C gene sequences for T and B cells
organizeGenes <- function(cellType, data2, chain1, chain2, chain3) {
    if(cellType %in% c("T-AB", "T-GD")) {
        data2 <- data2 %>% 
            mutate(TCR1 = ifelse(chain == chain1, paste(with(data2, 
                                                             interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(TCR2 = ifelse(chain == chain2, paste(with(data2, 
                                                             interaction(v_gene,  j_gene, d_gene, c_gene))), NA))
    }
    else {
        data2 <- data2 %>% 
            mutate(IGKct = ifelse(chain == "IGK", paste(with(data2, 
                                                             interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(IGLct = ifelse(chain == "IGL", paste(with(data2, 
                                                             interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(IGHct = ifelse(chain == "IGH", paste(with(data2, 
                                                             interaction(v_gene, j_gene, d_gene, c_gene))), NA))
    }
    return(data2)
    
}


completeBCR <- function(bcr.res){ # select BCR with heavy and light chain present
    unique.barcode <- unique(bcr.res$barcode)
    
    complete.barcode <- c()
    for (i in seq_along(unique.barcode)){
        barcode.i <- unique.barcode[i]
        bcr.res.i <- bcr.res[which(bcr.res$barcode==barcode.i),]
        if (nrow(bcr.res.i)>1){
            igh <- FALSE
            igl <- FALSE
            igk <- FALSE
            if ("IGH" %in% bcr.res.i$chain ){
                igh <- TRUE
            }
            if ("IGL" %in% bcr.res.i$chain ){
                igl <- TRUE
            }
            if ("IGK" %in% bcr.res.i$chain ){
                igk <- TRUE
            }
            
            if (igh & (igl | igk)){
                complete.barcode <- c(complete.barcode, barcode.i)
            }
        }
    }
    
    bcr.res.com <- bcr.res[which(bcr.res$barcode %in% complete.barcode),]
    return(bcr.res.com)
}


organizeBCR <- function(Con.df, bcr.res){
    for (i in seq_along(Con.df$barcode)){
        barcode.i <- Con.df$barcode[i]
        bcr.res.i <- bcr.res[which(bcr.res$barcode==barcode.i),]
        
        # heavy chain
        heavy.i <- bcr.res.i[which(bcr.res.i$chain=="IGH"),]
        # if more than two heavy chains, select one with most abundant umis
        heavy.i <- heavy.i[order(heavy.i$umis, decreasing = TRUE),]
        heavy.i <- heavy.i[1,]
        Con.df[which(Con.df$barcode==barcode.i),heavy_lines] <- heavy.i[,c("IGHct", "cdr3", "cdr3_nt", "v_gene")]
        
        # light chain
        light.i <- bcr.res.i[which(bcr.res.i$chain %in% c("IGL","IGK")),]
        # if more than two light chains, select one with most abundant umis
        light.i <- light.i[order(light.i$umis, decreasing = TRUE),]
        light.i <- light.i[1,]
        if (light.i$chain == "IGL"){
            Con.df[which(Con.df$barcode==barcode.i),light_lines] <- light.i[,c("IGLct", "cdr3", "cdr3_nt", "v_gene")]
        } else {
            Con.df[which(Con.df$barcode==barcode.i),light_lines] <- light.i[,c("IGKct", "cdr3", "cdr3_nt", "v_gene")]
        }
    }
    return(Con.df)
}


calClonalFreq <- function(tcr.res, cloneType="TCRgene"){
    tcr.res.sub <- tcr.res[,cloneType,drop=FALSE]
    tcr.res.sub$Frequency <- 1
    tcr.res.sub <- aggregate(Frequency ~ ., tcr.res.sub, FUN = sum)
    tcr.res.merge <- merge(tcr.res, tcr.res.sub, by = cloneType, all = TRUE)
    PreMeta <- unique(tcr.res.merge[,c("barcode",cloneType, "Frequency")])
    PreMeta$Proportion <- PreMeta$Frequency/nrow(PreMeta)
    return(PreMeta)
}



TCRgeneCompMatrix <- function(tcr.barcode){
    Vgene <- c()
    Jgene <- c()
    for (i in seq_along(tcr.barcode$TCRgene)){
        TCRgene.i <- tcr.barcode$TCRgene[i]
        TCRgene.i.gene <- unlist(strsplit(TCRgene.i,split = "\\_|\\.|\\;"))
        TCRgene.i.Vgene <- TCRgene.i.gene[grepl(pattern = "^TRAV|^TRBV", TCRgene.i.gene)]
        TCRgene.i.Jgene <- TCRgene.i.gene[grepl(pattern = "^TRAJ|^TRBJ", TCRgene.i.gene)]
        
        Vgene <- c(Vgene, TCRgene.i.Vgene)
        Jgene <- c(Jgene, TCRgene.i.Jgene)
    }
    
    Vgene <- unique(Vgene)
    Jgene <- unique(Jgene)
    
    Vgene <- Vgene[order(nchar(Vgene), Vgene)]
    Jgene <- Jgene[order(nchar(Jgene), Jgene)]
    
    gene.matrix <- as.data.frame(matrix(data = 0,nrow = length(Jgene),
                                        ncol = length(Vgene)))
    colnames(gene.matrix) <- Vgene
    rownames(gene.matrix) <- Jgene
    
    for (i in seq_along(tcr.barcode$TCRgene)){
        TCRgene.i <- tcr.barcode$TCRgene[i]
        TRA.i <- unlist(strsplit(TCRgene.i,split = "_"))[1]
        TRB.i <- unlist(strsplit(TCRgene.i,split = "_"))[2]
        
        # for each TRA
        TRA.i.all <- unlist(strsplit(TRA.i, split = ";"))
        for (j in seq_along(TRA.i.all)){
            TRA.i.j <- TRA.i.all[j]
            TRA.i.j.all <- unlist(strsplit(TRA.i.j, split = "\\."))
            TRA.i.j.Vgene <- TRA.i.j.all[grepl(pattern = "^TRAV", x = TRA.i.j.all)]
            TRA.i.j.Jgene <- TRA.i.j.all[grepl(pattern = "^TRAJ", x = TRA.i.j.all)]
            gene.matrix[TRA.i.j.Jgene, TRA.i.j.Vgene] <- gene.matrix[TRA.i.j.Jgene, TRA.i.j.Vgene] + 1
        }
        
        # for each TRB
        TRB.i.all <- unlist(strsplit(TRB.i, split = ";"))
        for (j in seq_along(TRB.i.all)){
            TRB.i.j <- TRB.i.all[j]
            TRB.i.j.all <- unlist(strsplit(TRB.i.j, split = "\\."))
            TRB.i.j.Vgene <- TRB.i.j.all[grepl(pattern = "^TRBV", x = TRB.i.j.all)]
            TRB.i.j.Jgene <- TRB.i.j.all[grepl(pattern = "^TRBJ", x = TRB.i.j.all)]
            gene.matrix[TRB.i.j.Jgene, TRB.i.j.Vgene] <- gene.matrix[TRB.i.j.Jgene, TRB.i.j.Vgene] + 1
        }
        
    }
    
    return(gene.matrix)
}



BCRgeneCompMatrix <- function(bcr.barcode){
    Vgene <- c()
    Jgene <- c()
    for (i in seq_along(bcr.barcode$BCRgene)){
        BCRgene.i <- bcr.barcode$BCRgene[i]
        BCRgene.i.gene <- unlist(strsplit(BCRgene.i,split = "\\_|\\.|\\;"))
        BCRgene.i.Vgene <- BCRgene.i.gene[grepl(pattern = "^IGHV|^IGKV|^IGLV", BCRgene.i.gene)]
        BCRgene.i.Jgene <- BCRgene.i.gene[grepl(pattern = "^IGHJ|^IGKJ|^IGLJ", BCRgene.i.gene)]
        
        Vgene <- c(Vgene, BCRgene.i.Vgene)
        Jgene <- c(Jgene, BCRgene.i.Jgene)
    }
    
    Vgene <- unique(Vgene)
    Jgene <- unique(Jgene)
    
    Vgene <- Vgene[order(nchar(Vgene), Vgene)]
    Jgene <- Jgene[order(nchar(Jgene), Jgene)]
    
    gene.matrix <- as.data.frame(matrix(data = 0,nrow = length(Jgene),
                                        ncol = length(Vgene)))
    colnames(gene.matrix) <- Vgene
    rownames(gene.matrix) <- Jgene
    
    for (i in seq_along(bcr.barcode$BCRgene)){
        BCRgene.i <- bcr.barcode$BCRgene[i]
        IGH.i <- unlist(strsplit(BCRgene.i,split = "_"))[1]
        IGLC.i <- unlist(strsplit(BCRgene.i,split = "_"))[2]
        
        # for each heavy chain
        IGH.i.all <- unlist(strsplit(IGH.i, split = ";"))
        for (j in seq_along(IGH.i.all)){
            IGH.i.j <- IGH.i.all[j]
            IGH.i.j.all <- unlist(strsplit(IGH.i.j, split = "\\."))
            IGH.i.j.Vgene <- IGH.i.j.all[grepl(pattern = "^IGHV", x = IGH.i.j.all)]
            IGH.i.j.Jgene <- IGH.i.j.all[grepl(pattern = "^IGHJ", x = IGH.i.j.all)]
            gene.matrix[IGH.i.j.Jgene, IGH.i.j.Vgene] <- gene.matrix[IGH.i.j.Jgene, IGH.i.j.Vgene] + 1
        }
        
        # for each light chain
        IGLC.i.all <- unlist(strsplit(IGLC.i, split = ";"))
        for (j in seq_along(IGLC.i.all)){
            IGLC.i.j <- IGLC.i.all[j]
            IGLC.i.j.all <- unlist(strsplit(IGLC.i.j, split = "\\."))
            IGLC.i.j.Vgene <- IGLC.i.j.all[grepl(pattern = "^IGKV|^IGLV", x = IGLC.i.j.all)]
            IGLC.i.j.Jgene <- IGLC.i.j.all[grepl(pattern = "^IGKJ|^IGLJ", x = IGLC.i.j.all)]
            gene.matrix[IGLC.i.j.Jgene, IGLC.i.j.Vgene] <- gene.matrix[IGLC.i.j.Jgene, IGLC.i.j.Vgene] + 1
        }
        
    }
    
    return(gene.matrix)
}

