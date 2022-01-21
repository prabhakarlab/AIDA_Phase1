tcr1_lines <- c("TCR1", "cdr3_aa1", "cdr3_nt1")
tcr2_lines <- c("TCR2", "cdr3_aa2", "cdr3_nt2")

TCRprofle <- function(tcr.res){
    chain1 <- "TRA"
    chain2 <- "TRB"
    cellType <- "T-AB"
    
    tcr.res <- subset(tcr.res, chain != "Multi")
    tcr.res <- subset(tcr.res, chain %in% c(chain1, chain2))
    tcr.res <- subset(tcr.res, productive %in% c(TRUE, "TRUE", "True","true"))
    
    tcr.res <- organizeGenes(cellType, tcr.res, chain1, chain2)
    # keeps cells at least one TRA and one TRB
    tcr.res <- completeTCR(tcr.res)
    
    
    unique_df <- unique(tcr.res$barcode)
    Con.df <- data.frame(matrix(NA, length(unique_df), 7))
    colnames(Con.df) <- c("barcode",tcr1_lines, tcr2_lines)
    
    Con.df$barcode <- unique_df
    Con.df <- organizeTCR(Con.df, tcr.res)
    Con.df$TCRgene <- paste(Con.df$TCR1, Con.df$TCR2, sep = "_")
    Con.df$cdr3_aa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep = "_")
    Con.df$cdr3_nt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep = "_")
    return(Con.df)
}




#Sorting the V/D/J/C gene sequences for T and B cells
organizeGenes <- function(cellType, data2, chain1, chain2) {
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

# make sure cells at least one TRA and one TRB
completeTCR <- function(tcr.res){
    unique.barcode <- unique(tcr.res$barcode)
    
    complete.barcode <- c()
    for (i in seq_along(unique.barcode)){
        barcode.i <- unique.barcode[i]
        tcr.res.i <- tcr.res[which(tcr.res$barcode==barcode.i),]
        tra.i <- tcr.res.i$TCR1[!is.na(tcr.res.i$TCR1)]
        trb.i <- tcr.res.i$TCR2[!is.na(tcr.res.i$TCR2)]
        if (length(tra.i)!=0 &  length(trb.i)!=0){
            tcr.res.i.sub <- tcr.res.i[(!(is.na(tcr.res.i$TCR1)) & !(is.na(tcr.res.i$TCR2))),]
            if (nrow(tcr.res.i.sub) == 0){
                complete.barcode <- c(complete.barcode, barcode.i)
            }
        }
    }
    
    tcr.res.com <- tcr.res[which(tcr.res$barcode %in% complete.barcode),]
    return(tcr.res.com)
}


organizeTCR <- function(Con.df, tcr.res){
    for (i in seq_along(Con.df$barcode)){
        barcode.i <- Con.df$barcode[i]
        tcr.res.i <- tcr.res[which(tcr.res$barcode==barcode.i),]
        
        tcr.res.i.tra <- tcr.res.i[!(is.na(tcr.res.i$TCR1)),]
        # order rows by gene names
        if (nrow(tcr.res.i.tra) >= 3){
            tcr.res.i.tra <- tcr.res.i.tra[order(tcr.res.i.tra$umis, decreasing = TRUE),]
            tcr.res.i.tra <- tcr.res.i.tra[c(1,2),]
        }
        tcr.res.i.tra <- tcr.res.i.tra[order(nchar(tcr.res.i.tra$TCR1), tcr.res.i.tra$TCR1),]
        Con.df[which(Con.df$barcode==barcode.i),"TCR1"] <- paste(tcr.res.i.tra$TCR1,collapse = ";")
        Con.df[which(Con.df$barcode==barcode.i),"cdr3_aa1"] <- paste(tcr.res.i.tra$cdr3,collapse = ";")
        Con.df[which(Con.df$barcode==barcode.i),"cdr3_nt1"] <- paste(tcr.res.i.tra$cdr3_nt,collapse = ";")
        
        tcr.res.i.trb <- tcr.res.i[!(is.na(tcr.res.i$TCR2)),]
        # select only TRB with top 
        tcr.res.i.trb <- tcr.res.i.trb[order(tcr.res.i.trb$umis, decreasing = TRUE),]
        tcr.res.i.trb <- tcr.res.i.trb[1,]

        Con.df[which(Con.df$barcode==barcode.i),"TCR2"] <- tcr.res.i.trb$TCR2
        Con.df[which(Con.df$barcode==barcode.i),"cdr3_aa2"] <- tcr.res.i.trb$cdr3
        Con.df[which(Con.df$barcode==barcode.i),"cdr3_nt2"] <- tcr.res.i.trb$cdr3_nt
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


getTCRCDR3fasta <- function(x){
    fasta.out <- c()
    for (i in 1:nrow(x)){
        x.i <- x[i,]
        barcode.i <- x.i[1,1]
        cdr3.i <- x.i[1,2]
        cdr3.i.all <- unlist(strsplit(cdr3.i, split = ";"))
        for (j in seq_along(cdr3.i.all)){
            barcode.i.j <- paste0(barcode.i,"-",j,split="")
            fasta.out <- c(fasta.out, paste0(">",barcode.i.j,""))
            fasta.out <- c(fasta.out, cdr3.i.all[j])
        }
    }
    return(fasta.out)
}
