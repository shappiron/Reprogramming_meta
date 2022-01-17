corplot <- function(data, type='full_intersection', method='spearman', top=500, pvplot=F,
                    label_size=1.2, value_size=0.8, title="Spearman correlation",
                    measure="logFC", criterion="FDR", plotw=13, ploth=13){
    #type = {full_intersection, pairwise_top}
    if (type == 'full_intersection'){
        #provided data must be a data.frame type
        cm <- round(cor(data, method='spearman'), 2)
        cellnote <- cm
    }else if (type == 'pairwise_top'){
        #pairwise intersection of top N genes or full in case of top==F
        #provided data must be a list of data frames
        cm <- matrix(0:0, nrow=length(data), ncol=length(data))
        pm <- matrix(0:0, nrow=length(data), ncol=length(data))
        for (i in 1:length(data)){
            for (j in 1:length(data)){
                if (top > 0) {
                    #intersection of all genes
                    I <- intersect(rownames(data[[i]]), rownames(data[[j]]))
                    d1 <- data[[i]][I,]
                    d2 <- data[[j]][I,]
                    #union of top genes
                    nx <- min(c(top, nrow(d1)))
                    ny <- min(c(top, nrow(d2)))
                    rx <- rownames(d1[order(d1[[criterion]]),][1:nx,])
                    ry <- rownames(d2[order(d2[[criterion]]),][1:ny,])
                    r <- union(rx, ry)
                    x <- d1[r,][measure]
                    y <- d2[r,][measure]
                }else{
                    x <- data[[i]][measure]
                    y <- data[[j]][measure]   
                }
                pair <- transform(merge(x, y, by="row.names", sort=FALSE), row.names=Row.names, Row.names=NULL)
                #pair <- pair[complete.cases(pair),]
                if (nrow(pair) > 1){
                    test <- cor.test(pair[[1]], pair[[2]], method=method)
                    sp <- test$estimate
                    pv <- test$p.value
                }else{
                    sp <- NA
                }
                if (is.na(sp)){
                    cm[i, j] <- 0
                    pm[i, j] <- 1.0
                }else{
                    cm[i, j] <- sp
                    pm[i, j] <- pv
                }    
            }
        }
        cm <- round(cm, 2)
        rownames(cm) <- names(data) 
        colnames(cm) <- names(data)
        rownames(pm) <- names(data) 
        colnames(pm) <- names(data)
        #prepare value annotation for heatmap
        FDR_table <- pm
        #extract upper triangular of the matrix and apply FDR correction
        tmp <- matrix(0, nrow=nrow(FDR_table), ncol=ncol(FDR_table))
        padj <- matrix(p.adjust(FDR_table[upper.tri(FDR_table)], method="BH"))
        tmp[upper.tri(tmp)] <- padj
        tmp <- tmp + t(tmp)
        diag(tmp) <- 0
        FDR_table <- as.data.frame(tmp, nrow=nrow(corr_matrix_list$P.Value))

        cellnote <- cm
        cellnote[FDR_table >= 0.05] <- NA

    }else{
        print("Please, choose a type")
    }
    if (pvplot==T){cm <- log10(pm)}
    
    #plotting
    options(repr.plot.width = plotw, repr.plot.height = ploth)
    par(mar=c(7,4,4,2)+0.1)
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 300)
    heatmap.2(cm, cellnote=cellnote, notecex=value_size, main=title,
            density.info="density", denscol="black", key.par=list(mar=c(7,4,4,2)+0.1), key.ylab=NA,
            notecol="black", srtCol=45, trace="none", 
            cexRow=label_size, cexCol=label_size, margins=c(18,18), col=my_palette,
            breaks=seq(-1,1, length.out=301))

    return(list(cm=cm))
}
