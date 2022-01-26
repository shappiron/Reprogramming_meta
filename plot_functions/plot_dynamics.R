#Plot top p-value gene dynamics
plot_dynamics <- function(rep, ultradf, top=5, interval=2, selected_genes=c(), fontsize=18){
    if (length(selected_genes)==0){
        #choose top genes
        pos_genes <- rep[rep$logFC > 0,] 
        pos_genes <- pos_genes[order(pos_genes$P.Value),][1:top,]
        neg_genes <- rep[rep$logFC < 0,] 
        neg_genes <- neg_genes[order(neg_genes$P.Value),][1:top,]
        genes <- rbind(pos_genes, neg_genes)
    }else{
        genes <- rep[rep$symbol %in% selected_genes,]
    }
        
    #construct dataset
    gdb <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c('data', 'time', 'dataset', 'gene', 'color'))))
    
    for (g in rownames(genes)){
        symbol <- as.character(rep[g,]$symbol)
        sign <- sign(rep[rep$symbol %in% c(symbol),]$logFC)
        color <- ifelse((sign > 0), 'red', 'blue')
        db = data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c('data', 'time', 'dataset', 'gene', 'color'))))
        for (name in names(ultradf)){
            data <- as.numeric(ultradf[[name]]$data[g,])
            time <- ultradf[[name]]$time
            data <- (data - mean(data)) / sd(data)
            tmp <- data.frame(list("data"=data, "time"=time, "dataset"=rep(name, length(time)), "gene"=NA))
            db <- rbind(db, tmp)
        }
        db['group_time'] <- db['time'] %/% interval * interval
        db$gene <- rep(symbol, nrow(db))
        db$color <- rep(color, nrow(db))
        gdb <- rbind(gdb, db)

    }
    gdb$color <- factor(gdb$color, levels=c('red', 'blue'))
    gdb$gene <- factor(gdb$gene, levels=as.character(genes$symbol))
    
    #plot
    g <- ggplot(gdb, aes(x = factor(group_time), y = data, fill = color)) + 
    geom_boxplot(alpha = 0.70) + 
    geom_point(aes(fill = color), size = 1, shape = 21, position = position_jitterdodge()) +
    geom_hline(yintercept=0, linetype="solid", color = "black", size=0.2, alpha=0.2)+
    facet_wrap(~reorder(gene, color), ncol=2, dir="v") + 
    scale_fill_manual(values = c("red2", "blue2")) + 
    theme(  axis.text=element_text(size=fontsize),
            axis.title=element_text(size=fontsize, face="bold"),
            strip.background = element_blank(),
            strip.placement = 'inside',
            strip.text = element_text(colour = 'black', size=fontsize, face="bold", hjust = 0.5),
            panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(1,1,1,0.5), "cm")) + 
            labs(x = "Time, days", y = "Normalized expression")
    return(g)
}