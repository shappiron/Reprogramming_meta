#
gene_barplot <- function(gsea, func_interest){
    #combine list of signatures into list of dataframes where each dataframe is the same column type
    SL <- list()
    cols <- c('NES', 'pval', 'padj')
    types <- c(rep('Reprogramming',3), rep('Interventions',6), rep('Aging', 7))
    
    for (col in cols){
        SL[[col]] <- data.frame(row.names=rownames(gsea[[1]]))
        for (name in names(gsea)){
                tmp <- gsea[[name]][col]
                colnames(tmp) <- c(name)
                SL[[col]] <- transform(merge(SL[[col]], tmp, 
                                                    by=0, sort=F,), 
                                                    row.names=Row.names, Row.names=NULL)
        }
    }

    tmpdf <- gsea[[1]]

    melter <- data.frame()
    for (f in func_interest){
            #func for mapping stars
            plab <- function(p )ifelse(p<0.001,'***', ifelse(p<0.01,'**', ifelse(p<0.05,'*',''))) 

            #rescale columns
            tmpLogFC = SL[['NES']] #/ apply(SL[['NES']], 2, sd)

            mel <- t(tmpLogFC[g,])
            mel <- merge(mel,  t(SL[['FDR']][g,])), by=0, sort=F)
            colnames(mel) <- c('Signature', 'Coef', 'SE', 'FDR')
            mel$Type <- types
            mel$Signature <- factor(as.character(mel$Signature), levels=unique(as.character(mel$Signature)))
            mel$Star <- plab(mel$FDR)
            mel$Vjust <- if_else(mel$Coef > 0, -0.6, 1.7)
            mel$Symbol <- rep(sym, nrow(mel))
            melter <- rbind(melter, mel)
    }
    melter$Symbol <- factor(melter$Symbol, levels=func_interest)

    #plot
    #suppressMessages(library(ggrepel))
    g <- ggplot(melter, aes(x=Signature, y=Coef, fill=Type))+ 
            geom_bar(stat="identity") +
            facet_wrap(~Symbol, ncol=1, scales = "free_y", dir="v") + 
            geom_text(aes(label=Star, vjust=melter$Vjust, hjust=0.5), fontface='bold', size=7)+
            geom_errorbar(aes(ymin=Coef-SE, ymax=Coef+SE), width=0.18) + 
            theme(axis.text=element_text(size=14, angle=45, hjust=1),
                    strip.background = element_blank(),
                    strip.placement = 'inside',
                    strip.text = element_text(colour = 'black', size=16, face="bold", hjust = 0.5),     
                    axis.title=element_text(size=18, face="bold"), 
                    plot.title=element_text(size = 20, face = "bold"),
                    legend.title=element_text(size=16, face = "bold"),
                    legend.text=element_text(size=14),
                    legend.position="top",
                    panel.background = element_rect(fill='white', colour='black', size=1.1, linetype='solid'),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.margin = unit(c(0.5,0.5,0.5,0.), "cm")) +
            labs(x = "", y = "Normalized meta slope")
    return(g)
}                