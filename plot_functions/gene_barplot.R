#
gene_barplot <- function(combined_full_rep, gene_interest, fontsize=18){
    #combine list of signatures into list of dataframes where each dataframe is the same column type
    SL <- list()
    cols <- c('logFC', 'SE', 'P.Value', 'FDR')
    types <- c(rep('Reprogramming',3), rep('Interventions',6), rep('Aging', 7))
    
    for (col in cols){
        SL[[col]] <- data.frame(row.names=rownames(combined_full_rep[[1]]))
        for (name in names(combined_full_rep)){
                tmp <- combined_full_rep[[name]][col]
                colnames(tmp) <- c(name)
                SL[[col]] <- transform(merge(SL[[col]], tmp, all=TRUE,
                                                    by=0, sort=F,), 
                                                    row.names=Row.names, Row.names=NULL)
        }
    }

    tmpdf <- combined_full_rep[[1]]
    gene_interest_etz <- rownames(tmpdf[tmpdf$symbol %in% gene_interest,])

    melter <- data.frame()
    for (g in gene_interest_etz){
            #func for mapping stars
            plab <- function(p )ifelse(p<0.001,'***', ifelse(p<0.01,'**', ifelse(p<0.05,'*',''))) 
            sym <- combined_full_rep[[1]][g,]$symbol

            #rescale columns
            tmpLogFC = SL[['logFC']] / apply(SL[['logFC']], 2, sd, na.rm=TRUE)
            tmpSE = SL[['SE']] / apply(SL[['logFC']], 2, sd, na.rm=TRUE)

            mel <- transform(merge(t(tmpLogFC[g,]), t(tmpSE[g,]), by=0, sort=F), row.names=Row.names, Row.names=NULL)
            mel <- merge(mel,  t(SL[['FDR']][g,]), by=0, sort=F)
            colnames(mel) <- c('Signature', 'Coef', 'SE', 'FDR')
            mel$Type <- types
            mel$Signature <- factor(as.character(mel$Signature), levels=unique(as.character(mel$Signature)))
            mel$Star <- plab(mel$FDR)
            mel$Vjust <- if_else(mel$Coef > 0, 0.3, 1.0)
            mel$Symbol <- rep(sym, nrow(mel))
            melter <- rbind(melter, mel)
    }
    melter$Symbol <- factor(melter$Symbol, levels=gene_interest)

    #plot
    #suppressMessages(library(ggrepel))
    g <- ggplot(melter, aes(x=Signature, y=Coef, fill=Type))+ 
            geom_bar(stat="identity", colour="black") +
            geom_hline(yintercept=0, linetype="solid", color = "black", size=0.4, alpha=0.8)+
            facet_wrap(~Symbol, ncol=1, scales = "free_y", dir="v") + 
            expand_limits(y=c(-3, 3))+
            geom_text(aes(label=Star, vjust=melter$Vjust, hjust=0.5), fontface='bold', size=7)+
            geom_errorbar(aes(ymin=Coef-SE, ymax=Coef+SE), width=0.18) + 
            theme(axis.text=element_text(size=fontsize, angle=45, hjust=1),
                    strip.background = element_blank(),
                    strip.placement = 'inside',
                    strip.text = element_text(colour = 'black', size=fontsize, face="bold", hjust = 0.5),     
                    axis.title=element_text(size=fontsize, face="bold"), 
                    plot.title=element_text(size=fontsize, face = "bold"),
                    legend.title=element_text(size=fontsize, face = "bold"),
                    legend.text=element_text(size=fontsize),
                    legend.position="top",
                    panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.margin = unit(c(0.5,0.5,0.5,0.), "cm")) +
            labs(x = "", y = "Normalized meta slope")
    return(g)
}                