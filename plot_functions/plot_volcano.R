#Build volcano plot for a signature with pointing pluripotency genes
plot_volcano <- function(db, pluripotency, fdr_thr=0.01, coef_thr=2.0,
                        xlim=c(-6, 6)){
    library(ggrepel)
    db$diffexpressed <- "NO"
    db$diffexpressed[(db$logFC > coef_thr) & (db$FDR < fdr_thr)] <- "UP"
    db$diffexpressed[(db$logFC < -coef_thr) & (db$FDR < fdr_thr)] <- "DOWN"

    db$delabel <- NA
    db$delabel[db$diffexpressed != "NO"] <- as.character(db$symbol[db$diffexpressed != "NO"])

    db[db$diffexpressed != "NO",][pluripotency,]$diffexpressed <- 'PLURI'

    mycolors <- c("blue", "red", "black", 'green')
    names(mycolors) <- c("DOWN", "UP", "NO", "PLURI")

    g <- ggplot(data=db, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=delabel)) + 
                geom_point(size=1) + theme_minimal() + geom_text_repel(size=4) +
                geom_vline(xintercept=c(-coef_thr, coef_thr), col="red", size=0.5) +
                geom_hline(yintercept=-log10(fdr_thr), col="red", size=0.5) + 
                xlim(xlim[[1]], xlim[[2]]) + 
                scale_colour_manual(values = mycolors) + 
                theme(plot.margin = unit(c(0,0,1,0), "cm"),
                      panel.background = element_rect(fill='white', colour='black', size=1.0, linetype='solid'),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.text=element_text(size=14),
                      axis.title=element_text(size=16, face="bold"),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),) +
                labs(x = "Meta-slope value", y = "-log10(adj.P-value)") +
                guides(colour=guide_legend(title="Type"))
    return(g)
}