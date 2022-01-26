#Build volcano plot for a signature with pointing pluripotency genes
plot_volcano <- function(db, pluripotency, fdr_thr=0.01, 
                        coef_thr=c(-2.0, 2.0),
                        xlim=c(-6, 6), fontsize=18){
      left_thr <- coef_thr[[1]]
      right_thr <- coef_thr[[2]]
      library(ggrepel)
      db$diffexpressed <- "No"
      db$diffexpressed[(db$logFC > right_thr) & (db$FDR < fdr_thr)] <- "Up"
      db$diffexpressed[(db$logFC < left_thr) & (db$FDR < fdr_thr)] <- "Down"

      db$delabel <- NA
      db$delabel[db$diffexpressed != "No"] <- as.character(db$symbol[db$diffexpressed != "No"])

      db[db$diffexpressed != "No",][pluripotency,]$diffexpressed <- 'Pluripotent'

      mycolors <- c("blue", "red", "black", 'darkgreen')
      names(mycolors) <- c("Down", "Up", "No", "Pluripotent")

      g <- ggplot(data=db, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=delabel)) + 
                  geom_point(size=1) + 
                  theme_minimal() + 
                  geom_text_repel(size=5) +
                  geom_vline(xintercept=c(left_thr, right_thr), col="black", size=0.5, linestyle='dashed') +
                  geom_hline(yintercept=-log10(fdr_thr), col="black", size=0.5, linestyle='dashed') + 
                  xlim(xlim[[1]], xlim[[2]]) + 
                  scale_colour_manual(breaks=c("Down", "Up", "", "Pluripotent"), values = mycolors) + 
                  theme(plot.margin = unit(c(0,1,0,1), "cm"),
                        panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.text=element_text(size=fontsize),
                        axis.title=element_text(size=fontsize, face="bold"),
                        legend.text=element_text(size=fontsize),
                        legend.title=element_text(size=fontsize, face="bold"),
                        legend.position='top') +
                  labs(x = "Meta-slope value", y = "-log10(adj.P-value)") +
                  guides(colour=guide_legend(title=""))
    return(g)
}