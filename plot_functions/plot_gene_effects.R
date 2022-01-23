plot_gene_effects <- function(gef, top=10, thr=0.05, fontsize=18){
        agef <- aggregate(.~symbol, gef, FUN="mean")[,c("symbol", "knockout_effect", "FDR")]
        gef$star <- as.character(gef$star)
        agef['star']<- aggregate(.~symbol, gef, FUN=function(x){x[[1]]})[,c("symbol", "star")]$star
        agef <- agef[agef$FDR < 0.05,]

        selected_agef <- agef[order(agef$knockout_effect, decreasing=T),][1:top,]
        selected <- selected_agef$symbol
        selected_mean <- selected_agef$knockout_effect
        selected_star <- selected_agef$star
        sgef <- gef[gef$symbol %in% selected,]
        sgef$symbol <- factor(sgef$symbol, levels=selected)

        p <- ggplot(data=sgef, aes(x=symbol, y=knockout_effect)) +
                geom_boxplot(aes(fill=symbol), alpha = 0.8, outlier.shape = NA, lwd=1.2) + 
                #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="H")+
                #scale_fill_manual(values = cols)+
                scale_fill_brewer(palette='Purples')+
                geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.4, alpha=0.8)+
                geom_point(aes(fill='black'), alpha=0.4, size = 1., position = position_jitterdodge(0.1))+
                geom_text(data=data.frame(), aes(label=selected_star, x=selected, y=selected_mean), 
                        vjust=-4.0, fontface='bold', size=7, color='black')+
                theme(  axis.text=element_text(size=fontsize),
                        axis.title=element_text(size=fontsize, face="bold"),
                        panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.position = "none",
                        plot.margin = unit(c(1,0.5,0,2.0), "cm")) +
                        labs(x = "", y = "Rejuvenation effect decreasing, %")
        return(p)
}