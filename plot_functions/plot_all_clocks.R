#plot all reprogramming datasets age dynamics separately
plot_all_clocks <- function(clockm, clockh, fontsize=18){
    #Data preparation
    clockm$BatchTreatment <- factor(clockm$BatchTreatment, levels=c('7F', 'OSKM', 'OSK', "Control"))
    clockm[clockm$Treatment == 'OSKM+dox_mef',]$Treatment <- "OSKM"
    clockm[clockm$Treatment == 'OSKM_WT-1',]$Treatment <- "OSKM"
    clockm['ID'] = paste(clockm$GEO, ':', clockm$Treatment, sep='')
    #clockm$Pseudotime <- clockm$Time %/% 2 * 2
    clockm$Treatment <- factor(clockm$Treatment, levels=sort(as.vector(unique(clockm$Treatment)), decreasing=T))
    clockm <- clockm[order(clockm$Treatment),]
    clockm$ID <- factor(clockm$ID, levels=unique(clockm$ID))
    clockm$Dataset <- gsub("α\\+", "a\\+", clockm$Dataset)
    clockm$Dataset <- gsub("α-", "a-", clockm$Dataset)
    clockm$Dataset <- gsub("dox_mef", "dox", clockm$Dataset)

    #clockh$Pseudotime <- clockh$Time %/% 2 * 2
    clockh$Tissue <- factor(clockh$Tissue, levels=sort(as.vector(unique(clockh$Tissue)), decreasing=T))
    clockh <- clockh[order(clockh$Tissue),]
    clockh$ID <- factor(clockh$ID, levels=unique(clockh$ID))

    #plot
    #clockm
    gA <- ggplot(data=clockm, aes(x=Time, y=AgeNorm)) + 
            geom_point(aes(fill=BatchTreatment), size=2.5, shape = 21, colour = "black")+
            #geom_point(shape = 1,size = 2.,colour = "black")+
            geom_hline(yintercept=0, linetype="solid", color = "black", size=0.2, alpha=0.2)+
            facet_wrap(~Dataset, ncol=6, dir="h") +
            theme(  axis.text=element_text(size=fontsize),
                axis.title=element_text(size=fontsize, face="bold"),
                legend.title=element_text(size=fontsize, face="bold"),
                legend.text=element_text(size=fontsize, ),
                strip.background = element_blank(),
                strip.placement = 'inside',
                strip.text = element_text(colour = 'black', size=fontsize-5, face="bold", hjust = 0.5),
                panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "top",
                plot.margin = unit(c(1,1,1,0.5), "cm")) +
                labs(x = "Time, days", y = "tAge", color='Treatment')

    #clockh
    gB <- ggplot(data=clockh, aes(x=Time, y=AgeNorm, col=Tissue)) + 
            geom_point(aes(fill=Tissue), size=2.5, shape = 21, colour = "black")+
            #geom_point(shape = 1,size = 2.,colour = "black")+
            geom_hline(yintercept=0, linetype="solid", color = "black", size=0.2, alpha=0.2)+
            facet_wrap(~Dataset, ncol=6, dir="h") +
            theme(  axis.text=element_text(size=fontsize),
                axis.title=element_text(size=fontsize, face="bold"),
                legend.title=element_text(size=fontsize, face="bold"),
                legend.text=element_text(size=fontsize, ),
                strip.background = element_blank(),
                strip.placement = 'inside',
                strip.text = element_text(colour = 'black', size=fontsize-5, face="bold", hjust = 0.5),
                panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "top",
                plot.margin = unit(c(1,1,1,0.5), "cm")) +
                labs(x = "Time, days", y = "tAge", color='Tissue') + 
                guides(color = guide_legend(nrow = 1))

    P <- cowplot::plot_grid(gA, gB, labels=c("A", "B"), label_size = 32, ncol=1, 
                            rel_heights = c(4, 2))
    return(P)
}
