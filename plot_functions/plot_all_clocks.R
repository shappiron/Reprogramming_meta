#plot all reprogramming datasets age dynamics separately
plot_all_clocks <- function(clockm, clockh){
    #Data preparation
    clockm$RepFactors <- ifelse(clockm$BatchTreatment == 'OSKM', 'YF', '7F')
    clockm$BatchTreatment <- factor(clockm$BatchTreatment, levels=c('7F', 'OSKM', 'OSK'))
    clockm[clockm$Treatment == 'OSKM+dox_mef',]$Treatment <- "OSKM"
    clockm[clockm$Treatment == 'OSKM.WT-1',]$Treatment <- "OSKM"
    clockm[clockm$Treatment == 'OSK',]$BatchTreatment <- "OSK"
    clockm['ID'] = paste(clockm$GEO, ':', clockm$Treatment, sep='')
    clockm$Pseudotime <- clockm$Time %/% 2 * 2
    clockm$Treatment <- factor(clockm$Treatment, levels=sort(as.vector(unique(clockm$Treatment)), decreasing=T))
    clockm <- clockm[order(clockm$Treatment),]
    clockm$ID <- factor(clockm$ID, levels=unique(clockm$ID))

    clockh$Pseudotime <- clockh$Time %/% 2 * 2
    clockh$Tissue <- factor(clockh$Tissue, levels=sort(as.vector(unique(clockh$Tissue)), decreasing=T))
    clockh <- clockh[order(clockh$Tissue),]
    clockh$ID <- factor(clockh$ID, levels=unique(clockh$ID))

    #plot
    #clockm
    options(repr.plot.width = 18, repr.plot.height = 15)
    gA <- ggplot(data=clockm, aes(x=Time, y=tAge_norm_mouseClock, col=BatchTreatment)) + 
            geom_point()+
            geom_hline(yintercept=0, linetype="solid", color = "black", size=0.2, alpha=0.2)+
            facet_wrap(~ID, nrow=3, dir="h") +
            theme(  axis.text=element_text(size=14),
                axis.title=element_text(size=14, face="bold"),
                legend.title=element_text(size=14, face="bold"),
                legend.text=element_text(size=12, ),
                strip.background = element_blank(),
                strip.placement = 'inside',
                strip.text = element_text(colour = 'black', size=11, face="bold", hjust = 0.5),
                panel.background = element_rect(fill='white', colour='black', size=0.3, linetype='solid'),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "top",
                plot.margin = unit(c(1,1,1,0.5), "cm")) +
                labs(x = "Time, days", y = "tAge", color='Treatment')

    #clockh
    gB <- ggplot(data=clockh, aes(x=Time, y=tAge_norm_globalclock, col=Tissue)) + 
            geom_point()+
            geom_hline(yintercept=0, linetype="solid", color = "black", size=0.2, alpha=0.2)+
            facet_wrap(~ID, nrow=2, dir="h") +
            theme(  axis.text=element_text(size=14),
                axis.title=element_text(size=14, face="bold"),
                legend.title=element_text(size=14, face="bold"),
                legend.text=element_text(size=12, ),
                strip.background = element_blank(),
                strip.placement = 'inside',
                strip.text = element_text(colour = 'black', size=14, face="bold", hjust = 0.5),
                panel.background = element_rect(fill='white', colour='black', size=0.3, linetype='solid'),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "top",
                plot.margin = unit(c(1,1,1,0.5), "cm")) +
                labs(x = "Time, days", y = "tAge", color='Tissue')

    cowplot::plot_grid(gA, gB, labels=c("A", "B"), label_size = 32, ncol=1, 
                            rel_heights = c(3, 2.7))
}
