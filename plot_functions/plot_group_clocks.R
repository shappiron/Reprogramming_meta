##plot reprogramming tAGE clock dynamics grouped by Treatment/Tissue
plot_group_clocks <- function(clockm, clockh){
	#Data preparation
	plab <- function(p )ifelse(p<0.001,'***', ifelse(p<0.01,'**', ifelse(p<0.05,'*',''))) 
	##Mouse
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
	#add stars
	treats <- unique(clockm$Treatment)
	clockm$star <- NaN
	for (t in treats){
		tsub <- clockm[clockm$Treatment == t,]
		test <- lm(tAge_norm_mouseClock~Time, data=tsub)
		pval <- summary(test)$coefficients['Time', 'Pr(>|t|)']
		star <- plab(pval)
		clockm[clockm$Treatment == t,]$star <- star
	}
	clockm$TreatmentStar <- paste0(clockm$Treatment, ' ', clockm$star)

	##Human
	clockh$Pseudotime <- clockh$Time %/% 2 * 2
	clockh$Tissue <- factor(clockh$Tissue, levels=sort(as.vector(unique(clockh$Tissue)), decreasing=T))
	clockh <- clockh[order(clockh$Tissue),]
	clockh$ID <- factor(clockh$ID, levels=unique(clockh$ID))
	#add stars
	treats <- unique(clockh$Tissue)
	clockh$star <- NaN
	for (t in treats){
		tsub <- clockh[clockh$Tissue == t,]
		test <- lm(tAge_norm_globalclock~Time, data=tsub)
		pval <- summary(test)$coefficients['Time', 'Pr(>|t|)']
		star <- plab(pval)
		clockh[clockh$Tissue == t,]$star <- star
	}
	clockh$TissueStar <- paste0(clockh$Tissue, ' ', clockh$star)
    
    #plot
    dataYF <- clockm[clockm$RepFactors == 'YF',]
    g1 <- ggplot(dataYF, aes(x=Time, y=tAge_norm_mouseClock, col=TreatmentStar)) + 
                geom_point()+
                geom_smooth(method=lm, se=F)+
                ggtitle('Yamanaka Factors')+
                geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.4, alpha=0.4)+
                theme(  plot.title = element_text(size = 20, face = "bold"),
                        axis.text=element_text(size=14),
                        axis.title=element_text(size=16, face="bold"),
                        panel.background = element_rect(fill='white', colour='black', size=0.3, linetype='solid'),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.text=element_text(size=14),
                        legend.title=element_text(size=16, face='bold'),
                        plot.margin = unit(c(0,0,0,0.5), "cm"),
                        aspect.ratio=1) +
                        labs(x = "Time, days", y = "tAge", color='Treatment')
                        

    data7F <- clockm[clockm$RepFactors == '7F',]
    g2 <- ggplot(data7F, aes(x=Time, y=tAge_norm_mouseClock, col=TreatmentStar)) + 
                geom_point()+
                geom_smooth(method=lm, se=F)+
                ggtitle('7 Factors (GSE127927)')+
                geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.4, alpha=0.4)+
                theme(  plot.title = element_text(size = 20, face = "bold"),
                        axis.text=element_text(size=14),
                        axis.title=element_text(size=16, face="bold"),
                        panel.background = element_rect(fill='white', colour='black', size=0.3, linetype='solid'),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.text=element_text(size=14),
                        legend.title=element_text(size=16, face='bold'),
                        plot.margin = unit(c(0,0,0,0.5), "cm"),
                        aspect.ratio=1) +
                        labs(x = "Time, days", y = "tAge", color='Treatment')


    dataH <- clockh
    g3 <- ggplot(dataH, aes(x=Time, y=tAge_norm_globalclock, col=TissueStar)) + 
                geom_point()+
                geom_smooth(method=lm, se=F)+
                ggtitle('Human OSKM')+
                geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.4, alpha=0.4)+
                theme(  plot.title = element_text(size = 20, face = "bold"),
                        axis.text=element_text(size=14),
                        axis.title=element_text(size=16, face="bold"),
                        panel.background = element_rect(fill='white', colour='black', size=0.3, linetype='solid'),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.text=element_text(size=14),
                        legend.title=element_text(size=16, face='bold'),
                        plot.margin = unit(c(0,0,0,0.5), "cm"),
                        aspect.ratio=1) +
                        labs(x = "Time, days", y = "tAge", color='Tissue')
                    

    G <- cowplot::plot_grid(g1, g2, g3, labels=c("C", "", 'D'), label_size = 32, ncol=1, align='v', 
                            rel_heights = c(3, 3, 3))
    return(G)
}