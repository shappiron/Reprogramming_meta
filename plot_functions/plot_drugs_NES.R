plot_drugs_NES <- function(sd, dh, fontsize=18){
        dh <- dh[!duplicated(dh[ , c("drug", "OLS_intercept")]),]
        df <- merge(sd, dh, by.x=c('gene', 'ptime'), by.y=c('drug_name', 'pert_time'))
        df$Treatment <- paste(df$gene, df$ptime, df$pert_type)
        library(ggrepel)

        bottom <- ggplot(df, aes(x=OLS_intercept, y=avgNES))+
                geom_point(aes(fill=Treatment), colour="black", size=8, shape=21) + 
                #geom_text_repel(size=6) +
                geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha=0.5)+
                geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha=0.5)+
                geom_errorbarh(aes(xmax = CI_high, xmin = CI_low),  height = 0., size=0.3)+
                geom_errorbar(aes(ymax = CIH, ymin = CIL),  width = 0., size=0.3)+
                theme(  axis.text=element_text(size=fontsize, color='black'),
                        axis.title=element_text(size=fontsize, face="bold"),
                        legend.title=element_text(size=fontsize, color='black', face="bold"),
                        legend.text=element_text(size=fontsize, color='black'),
                        #legend.position='none',
                        panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                        panel.grid.major = element_blank(),
                        #panel.grid.minor = element_blank(),
                        plot.margin = unit(c(1,1,1,2.5), "cm")) + 
                labs(x = "Relative age change (Treatment coefficient in model)", 
                        y = "Pluripotency induction \n(Average NES)", fill='Treatment')
        return(bottom)
}
#######