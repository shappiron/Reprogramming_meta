plot_drug_effects <- function(drugHuman, fontsize = 18, top=14){
############
df <- drugHuman[!grepl("BRD", drugHuman$drug_name),]
df <- df[!duplicated(df[ , c("drug","OLS_intercept")]),]
df$drug <- as.character(df$drug)
df$drug <- gsub('oe', '', df$drug)
df$drug <- gsub('sh', '', df$drug)
df$drug <- gsub('xpr', '', df$drug)
df$drug <- gsub('cp', '', df$drug)
df$drug <- gsub('__', ' ', df$drug)
df$drug <- gsub('_', ' ', df$drug)
df$drug <- gsub('  ', ' ', df$drug)
df$drug <- factor(df$drug, levels = df[order(abs(df$OLS_pvalue), decreasing=T),]$drug)

df <- rbind(df[df$pert_type == 'cp',][1:top,], 
            df[df$pert_type == 'oe',][1:top,],
            df[df$pert_type == 'xpr',][1:top,],
            df[df$pert_type == 'sh',][1:top,])

labeller <- c(
  'cp'="Chemical compound",
  'oe'="Overexpression",
  'xpr'="CRISPR knock-out",
  'sh'="shRNA knock-out"
)

top <- ggplot(df, aes(x=OLS_intercept, y=drug))+#color=scores, shape=coincide) ) + 
                geom_point(aes(shape=coincide, fill=scores), colour="black", size=4) + 
                geom_errorbarh(aes(xmax = CI_high, xmin = CI_low),  height = 0.)+
                geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha=0.5)+
                facet_wrap(~pert_type, ncol=2, dir="v", scales='free_y', labeller=as_labeller(labeller)) +
                scale_shape_manual(values=c(23, 21))+
                scale_fill_gradient2(high="red2", low="blue2", mid='grey', 
                                limits=c(-4, 4), midpoint=0, oob = scales::squish)+
                theme(  axis.text=element_text(size=fontsize, color='black'),
                        axis.title=element_text(size=fontsize, face="bold"),
                        strip.background = element_blank(),
                        strip.placement = 'inside',
                        strip.text = element_text(colour = 'black', size=fontsize, face="bold", hjust = 0.5),
                        legend.title=element_text(size=fontsize, color='black'),
                        legend.text=element_text(size=fontsize, color='black'),
                        panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                        panel.grid.major = element_blank(),
                        #panel.grid.minor = element_blank(),
                        plot.margin = unit(c(1,1,1,0.5), "cm")) + 
                        labs(x = "Relative age change (Treatment coefficient in model)", 
                             y = "", shape='Coincide', fill='Significance \nscore')
  return(top)
##############   
}