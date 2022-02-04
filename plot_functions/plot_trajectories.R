plot_trajectories <- function(ma_outputs, fontsize=18){
    limits <- aes(ymin=Mean-SE,ymax=Mean+SE)
    t <- ggplot(ma_outputs,aes(x=Time,y=Mean,fill=Species,color=Species))+
            theme_bw(base_size = 16)+
            scale_fill_brewer(palette='Dark2')+
            scale_color_brewer(palette='Dark2')+
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.4, alpha=0.4)+
            geom_line(aes(color=Species), lwd=1)+
            geom_point(colour="black",size=4,shape=21,alpha=0.9)+
            geom_errorbar(limits,width=0.5,size=0.75)+
            theme(plot.title = element_text(size=fontsize, hjust = 0.,face="bold", color="black"),
                    axis.text.x = element_text(angle = 0, hjust = 0.5,size=fontsize,color="black"),
                    axis.text.y = element_text(size=fontsize,color="black"),
                    axis.title=element_text(size=fontsize, face='bold'),
                    legend.title=element_text(size=fontsize,color="black", face='bold'),
                    legend.text=element_text(size=fontsize,color="black"),
                    legend.key.size = unit(1.5, 'lines'),
                    legend.position = 'bottom',
                    panel.grid = element_blank(),
                    panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                    strip.background =element_rect(fill="white"),
                    strip.text = element_text(size=12,face="bold"))+
            labs(x="Time, days",y="tAge, normalized")
return(t)        
}