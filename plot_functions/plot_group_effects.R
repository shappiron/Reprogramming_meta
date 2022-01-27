plot_group_effects <- function(gref, fontsize=18){
        library(tidyverse)
        t_tests = gref %>%
        group_by(variable) %>%
        summarise(p = t.test(value, mu = 0)$p.value,
                Sig = ifelse(p<0.001,'***', ifelse(p<0.01,'**', ifelse(p<0.05,'*',''))),
                MaxWidth = 75)

        p<-ggplot(gref, aes(x = variable, y = value)) +
                geom_boxplot(aes(fill=variable), size=1.2) +
                scale_fill_brewer(palette='Dark2')+
                geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.4, alpha=0.8)+
                geom_text(aes(label = Sig, y = MaxWidth), size = 12, data = t_tests)+
                theme_minimal()+
                theme(  axis.text=element_text(size=fontsize),
                        axis.title=element_text(size=fontsize, face="bold"),
                        panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.position = "none",
                        plot.margin = unit(c(1,1.5,0,0.5), "cm")) +
                labs(x = "", y = "Rejuvenation effect, %")
        return(p)
}