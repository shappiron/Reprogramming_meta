plot_clock_effects <- function(df){
    names <- colnames(df)[3:6]
    top <- 4

    effect_table <- data.frame()
    for (n in names){
        effect <- df[n] * df['w']
        rej_effect <- effect[effect<0,, drop=FALSE]
        age_effect <- effect[effect>0,, drop=FALSE]
        rej_effect <- rej_effect / sum(rej_effect)
        age_effect <- age_effect / sum(age_effect)
        rej_effect <- rej_effect[order(-rej_effect[[n]]),, drop=F][1:top,, drop=F]
        age_effect <- age_effect[order(-age_effect[[n]]),, drop=F][1:top,, drop=F]
        effect <- rbind(age_effect, rej_effect)
        effect$effect_name <- rep(c('Aging', 'Rejuvenation'), c(top, top))
        colnames(effect) <- c('effect', 'effect_name')
        effect$symbol <- as.character(df[rownames(effect), 'symbol'])
        effect$dataset <- rep(n, nrow(effect))
        effect$entrez <- rownames(effect)
        rownames(effect) <- NULL
        effect_table <- rbind(effect_table, effect)
    }
    effect_table$effect <- effect_table$effect * 100 #make percents
    effect_table$symbol <- factor(effect_table$symbol, levels=unique(effect_table$symbol))
    effect_table$symbol <- reorder(effect_table$symbol, effect_table$symbol)

    g <- ggplot(effect_table, aes(x=symbol, y=effect, fill=dataset))+ 
            geom_bar(stat="identity") +
            facet_wrap(~effect_name, nrow=1, scales = "free_x", dir="v") + 
            expand_limits(y=c(-3, 3))+
            theme(axis.text=element_text(size=16, angle=45, hjust=1),
                    strip.background = element_blank(),
                    strip.placement = 'inside',
                    strip.text = element_text(colour = 'black', size=16, face="bold", hjust = 0.5),     
                    axis.title=element_text(size=16, face="bold"), 
                    legend.title=element_text(size=16, face = "bold"),
                    legend.text=element_text(size=16),
                    legend.position="top",
                    panel.background = element_rect(fill='white', colour='black', size=1.1, linetype='solid'),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.margin = unit(c(0.5,2.,0.,2.), "cm")) +
            labs(x = "", y = "Effect size, %", fill='Dataset')
    return(g)
}