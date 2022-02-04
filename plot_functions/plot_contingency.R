plot_contingency <- function(combined_full_rep, suf1='Interventions', suf2='Reprogramming', 
                             lims=c(), textcol='white', fontsize=18){
    #combine list of signatures into list of dataframes where each dataframe is the same column type
    SLR <- list()
    SLI <- list()
    SLA <- list()
    cols <- c('logFC', 'SE', 'P.Value', 'FDR')

    for (col in cols){
        SLR[[col]] <- data.frame(row.names=rownames(combined_full_rep[['Reprogramming:Mouse']]))
        SLI[[col]] <- data.frame(row.names=rownames(combined_full_rep[['Interventions:Median_lifespan']]))
        SLA[[col]] <- data.frame(row.names=rownames(combined_full_rep[['Aging:Brain']]))
        for (name in names(combined_full_rep)){
            tmp <- combined_full_rep[[name]][col]
            colnames(tmp) <- c(name)
            if (grepl("Aging", name, fixed=T)){    
                SLA[[col]] <- transform(merge(SLA[[col]], tmp, 
                                                    by=0, sort=F,), 
                                                    row.names=Row.names, Row.names=NULL)
            }else if (grepl("Interventions", name, fixed=T)){
                SLI[[col]] <- transform(merge(SLI[[col]], tmp, 
                                                    by=0, sort=F,), 
                                                    row.names=Row.names, Row.names=NULL)
            }else{
                SLR[[col]] <- transform(merge(SLR[[col]], tmp, 
                                        by=0, sort=F,), 
                                        row.names=Row.names, Row.names=NULL) 
            }
        }
    }
    agg_func <- function(x) 1 / mean(1 / x)    #harmonic mean
    #agg_func <- function(x) exp(mean(log(x))) #geometric mean
    slr <- data.frame(cbind(rowMeans(SLR[['logFC']]), apply(SLR[['P.Value']], 1, agg_func)))
    colnames(slr) <- c('logFC', 'P.Value')
    slr$FDR <- p.adjust(slr$P.Value, method="BH")

    sli <- data.frame(cbind(rowMeans(SLI[['logFC']]), apply(SLI[['P.Value']], 1, agg_func)))
    colnames(sli) <- c('logFC', 'P.Value')
    sli$FDR <- p.adjust(sli$P.Value, method="BH")

    sla <- data.frame(cbind(rowMeans(SLA[['logFC']]), apply(SLA[['P.Value']], 1, agg_func)))
    colnames(sla) <- c('logFC', 'P.Value')
    sla$FDR <- p.adjust(sla$P.Value, method="BH")

    subri = merge(slr, sli, by=0, suffixes=c(".Reprogramming", ".Interventions"))
    subra = merge(slr, sla, by=0, suffixes=c(".Reprogramming", ".Aging"))

    if (suf1 == 'Interventions'){
        sub <- subri[(subri$FDR.Interventions < 0.05) & (subri$FDR.Reprogramming < 0.05),]
    }else{
        sub <- subra[(subra$FDR.Aging < 0.05) & (subra$FDR.Reprogramming < 0.05),]
    }

    col1 <- paste0('logFC', '.', suf1)
    col2 <- paste0('logFC', '.', suf2)

    urua <- length(sub[(sub[[col1]] > 0) & (sub[[col2]] > 0),]$Row.names)
    urda <- length(sub[(sub[[col1]] < 0) & (sub[[col2]] > 0),]$Row.names)
    drua <- length(sub[(sub[[col1]] > 0) & (sub[[col2]] < 0),]$Row.names)
    drda <- length(sub[(sub[[col1]] < 0) & (sub[[col2]] < 0),]$Row.names)


    if (all(c(urua, urda, drua, drda)) == T){
        #test for distribution randomness
        cont_table <- cbind(c(urua, drua), 
                            c(urda, drda))
        test <- chisq.test(cont_table)
        pval <- test$p.val
    }

    data <- as.matrix(test$observed)
    #ex <- as.matrix(test$exprected)
    rownames(data) <- c('Up', 'Down')
    colnames(data) <- c('Up', 'Down')

    suppressMessages(library(reshape2))
    data <- melt(data)
    colnames(data) <- c(suf1, suf2, 'value')
    data$exp <- melt(test$expected)$value
    data$OE <- data$value - data$exp
    mindata <- min(data$OE)
    maxdata <- max(data$OE)
    avgdata <- mean(c(mindata, maxdata)) #exp(mean(log(data$value))) 
    pval = format.pval(test$p.value, digits = 3)

    if (length(lims) == 0){
    lims = c(mindata, maxdata)
    }

    library(reshape2)
    library(scales)
#    options(repr.plot.width = 5, repr.plot.height = 5)
    g <- ggplot(data, aes(x=data[[suf1]], y=data[[suf2]], fill=OE)) + 
            geom_tile(colour = "black", size=0.0) + 
            scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) + 
            geom_text(data=data, aes(label=value), size=8, col=textcol) + 
            labs(tag = paste0("P-value=", pval), x=suf1, y=suf2) +
            scale_fill_gradient2(high="red2", low="blue2", mid='white', limits=lims, midpoint=avgdata)+
            guides(fill=guide_colorbar(title.position = "top", title='#Observed - #Expected')) + 
            theme(legend.position = "top", 
                    legend.key.width = unit(1.5, "cm"),
                    axis.text=element_text(size=fontsize), 
                    axis.title=element_text(size=fontsize, face="bold"),
                    legend.text=element_text(size=fontsize), 
                    legend.title=element_text(size=fontsize, face="bold"), 
                    legend.title.align=0.5,
                    panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                    plot.tag.position = c(0.02, 0.02), 
                    plot.tag = element_text(size=fontsize),
                    plot.margin = unit(c(0.,2.0, 0.5, 2.0), "cm")) 
    return(g)
}    