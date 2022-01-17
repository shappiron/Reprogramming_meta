plot_contingency_pair <- function(X, Y, suf1='Mouse', suf2='Human', thr=0.05, lims=c(), textcol='white'){
    X <- X[X$FDR < thr,]
    Y <- Y[Y$FDR < thr,]

    sub <- merge(X, Y, by=0)
    col1 <- "logFC.x"
    col2 <- "logFC.y"
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
            geom_tile(colour = "black", size=0.5) + 
            scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) + 
            geom_text(data=data, aes(label=value), size=10, col=textcol) + 
            labs(tag = paste0("pval=", pval), x=suf1, y=suf2) +
            scale_fill_gradient2(high="red2", low="blue2", mid='white', limits=lims, midpoint=avgdata)+
            guides(fill=guide_colorbar(title.position = "top", title='#Observed - #Expected')) + 
            theme(legend.position = "top", legend.key.width = unit(1.8, "cm"),
                    axis.text=element_text(size=16), axis.title=element_text(size=16, face="bold"),
                    legend.text=element_text(size=12), 
                    legend.title=element_text(size=14, face="bold"), 
                    legend.title.align=0.5,
                    plot.tag.position = c(0.15, 0.05), plot.tag = element_text(size=16),
                    plot.margin = unit(c(0.,5.0, 1.5, 5.0), "cm")) 
    return(g)
}    