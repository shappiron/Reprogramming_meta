#' Generate a ggplot2 heatmap with row and column dendrograms
#'
#' @param dataMatrix A data.frame containing the input data.
#' @param orderCol Reorder the columns (default=T)
#' @param orderRow Reorder the rows (default=T)
#' @param dendroLineSize Size of the dendrogram lines (default=0.5)
#' @param fontsize Font size (default=20)
#' @param colorPalette Color palette (default='Spectral')
#' @param revColors Invert color scale
#' @param scaleName Name of the colorscale (default='value')
#' @param distMethod Distance method (default='euclidean', see ?dist)
#' @param clustMethod Clustering method (default='complete', see ?hclust)
#' @examples ggheatmap(mtcars)
#' @importFrom magrittr %>%
#' @export 
ggheatmap <- function(data, 
    type='full_intersection', method='spearman', top=500, pvplot=F,
    label_size=1.2, value_size=0.8, title="Spearman correlation",
    measure="logFC", criterion="FDR",
    orderCol = T, orderRow = T, dendroLineSize = 0.5, 
    fontsize = 20, color = "Spectral", scaleName = "value", distMethod = "euclidean", 
    clustMethod = "complete", revColors=F, scale_name="Spearman \ncorrelation ") {
    
    names(data) <- gsub("\\$", ":", names(data))
    names(data) <- gsub("_", " ", names(data))
    names(data) <- gsub("\\+dox_mef", "", names(data))
    names(data) <- gsub(" liver", "", names(data))
    names(data) <- gsub("All", "Global", names(data))

    #type = {full_intersection, pairwise_top}
    if (type == 'full_intersection'){
        #provided data must be a data.frame type
        cm <- round(cor(data, method='spearman'), 2)
        cellnote <- cm
    }else if (type == 'pairwise_top'){
        #pairwise intersection of top N genes or full in case of top==F
        #provided data must be a list of data frames
        cm <- matrix(0:0, nrow=length(data), ncol=length(data))
        pm <- matrix(0:0, nrow=length(data), ncol=length(data))
        for (i in 1:length(data)){
            for (j in 1:length(data)){
                if (top > 0) {
                    #intersection of all genes
                    I <- intersect(rownames(data[[i]]), rownames(data[[j]]))
                    d1 <- data[[i]][I,]
                    d2 <- data[[j]][I,]
                    #union of top genes
                    nx <- min(c(top, nrow(d1)))
                    ny <- min(c(top, nrow(d2)))
                    rx <- rownames(d1[order(d1[[criterion]]),][1:nx,])
                    ry <- rownames(d2[order(d2[[criterion]]),][1:ny,])
                    r <- union(rx, ry)
                    x <- d1[r,][measure]
                    y <- d2[r,][measure]
                }else{
                    x <- data[[i]][measure]
                    y <- data[[j]][measure]   
                }
                pair <- transform(merge(x, y, by="row.names", sort=FALSE), row.names=Row.names, Row.names=NULL)
                #pair <- pair[complete.cases(pair),]
                if (nrow(pair) > 1){
                    test <- cor.test(pair[[1]], pair[[2]], method=method)
                    sp <- test$estimate
                    pv <- test$p.value
                }else{
                    sp <- NA
                }
                if (is.na(sp)){
                    cm[i, j] <- 0
                    pm[i, j] <- 1.0
                }else{
                    cm[i, j] <- sp
                    pm[i, j] <- pv
                }    
            }
        }
        cm <- round(cm, 2)
        rownames(cm) <- names(data) 
        colnames(cm) <- names(data)
        rownames(pm) <- names(data) 
        colnames(pm) <- names(data)
        #prepare value annotation for heatmap
        FDR_table <- pm
        #extract upper triangular of the matrix and apply FDR correction
        tmp <- matrix(0, nrow=nrow(FDR_table), ncol=ncol(FDR_table))
        padj <- matrix(p.adjust(FDR_table[upper.tri(FDR_table)], method="BH"))
        tmp[upper.tri(tmp)] <- padj
        tmp <- tmp + t(tmp)
        diag(tmp) <- 0
        FDR_table <- as.data.frame(tmp, nrow=nrow(corr_matrix_list$P.Value))

        cellnote <- cm
        cellnote[FDR_table >= 0.05] <- NA
    }else{
        print("Please, choose a type")
    }
    dataMatrix <- as.data.frame(cm)
    pvalMatrix <- as.data.frame(cellnote)
    data_m <- tibble::rownames_to_column(dataMatrix) %>% reshape2::melt()
    data_p <- tibble::rownames_to_column(pvalMatrix) %>% reshape2::melt()
    data_m <- merge(data_m, data_p, by=c('rowname', 'variable'), suffixes=c('', '.pass'))

    # Cluster rows
    if (orderRow) {
        dd.row <- as.dendrogram(hclust(dist(dataMatrix, method = distMethod), method = clustMethod))
        row.ord <- order.dendrogram(dd.row)
        ordered_row_names <- row.names(dataMatrix[row.ord, ])
        data_m$rowname <- factor(data_m$rowname, levels = ordered_row_names)
    }
    
    # Cluster columns
    if (orderCol) {
        dd.col <- as.dendrogram(hclust(dist(t(dataMatrix), method = distMethod), 
                      method = clustMethod))
        col.ord <- order.dendrogram(dd.col)
        ordered_col_names <- colnames(dataMatrix[, col.ord])
        data_m$variable <- factor(data_m$variable, levels = ordered_col_names)
    }
    
    heat_plot <- ggplot(data_m, aes(x = variable, y = rowname, fill = value)) + 
                    geom_tile(colour = "black", size=0.) + 
                    #geom_text(aes(label=value.pass)) +
                    theme_minimal() + 
                    theme(
                        axis.line = element_line(size = 0),
                        text = element_text(size = fontsize, color='black'),
                        axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, color='black'),
                        axis.text.y = element_text(color='black'),
                        legend.position="top", 
                        legend.key.width = unit(2.5, "cm"),
                        legend.title = element_text(size = fontsize, face='bold'),
                        plot.margin = unit(c(1, 0, 0, 1.5), "cm"),
                        panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        ) + 
                    #guides(fill=guide_legend(title="Spearman correlation", size=fontsize)) +
                    scale_y_discrete(position = "right") + 
                    xlab("") + 
                    ylab("") + 
                    scale_fill_gradient2(name=scale_name, 
                            high="purple3", low="orange3", mid='white', 
                            limits=c(-1, 1), 
                            midpoint=0)
    
    final_plot <- heat_plot

    if (orderRow) {
        dendro_data_row <- ggdendro::dendro_data(dd.row, type = "rectangle")
        dendro_row <- cowplot::axis_canvas(heat_plot, axis = "y", coord_flip = TRUE) + 
                          geom_segment(data = ggdendro::segment(dendro_data_row), aes(y = -y, 
	                  x = x, xend = xend, yend = -yend), size = dendroLineSize) + 
                          coord_flip()
        final_plot <- cowplot::insert_yaxis_grob(final_plot, dendro_row, grid::unit(0.2, 
                          "null"), position = "left")
    }
    
    # if (orderCol) {
    #     dendro_data_col <- ggdendro::dendro_data(dd.col, type = "rectangle")
    #     dendro_col <- cowplot::axis_canvas(heat_plot, axis = "x") + 
	# 	      geom_segment(data = ggdendro::segment(dendro_data_col), aes(x = x, y = y, xend = xend, yend = yend), size = dendroLineSize)
    #     final_plot <- cowplot::insert_xaxis_grob(final_plot, dendro_col, grid::unit(0.2, 
    #                       "null"), position = "top")
    # }
    
    return(final_plot)
    
}