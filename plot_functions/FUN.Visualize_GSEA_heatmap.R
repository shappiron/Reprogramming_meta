#This function visualizes GSEA output for selected functions in a heatmap format

Visualize_GSEA_heatmap <- function(GSEA_output,
                                   selected_functions,
                                   Gene_ratio_table=F,
                                   sort_mean=T,
                                   scale=F,max_value=F,min_value=F,
                                   order_rows=T,order_cols=T,
                                   clustering_method="spearman",
                                   linkage="complete",color_two_sided=T,
                                   background_color="white",
                                   max_size=7,
                                   as.circles=F,flip=F,
                                   labels_nogroup=T,
                                   color_labels_datasets=F){

  temp_output <- GSEA_output[selected_functions,]
  if (class(Gene_ratio_table)!="logical"){
    temp_ratio_output <- Gene_ratio_table[selected_functions,]
    temp_ratio_output[temp_output==0] <- 0
    temp_ratio_output <- as.data.frame(temp_ratio_output)
  }

  if (labels_nogroup==T){
    colnames(temp_output) <- sapply(colnames(temp_output),
                                    function(x){strsplit(x,":")[[1]][2]})
  }
  
  
    
  library(scales)
  if (scale==T){
    temp_output <- t(scale(t(temp_output),center = F,
                         scale = apply(t(temp_output), 2, sd, na.rm = TRUE)))
  }
  
  temp_output <- as.data.frame(temp_output)
  
  #Clustering
  if (order_rows==T){
    if (clustering_method=="spearman"){
      order <- hclust(as.dist((1-cor(t(temp_output),
                                     method="spearman",use = "complete.obs"))/2), 
                      method = linkage )$order
      
    }else if(clustering_method=="pearson"){
      order <- hclust(as.dist((1-cor(t(temp_output),
                                     method="pearson",use = "complete.obs"))/2), 
                      method = linkage )$order
    }else if(clustering_method=="euclidean"){
      order <- hclust(dist(temp_output,method="euclidean"), 
                      method = linkage)$order
    }else if(clustering_method=="manhattan"){
      order <- hclust(dist(temp_output,method="manhattan"), 
                      method = linkage)$order
    }else{
      stop("Invalid clustering method")
    }
    temp_output <- temp_output[order,]
    if (class(Gene_ratio_table)!="logical"){
      temp_ratio_output <- temp_ratio_output[order,]
    }
    
  }
  if (order_cols==T){
    if (clustering_method=="spearman"){
      order <- hclust(as.dist((1-cor(temp_output,method="spearman",use="complete.obs"))/2), 
                    method = linkage)$order
    }else if(clustering_method=="pearson"){
      order <- hclust(as.dist((1-cor(temp_output,method="pearson",use="complete.obs"))/2), 
                      method = linkage)$order
    }else{
      order <- hclust(dist(t(temp_output),method="euclidean"), 
                      method = linkage)$order
    }
    temp_output <- temp_output[,order]
    if (class(Gene_ratio_table)!="logical"){
      temp_ratio_output <- temp_ratio_output[,order]
    }
    
    if (color_labels_datasets!=F){
      color_labels_datasets <- color_labels_datasets[order]
    }
    
    
  }
  
  if (sort_mean==T){
    new_order <- order(rowMeans(temp_output),decreasing = T)
    temp_output <- temp_output[new_order,]

    if (class(Gene_ratio_table)!="logical"){
      temp_ratio_output <- temp_ratio_output[new_order,]
    }
    if (color_labels_datasets!=F){
      color_labels_datasets <- color_labels_datasets[new_order]
    }
    
    
  }
  #order2 <- hclust(dist(t(temp_output), method = "euclidean"), 
  #                method = "complete" )$order
  #temp_output <- temp_output[,order2]
  
    
  #Updating names
  temp_output$Term <- tolower(rownames(temp_output))
  temp_output$Term <- gsub("_"," ",rownames(temp_output))
  temp_output$Term <- gsub("[:]"," ",temp_output$Term)
  temp_output$Term <- paste(sapply(temp_output$Term,function(x) {strsplit(x," ")[[1]][1]}),
                            sapply(temp_output$Term,function(x) {paste(toupper(substr(paste(strsplit(x," ")[[1]][-1],collapse =" "), 1, 1)), 
                                                                       tolower(substr(paste(strsplit(x," ")[[1]][-1],collapse =" "), 2, nchar(x))), sep="")}),
                            sep=" ")
  temp_output$Term <- gsub("p450","P450",temp_output$Term)
  temp_output$Term <- gsub("Interleukin","IL",temp_output$Term)
  temp_output$Term <- gsub("interleukin","IL",temp_output$Term)
  temp_output$Term <- gsub("sars","SARS",temp_output$Term)
  temp_output$Term <- gsub("Sars","SARS",temp_output$Term)
  temp_output$Term <- gsub("cov","CoV",temp_output$Term)
  temp_output$Term <- gsub("amp","AMP",temp_output$Term)
  temp_output$Term <- gsub("atp","ATP",temp_output$Term)
  temp_output$Term <- gsub("Kegg","KEGG",temp_output$Term)
  temp_output$Term <- gsub("Reactome","REACTOME",temp_output$Term)
  temp_output$Term <- gsub("Go ","GO ",temp_output$Term)
  temp_output$Term <- gsub("trna","tRNA",temp_output$Term)
  temp_output$Term <- gsub("tca","TCA",temp_output$Term)
  temp_output$Term <- gsub("nmd","NMD",temp_output$Term)
  temp_output$Term <- gsub("rna","RNA",temp_output$Term)
  temp_output$Term <- gsub("dna","DNA",temp_output$Term)
  temp_output$Term <- gsub("Dna","DNA",temp_output$Term)
  temp_output$Term <- gsub("Rna","RNA",temp_output$Term)
  temp_output$Term <- gsub("Foxo","FOXO",temp_output$Term)
  temp_output$Term <- gsub("P53","p53",temp_output$Term)
  temp_output$Term <- gsub("tcr","TCR",temp_output$Term)
  temp_output$Term <- gsub("Vegf","VEGF",temp_output$Term)
  temp_output$Term <- gsub("Ppar","PPAR",temp_output$Term)
  temp_output$Term <- gsub("alzheimers","Alzheimers",temp_output$Term)
  temp_output$Term <- gsub("parkinsons","Parkinsons",temp_output$Term)
  temp_output$Term <- gsub("huntingtons","Huntingtons",temp_output$Term)
  temp_output$Term <- gsub("Citrate cycle TCA cycle","Citrate cycle (TCA cycle)",temp_output$Term)
  temp_output$Term <- gsub("mapk","MAPK",temp_output$Term)
  temp_output$Term <- gsub("REAC ","REACTOME ",temp_output$Term)
  temp_output$Term <- gsub("Nf-kappa b","NF-kB",temp_output$Term)
  temp_output$Term <- gsub("Jak-stat","JAK-STAT",temp_output$Term)
  temp_output$Term <- gsub("Hif-1","HIF-1",temp_output$Term)
  temp_output$Term <- gsub("Mapk","MAPK",temp_output$Term)
  temp_output$Term <- gsub("Bp ","",temp_output$Term)
  temp_output$Term <- gsub("Tnf","TNF",temp_output$Term)
  temp_output$Term <- gsub("The citric acid [(]TCA[)]","TCA",temp_output$Term)
  
  print(head(temp_output))
  library(reshape2)
  library(ggplot2)
  library(colorspace)
  melted_output <- melt(temp_output,id.vars="Term",measure.vars=c(1:(ncol(temp_output)-1)))
  colnames(melted_output) <- c("Function","Metric","Score")

  if (class(Gene_ratio_table)!="logical"){
    temp_ratio_output$Term <- rownames(temp_ratio_output)
    melted_ratio_output <- melt(temp_ratio_output,id.vars="Term",measure.vars=c(1:(ncol(temp_ratio_output)-1)))
    colnames(melted_ratio_output) <- c("Function","Metric","Ratio")
  }
  
  melted_output$Function <- as.character(melted_output$Function)
  if (flip==T){
    melted_output$Function <- factor(melted_output$Function, levels=unique(as.character(temp_output$Term)))
  }else{
    melted_output$Function <- factor(melted_output$Function, levels=rev(unique(as.character(temp_output$Term))))
  }
  melted_output$Metric <- factor(melted_output$Metric, 
                                 levels=unique(as.character(colnames(temp_output))))
  if (max_value!=F){
    if (sum(melted_output$Score>max_value)>0){
      melted_output[which(melted_output$Score>max_value),]$Score <- max_value
    }
    if (sum(melted_output$Score<(-max_value))>0){
      melted_output[which(melted_output$Score<(-1)*max_value),]$Score <- (-1)*max_value
    }
    max_abs_score <- abs(max_value)
  }else{
    max_abs_score <- max(abs(melted_output$Score),na.rm=T)
  }
  

  if (as.circles==T){
    if (class(Gene_ratio_table)!="logical"){
      melted_output$Ratio <- melted_ratio_output$Ratio
      max_ratio <- max(melted_output$Ratio)
      print(max_ratio)
      g <- ggplot(melted_output,aes(Metric,Function))+
        geom_point(aes(size=Ratio,fill=Score),shape=21,color="black")+
        scale_size_continuous(limits=c(0,0.2083),range = c(0,max_size))+
        scale_fill_gradient2(high="red3",low="blue3",mid=background_color,
                               guide = "colorbar",
                               limits=c(-max_abs_score,max_abs_score))+
        labs(x = "",y = "")+
        theme_bw(base_size = 18)+
        theme(axis.text.x = element_text(angle = 45, hjust=1,size=14,colour="black"),
              axis.text.y = element_text(size=15,colour="black"),
              axis.title=element_text(size=18,colour="black"))
    }else{
      melted_output$Ratio=0
      melted_output[,]
      g <- ggplot(melted_output,aes(Metric,Function))+
        geom_point(aes(size=abs(Score),fill=Score),shape=21,color="black")+
        scale_fill_gradient2(high="red3",low="blue3",mid=background_color,
                               guide = "colorbar",
                               limits=c(-max_abs_score,max_abs_score))+
        labs(x = "",y = "")+
        theme_bw(base_size = 18)+
        theme(axis.text.x = element_text(angle = 45,hjust=1,size=14,colour="black"),
              axis.text.y = element_text(size=15,colour="black"),
              axis.title=element_text(size=18,colour="black"))
    }
    
  }else{
    
  if (color_two_sided==T){
  g <- ggplot(melted_output, aes(Metric,Function,fill=Score)) + 
    geom_tile(colour = "black",size=0.1)+
    scale_fill_gradient2(high="red3",low="blue3",mid=background_color,
                         guide = "colorbar",
                         limits=c(-max_abs_score,max_abs_score),
                         oob=squish)+
    labs(x = "",y = "")+
    scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0))+
    theme_bw(base_size = 18)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1,size=14,colour="black"),
          axis.text.y = element_text(size=15,colour="black"),
          axis.title=element_text(size=18,colour="black"),
          panel.grid.major = element_line(colour="black"), panel.grid.minor = element_line(colour="black"))
  
  }else{
    if (min_value==F){
      g <- ggplot(melted_output, aes(Metric,Function,fill=Score)) + 
        geom_tile(colour = "black",size=0.1)+
        scale_fill_gradient(high="yellow2",low=background_color,guide = "colorbar",
                             limits=c(0,max_abs_score),
                             oob=squish)+
        labs(x = "",y = "")+
        scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0))+
        theme_bw(base_size = 18)+
        theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1,size=14,colour="black"),
              axis.text.y = element_text(size=15,colour="black"),
              axis.title=element_text(size=18,colour="black"),
              panel.grid.major = element_line(colour="black"), panel.grid.minor = element_line(colour="black"))
    }else{
      g <- ggplot(melted_output, aes(Metric,Function,fill=Score)) + 
        geom_tile(colour = "black",size=0.1)+
        scale_fill_gradient(high="yellow2",low=background_color,guide = "colorbar",
                             limits=c(min_value,max_abs_score),
                             oob=squish)+
        labs(x = "",y = "")+
        scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0))+
        theme_bw(base_size = 18)+
        theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1,size=14,colour="black"),
              axis.text.y = element_text(size=15,colour="black"),
              axis.title=element_text(size=18,colour="black"),
              panel.grid.major = element_line(colour="black"), panel.grid.minor = element_line(colour="black"))
    }
  }
  }
  
  if (color_labels_datasets!=F){
    g <- g+theme(axis.text.x=element_text(color=color_labels_datasets))
  }
  
  return(list(g,temp_output,levels(melted_output$Function),
              levels(melted_output$Metric)))
}