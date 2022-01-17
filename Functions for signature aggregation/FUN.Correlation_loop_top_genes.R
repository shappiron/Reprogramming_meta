#This script runs loop for the best number of top genes
Correlation_loop_top_genes <- function(logFC_table,
                                       pvalue_table,
                                       vector_number_genes,
                                       function_path,
                                       method="pearson",
                                       FDR_thresh=0.05,
                                       corr_thresh=0.1,
                                       samplenames=F){
  
  #Select samples for correlation matrix ---- must be deprecated in future
  if (samplenames==F){
    pvalue_table <- pvalue_table[,colnames(logFC_table)]
  }else{
    logFC_table <- logFC_table[,samplenames]
    pvalue_table <- pvalue_table[,samplenames]
  }
  
  #Run correlation matrix loop for different number of top genes
  source(paste0(function_path,"FUN.Calculate_paired_correlation_from_matrix.R"))
  #get pairwise correlation matrix and some additional information
  complete_cor_matrix <- Calculate_paired_correlation_from_matrix(logFC_table,
                                                              pvalue_table,
                                                              number_top = F, method = method)
  
  significant_pairs <- list(FDR=c(), Coeff=c(), Both=c())
  for (number_genes in vector_number_genes){
    #get pairwise correlation matrix and some additional information for top N significant genes
    temp_corr_list <- Calculate_paired_correlation_from_matrix(logFC_table,
                                                                    pvalue_table,
                                                                    number_top = number_genes,
                                                                    method = method)
    
    significant_pairs$FDR <- c(significant_pairs$FDR,
                               sum(temp_corr_list$FDR < FDR_thresh) - nrow(temp_corr_list$FDR)) #some problems can be here
    significant_pairs$Coeff <- c(significant_pairs$Coeff,
                               sum(abs(temp_corr_list$corr) > corr_thresh) - nrow(temp_corr_list$FDR))
    significant_pairs$Both <- c(significant_pairs$Both,
                                sum(temp_corr_list$FDR < FDR_thresh & abs(temp_corr_list$corr) > corr_thresh)-nrow(temp_corr_list$FDR))
  }
  
  top_values <- data.frame(Top=vector_number_genes, FDR=0, Coeff=0, Both=0)
  for (a in names(significant_pairs)){
    top_values[,a] <- significant_pairs[[a]]
  }
  print("Optimal number of top genes:")

  #print(top_values[which.max(rev(top_values$Both)),])
  #invert the table to get maximum value of Top column by maximum value of Both column
  #elsewise the first (least) maximum will be returned
  top_gene_number <- as.numeric(as.character(top_values[seq(nrow(top_values),1,by=-1),][which.max(rev(top_values$Both)),]$Top))
  print(top_gene_number)
    
  #PLOTTING
  library(reshape2)
  top_values2 <- melt(top_values,measure.vars = 2:4)
  library(ggplot2)
  print(ggplot(top_values2,aes(Top,value,color=variable))+
    geom_line(lwd=1.5)+
    geom_vline(xintercept = 0,lwd=1)+
    geom_vline(xintercept=top_gene_number,col="Blue",lwd=1,lty=2)+
    geom_hline(yintercept = 0,lwd=1)+
    geom_hline(yintercept = nrow(temp_corr_list$FDR)*(nrow(temp_corr_list$FDR)-1),
               lwd=1,lty=2)+
    xlab("# top genes")+
    ylab("# significant correlations")+
    scale_x_continuous(breaks=seq(0,2000,by=100))+
    scale_y_continuous(breaks=seq(0,400,by=50))+
    scale_color_discrete(name="Threshold",
                         labels=c(paste0("P.adjust < ", FDR_thresh),
                                  paste0("Coefficient > ", corr_thresh),
                                  "Both"))+
    theme_bw()+
    theme(strip.text.x = element_text(size=18,colour="black",face="bold"),
          strip.background = element_rect(fill="white"),
          axis.text.x = element_text(size=12,colour="black"),
          axis.text.y = element_text(size=12,colour="black"),
          axis.title=element_text(size=22),
          title=element_text(size = 18),
          legend.position = c(0.88,0.74),
          legend.justification = c(0.73,0.15),
          legend.title=element_text(size=27,face="bold"),
          legend.text = element_text(size=23),
          plot.title = element_text(hjust=0.5)))
  
  optimal_correlation <- Calculate_paired_correlation_from_matrix(logFC_table,
                                                                  pvalue_table,
                                                                  number_top = top_gene_number,
                                                                  method = method)
  
  return(list(Complete=complete_cor_matrix,Optimal=optimal_correlation,Number_top_genes=top_gene_number))
  
}