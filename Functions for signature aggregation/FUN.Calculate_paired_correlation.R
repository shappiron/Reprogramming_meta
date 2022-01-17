# This function builds paired correlation matrix of logFC values
# for the specified number of top significant genes

Calculate_paired_correlation <- function(logFC_list, number_top=F, method="spearman"){
  
  correlation_matrix <- matrix(nrow=length(logFC_list), ncol=length(logFC_list))
  rownames(correlation_matrix) = colnames(correlation_matrix) = names(logFC_list)
  diag(correlation_matrix) <- 1
    
  corr_matrix_list <- list(corr=correlation_matrix,
                           P.Value=correlation_matrix,
                           FDR=correlation_matrix,
                           Zscore=correlation_matrix)
  #diag(corr_matrix_list$Zscore) <- 1
  for (i in 1:(nrow(correlation_matrix)-1)){
    name1 <- rownames(correlation_matrix)[i]
    #print(name1)
    signature_table1 <- logFC_list[[name1]]
    for (j in (i+1):(ncol(correlation_matrix))){
      name2 <- rownames(correlation_matrix)[j]
      signature_table2 <- logFC_list[[name2]]
      common_genes <- intersect(rownames(signature_table1),rownames(signature_table2))
      print(length(common_genes))
      data1 <- signature_table1[common_genes,]
      data2 <- signature_table2[common_genes,]
      
      if (number_top!=F & number_top!=0){
        data1 <- data1[order(data1$P.Value,decreasing = F),]
        data1_top <- rownames(data1)[1:number_top]
        data2 <- data2[order(data2$P.Value,decreasing = F),]
        data2_top <- rownames(data2)[1:number_top]
        top_genes <- union(data1_top,data2_top)
        data1 <- data1[top_genes,]
        data2 <- data2[top_genes,]
      }     
       
      #redundant line as before
      temp_corr <- cor(data1$logFC, data2$logFC, method=method) #????
      temp_corr_test <- cor.test(data1$logFC,data2$logFC,method=method)
      temp_corr <- temp_corr_test$estimate
      temp_corr_pval <- temp_corr_test$p.value
      if (temp_corr_pval==0){
        temp_corr_pval=10^(-100)
      }
      
      corr_matrix_list$corr[i,j] = corr_matrix_list$corr[j,i] = temp_corr
      corr_matrix_list$P.Value[i,j] = corr_matrix_list$P.Value[j,i] = temp_corr_pval
      corr_matrix_list$Zscore[i,j] = corr_matrix_list$Zscore[j,i] = -log10(temp_corr_pval)*sign(temp_corr)
      
    }
    #print(" ")
  }
  diag(corr_matrix_list$Zscore) <- max(corr_matrix_list$Zscore)
  corr_matrix_list$Zscore_sign <- corr_matrix_list$Zscore
  for (i in (1:(nrow(corr_matrix_list$Zscore) - 1))){
    for (j in ((i+1):nrow(corr_matrix_list$Zscore))){
      if (abs(corr_matrix_list$Zscore_sign[i,j])<1){
        corr_matrix_list$Zscore_sign[i,j] = corr_matrix_list$Zscore_sign[j,i] = 0
      }
    }
  }
  
  
  FDR_table <- corr_matrix_list$P.Value
  #extract upper triangular of the matrix and apply FDR correction
  tmp <- matrix(0, nrow=nrow(FDR_table), ncol=ncol(FDR_table))
  padj <- matrix(p.adjust(FDR_table[upper.tri(FDR_table)], method="BH"))
  tmp[upper.tri(tmp)] <- padj
  tmp <- tmp + t(tmp)
  diag(tmp) <- 1
  FDR_table <- as.data.frame(tmp, nrow=nrow(corr_matrix_list$P.Value))
  ###
  rownames(FDR_table) <- rownames(corr_matrix_list$P.Value)
  colnames(FDR_table) <- colnames(corr_matrix_list$P.Value)
  
  corr_matrix_list[["FDR"]] <- FDR_table
  
  return(corr_matrix_list)
}
