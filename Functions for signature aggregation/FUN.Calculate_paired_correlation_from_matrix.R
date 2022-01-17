# This function builds paired correlation matrix of logFC values
# for the specified number of top significant genes from matrices

Calculate_paired_correlation_from_matrix <- function(logFC_matrix,
                                                     pvalue_matrix,
                                                     number_top=F,
                                                     method="spearman"){
    
  #pvalue_matrix <- pvalue_matrix[,colnames(logFC_matrix)]
  #Contruct matrix template
  correlation_matrix <- matrix(nrow=ncol(logFC_matrix),ncol=ncol(logFC_matrix))
  rownames(correlation_matrix) = colnames(correlation_matrix) = colnames(logFC_matrix)
  diag(correlation_matrix) <- 1
  corr_matrix_list <- list(corr=correlation_matrix,
                           P.Value=correlation_matrix,
                           FDR=correlation_matrix,
                           Zscore=correlation_matrix)
  #diag(corr_matrix_list$Zscore) <- 1
  for (i in 1:(nrow(correlation_matrix)-1)){
    name1 <- rownames(correlation_matrix)[i]
    #print(name1)
    for (j in (i+1):(ncol(correlation_matrix))){
      name2 <- rownames(correlation_matrix)[j]
      temp_data_logFC <- logFC_matrix[,c(name1, name2)] #take two datasets
      temp_data_logFC <- temp_data_logFC[complete.cases(temp_data_logFC),] #drop NaN rows
      temp_data_pvalue <- pvalue_matrix[rownames(temp_data_logFC),c(name1,name2)]
      
      #print(nrow(temp_data_logFC))
      
      # extract top genes from each of two dataset -> union them -> slice this union from the original one
      if (number_top!=F & number_top!=0){
        temp_data_pvalue <- temp_data_pvalue[order(temp_data_pvalue[,name1],decreasing = F),]
        data1_top <- rownames(temp_data_pvalue)[1:number_top]
        temp_data_pvalue <- temp_data_pvalue[order(temp_data_pvalue[,name2],decreasing = F),]
        data2_top <- rownames(temp_data_pvalue)[1:number_top]
        top_genes <- union(data1_top,data2_top)
        temp_data_logFC <- temp_data_logFC[top_genes,]
      }     
        
      #??the first row here is redundant
      temp_corr <- cor(temp_data_logFC[,1], temp_data_logFC[,2], method=method)
      temp_corr_test <- cor.test(temp_data_logFC[,1], temp_data_logFC[,2], method=method)
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
  ##? looks like redundant block  
  diag(corr_matrix_list$Zscore) <- max(corr_matrix_list$Zscore)
  corr_matrix_list$Zscore_sign <- corr_matrix_list$Zscore
  for (i in (1:(nrow(corr_matrix_list$Zscore) - 1))){
    for (j in ((i+1):nrow(corr_matrix_list$Zscore))){
      #if p-value <0.1 (corresponds -log10(p) = 1) then zero this Z-score
      if (abs(corr_matrix_list$Zscore_sign[i,j]) < 1){
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
