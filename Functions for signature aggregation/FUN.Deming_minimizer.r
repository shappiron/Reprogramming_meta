#This functions performs multiple Deming regression using logFC matrix and parameters for filltering:
#number of genes, threshold for FDR and threshold for correlation coefficient


Deming_minimizer = function(logFC_list, 
                            number_top=0,
                            coeff_threshold=0.1, #correlation coefficient
                            FDR_threshold=0.05,
                            function_path="D:/Documents/Science/PhD/Functions/"){
  a = Sys.time()
  source(paste(function_path,"FUN.Calculate_paired_correlation.R",sep=""))
  invisible(capture.output(temp_corr_list <- Calculate_paired_correlation(logFC_list, 
                                                                          number_top, 
                                                                          method="pearson")))
  
  temp_matrix <- temp_corr_list$corr
  temp_matrix <- (temp_corr_list$FDR<FDR_threshold & abs(temp_matrix)>coeff_threshold) * sign(temp_matrix)

    #this block create a double key dictionary with union of top genes between datasets
  paired_genes <- list()
  for (i in 1:(length(logFC_list)-1)){
    namei = names(logFC_list)[i]
    for (j in (i+1):length(logFC_list)){
      namej = names(logFC_list)[j]
      common_genes <- intersect(rownames(logFC_list[[namei]]), rownames(logFC_list[[namej]]))
      data1 <- logFC_list[[namei]][common_genes,]
      data2 <- logFC_list[[namej]][common_genes,]
      
      if (number_top==F | number_top==0){
        top_genes <- common_genes
      }else{
        data1 <- data1[order(data1$P.Value, decreasing = F),]
        data1_top <- rownames(data1)[1:number_top]
        data2 <- data2[order(data2$P.Value, decreasing = F),]
        data2_top <- rownames(data2)[1:number_top]
        top_genes <- union(data1_top,data2_top)
      }
      paired_genes[[namei]][[namej]] <- top_genes
    }
  }

    #This function searches for optimal normalization coefficients over all datasets having common genes
  fn = function(k_no_first){
    count = 0
    k = c()
    k[1] = 1
    k[2:length(logFC_list)] = k_no_first
    res = 0
    for (i in 1:(length(logFC_list)-1)){
      namei = names(logFC_list)[i]
      for (j in (i + 1):length(logFC_list)){
        namej = names(logFC_list)[j]
        if (temp_matrix[namei, namej] == 0){ #if both FDR and Corr conditions were not satisfied - skip
          next
        }
        totalrownames = paired_genes[[namei]][[namej]]
        ai = logFC_list[[namei]][totalrownames, ]$logFC
        aj = logFC_list[[namej]][totalrownames, ]$logFC
        if (temp_matrix[namei, namej] == -1){
          aj = (-1)*aj
        }
        res = res + sum(
          (((aj - (k[j]/k[i])*ai)^2)*((ai - (k[i]/k[j])*aj)^2))/
            (((aj - (k[j]/k[i])*ai)^2)+((ai - (k[i]/k[j])*aj)^2))) / length(totalrownames)
        count = count+1
      }
    }
  res = res/count
  return(res)
  }
  
  kvec = rnorm(length(logFC_list) - 1, 1, 1) #generate initial coefs around 1 as mean
  ptm <- proc.time()
  #optimized = optim(kvec, fn, lower = 0.01, upper = 100, method = "L-BFGS-B", control = list(factr = 1e3))
  optimized = optim(kvec, fn, lower = 0.01, upper = 100, method = "L-BFGS-B")
  proc.time() - ptm
  
  kres = c(1, optimized$par)
  minimum = optimized$value
  bigres = list(kres, minimum)
  names(bigres) = c("coefs", "minimum")
  b = Sys.time()
  print(b-a)
  return(bigres)
  
}