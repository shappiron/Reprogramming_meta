####################################################################
#This function runs GSEA test between changes induced by compounds and
#corresponding longevity signatures genes
#Input is zscore matrix from single dataset (single permutation test)

Association_test <- function(signatures_genes_list, 
                             compounds_zscore_matrix,
                             number_permutations=1000,
                             function_path="D:/Documents/Science/PhD/Functions/"){
    
  source(paste(function_path,"FUN.GSEA.R",sep=""))
  source(paste(function_path,"FUN.GSEA_with_info.R",sep=""))

  rownames(compounds_zscore_matrix) <- toupper(rownames(compounds_zscore_matrix))
  
  #1) Prepare output matrix
  out_GSEA_matrix <- as.data.frame(matrix(nrow=length(signatures_genes_list), ncol=ncol(compounds_zscore_matrix)+1))
  rownames(out_GSEA_matrix) <- names(signatures_genes_list)
  colnames(out_GSEA_matrix) <- c("Signature", colnames(compounds_zscore_matrix))
  out_GSEA_matrix[,1] <- names(signatures_genes_list)
  out_pvalue_matrix <- out_GSEA_matrix
  
  direction_list <- list()
  gene_tables_list <- list()
  top_list <- list()
  es_list <- list()

  for (i in names(signatures_genes_list)){
    aaa <- Sys.time()
    
    print(i)
    #2) Select significant signature genes
    temp_signature_up <- intersect(toupper(signatures_genes_list[[i]]$Up), 
                                   toupper(rownames(compounds_zscore_matrix)))
    temp_signature_down <- intersect(toupper(signatures_genes_list[[i]]$Down),
                                     toupper(rownames(compounds_zscore_matrix)))
    
    print(paste("Number of significant genes in this signature: ",
                length(signatures_genes_list[[i]]$Up)+length(signatures_genes_list[[i]]$Down),sep=""))
    print(paste("Number of identified genes: ", length(temp_signature_up)+length(temp_signature_down),sep=""))
    
    #3) Run resampling cycle for pvalue calculation
    if (length(temp_signature_up) + length(temp_signature_down)==0){
      next
    }
    temp_profile <- compounds_zscore_matrix[,1]
    names(temp_profile) <- rownames(compounds_zscore_matrix)
    temp_profile <- temp_profile[complete.cases(temp_profile)]
    resampling_scores <- vector(length=number_permutations)
    resampling_scores_up <- resampling_scores
    resampling_scores_down <- resampling_scores
    for (t in 1:number_permutations){
        temp_genes_up <- sample(names(temp_profile), length(temp_signature_up), replace=F)
        temp_genes_down <- sample(setdiff(names(temp_profile), temp_genes_up), length(temp_signature_down), replace=F)
        resampling_scores_up[t] <- GSEA(temp_profile, temp_genes_up)
        resampling_scores_down[t] <- GSEA(temp_profile, temp_genes_down)
    }
    
    #4) Normalize resampling scores by standard deviations
#    print(plot(density((resampling_scores_up-resampling_scores_down)/2)))
#    print(abline(v = (similarity_up-similarity_down)/2,col="blue"))
#    print(plot(density(resampling_scores_up)))
#    print(lines(density(resampling_scores_down),col="red"))
    sd_up <- sd(resampling_scores_up)
    sd_down <- sd(resampling_scores_down)
    resampling_scores_up = resampling_scores_up/sd_up
    resampling_scores_down = resampling_scores_down/sd_down
    resampling_scores = (resampling_scores_up-resampling_scores_down)/2
#    print(plot(density(resampling_scores_up)))
#    print(lines(density(resampling_scores_down),col="red"))
    

    #5) Run GSEA for every compound and dose
    for (j in colnames(compounds_zscore_matrix)){
      #6) Prepare ranked profile
      temp_profile <- compounds_zscore_matrix[,j]
      names(temp_profile) <- rownames(compounds_zscore_matrix)
      temp_profile <- temp_profile[complete.cases(temp_profile)]
      temp_profile <- sort(temp_profile,decreasing = T)
      
      #7) Run GSEA
      temp_signature_up <- intersect(toupper(signatures_genes_list[[i]]$Up),
                                     toupper(names(temp_profile)))
      temp_signature_down <- intersect(toupper(signatures_genes_list[[i]]$Down),
                                       toupper(names(temp_profile)))
      
      similarity_up_list <- GSEA_with_info(temp_profile,temp_signature_up)
      similarity_down_list <- GSEA_with_info(temp_profile,temp_signature_down)
      #      print(plot(density(similarity_up_list$Genes$Index)))
      #      print(plot(density(similarity_down_list$Genes$Index)))
      gene_tables_list[[i]][[j]] <- list(Up=similarity_up_list$Genes,Down=similarity_down_list$Genes)
      direction_list[[i]][[j]] <- list(Up=similarity_up_list$Direction,Down=similarity_down_list$Direction)
      top_list[[i]][[j]] <- list(Up=similarity_up_list$Top,Down=similarity_down_list$Top)
      es_list[[i]][[j]] <- list(Up=similarity_up_list$ES,Down=similarity_down_list$ES)
      
      similarity_up <- similarity_up_list[[1]]
      similarity_down <- similarity_down_list[[1]]
      
      #8) Normalize scores by standard deviation of resampling scores
      es_list[[i]][[j]]$Up <- es_list[[i]][[j]]$Up/sd_up
      es_list[[i]][[j]]$Down <- es_list[[i]][[j]]$Down/sd_down
      NES_up <- similarity_up/sd_up
      NES_down <- similarity_down/sd_down
      #    print(plot(density((resampling_scores_up-resampling_scores_down)/2)))
      #    print(abline(v=(NES_up-NES_down)/2,col="blue"))
      
      #9) Calculate final similarity scores
      similarity_score <- (NES_up - NES_down) / 2
      out_GSEA_matrix[i,j] <- similarity_score
        
    }
    out_pvalue_matrix[i,-1] <- sapply(as.numeric(as.character(out_GSEA_matrix[i,-1])),function(x){sum(abs(resampling_scores)>abs(x))/number_permutations})

    bbb <- Sys.time()
    print(bbb-aaa)
    print("Resampling pvalue calculation completed")
    print(" ")
  }
  out_pvalue_matrix <- out_pvalue_matrix[complete.cases(out_pvalue_matrix),]
  out_GSEA_matrix <- out_GSEA_matrix[rownames(out_pvalue_matrix),]
  
  #10) Update p-values equal to 0
  for (i in 1:nrow(out_pvalue_matrix)){
    for (j in 1:ncol(out_pvalue_matrix)){
      if (out_pvalue_matrix[i,j]==0){
        out_pvalue_matrix[i,j]=1/number_permutations
      }
    }
  }
  
  return(list(GSEA=out_GSEA_matrix,pvalue=out_pvalue_matrix,
              Gene_table=gene_tables_list,Top=top_list,Direction=direction_list,NES_vector=es_list))
}
