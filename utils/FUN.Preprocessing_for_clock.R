preprocessing_for_clock <- function(exprs_data,na_fill_data,
                                    map_to_orthologs=F,map_from="Homo sapiens",
                                    map_to="Mus musculus",map_type="Entrez",
                                    function_path="/utils/"){
  
  #1. Map to orthologs (if needed)
  # source(paste0(function_path,"FUN.Align_to_normal_dist.R"))
  # source(paste0(function_path,"FUN.Map_to_orthologs.R"))

  
  exprs_data <- t(exprs_data)
  # if (map_to_orthologs==T){
  #   print("Mapping to orthologs...")
  #   source(paste0(function_path,"FUN.Map_to_orthologs.R"))
  #   mapping_table <- Map_to_orthologs(colnames(exprs_data),
  #                               species_from = map_from,species_to = map_to,
  #                               type=map_type)
  #   exprs_data <- exprs_data[,!is.na(colnames_map[colnames(exprs_data),]$To_Entrez)]
  #   colnames(exprs_data) <- as.character(colnames_map[colnames(exprs_data),]$To_Entrez)
  # }
  
  #2. Scale and transform to norm distribution
  #Scaled data
  print("Scaling...")
  exprs_scale <- t(scale(t(exprs_data)))
  # plot(density(as.numeric(exprs_scale[1,]),na.rm = T))
  # for (i in 1:nrow(exprs_scale)){
  #   lines(density(as.numeric(exprs_scale[i,]),na.rm=T),col=i)
  # }
  
  #Normally distributed data
  print("Transforming to normal distribution...")
  temp <- Align_to_normal_dist(exprs_data,by="rows")
  exprs_norm <- temp
  # plot(density(as.numeric(exprs_norm[1,]),na.rm = T))
  # for (i in 1:nrow(exprs_norm)){
  #   lines(density(as.numeric(exprs_norm[i,]),na.rm=T),col=i)
  # }
  # dim(exprs_norm)

  #3. Fill NAs with average values
  print("Filling NAs with mean values...")
  common_genes <- intersect(rownames(na_fill_data),colnames(exprs_norm))
  extra_genes <- setdiff(rownames(na_fill_data),colnames(exprs_norm))
  print("Total number of genes in the clock:")
  print(nrow(na_fill_data))
  print("Number of genes presented in the data:")
  print(length(common_genes))

  exprs_scale_final = exprs_norm_final = matrix(nrow=nrow(exprs_scale),
                                                ncol=nrow(na_fill_data))
  rownames(exprs_norm_final) = rownames(exprs_scale_final) = rownames(exprs_scale)
  colnames(exprs_norm_final) = colnames(exprs_scale_final) = rownames(na_fill_data)

  exprs_scale_final[,extra_genes] <- rep(na_fill_data[extra_genes,"Scaled_value"],
                                         each=nrow(exprs_scale_final))
  exprs_norm_final[,extra_genes] <- rep(na_fill_data[extra_genes,"Norm_value"],
                                         each=nrow(exprs_norm_final))
  exprs_scale_final[,common_genes] <- exprs_scale[,common_genes]
  exprs_norm_final[,common_genes] <- exprs_norm[,common_genes]
  missing_genes <- colnames(exprs_scale_final)[which(!complete.cases(t(exprs_scale_final)))]
  #return(missing_genes)
  if (length(missing_genes)>0){
    for (temp_gene in missing_genes){
      exprs_norm_final[which(is.na(exprs_norm_final[,temp_gene])),temp_gene] <- na_fill_data[temp_gene,"Norm_value"]
      exprs_scale_final[which(is.na(exprs_scale_final[,temp_gene])),temp_gene] <- na_fill_data[temp_gene,"Scaled_value"]
    }
  }

  return(list(Scaled_data = exprs_scale_final,
              Normal_data = exprs_norm_final))
}