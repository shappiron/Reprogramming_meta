Run_edgeR_for_dataset <- function(exprs_data,phenodata,
                                  control_group,batch_col=""){
  

        #make some corrections (if needed)
        samplename_matrix  <- phenodata
        samplename_matrix = as.matrix(samplename_matrix)
  
        #create output matrices
        phenotype_pairs <- levels(factor(samplename_matrix[,"Group"]))
        phenotype_pairs <- setdiff(phenotype_pairs,control_group)
        output_list <- list()
  
        for (i in seq_along(phenotype_pairs)){
          temp_phenotype <- phenotype_pairs[i]
          temp_phenodata <- phenodata[phenodata$Group %in% c(temp_phenotype,control_group),]
          data_for_analysis = exprs(exprs_data)
          data_for_analysis = data_for_analysis[,rownames(temp_phenodata)]
    
          control = rep(1,ncol(data_for_analysis))
          design = model.matrix(~control-1)
          rownames(design) <- rownames(temp_phenodata)
          design <- as.data.frame(design)
    
          if (batch_col!=""){
            temp_phenodata$Batch <- temp_phenodata[,batch_col]
            design = model.matrix(~temp_phenodata$Batch)
            rownames(design) <- rownames(temp_phenodata)
            design <- as.data.frame(design)
          }
          design[,"Group"] <- as.numeric(temp_phenodata$Group==temp_phenotype)
    
          result <- Run_edgeR(data_for_analysis,design,c(rep(0,ncol(design)-1),1))
          output_list[[temp_phenotype]] <- result
        }
       
        return(output_list)
}