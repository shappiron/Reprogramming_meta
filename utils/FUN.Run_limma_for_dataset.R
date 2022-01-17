#function adjusted_pvalue_calculator creates matrix where every gene from certain list of genes
#corresponds to the row and contains p-value and adjusted pvalue for its Fold Change for every phenotype pair
#or NAs (if this gene is not measured in the dataset or is filtered by sds).
#Names of columns are taken from samplename_matrix (input).
#Samplename_matrix must contain names of samples used in matrix
#in first column and corresponding samples from expression set object.
#Controls are taken from samples not considered in samplename_matrix.

Run_limma_for_dataset <- function(exprs_data,phenodata,
                                  control_group,batch_col=""){

        library(limma)
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

          fit <- lmFit(data_for_analysis,design)
          fit <- eBayes(fit)
          
          result <- topTable(fit,coef="Group",adjust.method = "BH", sort.by="none", 
                             n=Inf, lfc=log2(1),confint = T)
          result$PValue <- result$P.Value
          result$FDR <- result$adj.P.Val
          result$SE <- (result[,3]-result[,1])/1.96
          output_list[[temp_phenotype]] <- result
        }
       
        return(output_list)
}