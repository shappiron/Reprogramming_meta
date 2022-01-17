#################################################
#This function is designed to run edgeR on specified design matrix

Run_edgeR <- function(expression_matrix,design_matrix,contrast_vector){
        library(edgeR)

        #Create edgeR object and calculate variance
        edger_matrix = DGEList(expression_matrix)
        edger_matrix = calcNormFactors(edger_matrix, method = "RLE")
        edger_matrix = estimateDisp(edger_matrix, design_matrix,robust=T)
        
        #Fit design matrix and calculate result
        glm.st_matrix = glmFit(edger_matrix, design_matrix)
        result <- topTags(glmLRT(glm.st_matrix,contrast=contrast_vector),n=Inf,adjust.method = "BH")$table
        result$P.Value <- result$PValue
        
        #Save logFC, pvalue and FDR
        #output_matrix <- result[,c("logFC","PValue","FDR")]
        output_matrix <- result
        return(output_matrix)        
        
}