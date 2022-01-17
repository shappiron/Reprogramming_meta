###############################################
#This script performs RLE normalization and returns normalized dataset

RLE_normalization  <- function(original_dataset){
         library(edgeR)
        #First, we need to get library sizes vector
        lib.size = apply(exprs(original_dataset), 2, sum)
        
        #Now we should obtain normalization factors for each samples using RLE technique
        edger.rle = lib.size * calcNormFactors(exprs(original_dataset), method="RLE")
        
        #Finally we should divide expression values on normalization factors
        RLE_dataset <- original_dataset
        exprs(RLE_dataset) = sweep(exprs(original_dataset), 2, edger.rle, "/")*10^7
        
        return(RLE_dataset)
}


