######################################
#This script uses previously found microarray to Entrez dictionary and transforms
#microarray data into Entrez

Microarray_to_Entrez <- function(dataset_Entrez,dataset_Entreztoprobe,dataset_probetoEntrez){


i <- 0
for (v in ls(dataset_Entreztoprobe))
{
        i <- i+1
        
        if (i==1)
        {
                if (length(dataset_Entreztoprobe[[v]])==1)
                {
                        Entrez_dataset_exprs <- t(as.data.frame(exprs(dataset_Entrez)[dataset_Entreztoprobe[[v]],]))
                }
                else
                {
                        t <- rep(0,length(dataset_Entreztoprobe[[v]]))
                        for (r in 1:length(dataset_Entreztoprobe[[v]]))
                        {
                                if (length(dataset_probetoEntrez[[dataset_Entreztoprobe[[v]][r]]])==1)
                                {
                                        t[r] <- 1
                                }
                        }
                        if (identical(t,rep(0,length(dataset_Entreztoprobe[[v]]))))
                        {
                                temp_Entrez <- apply(dataset_Entrez[dataset_Entreztoprobe[[v]],],2,FUN=mean)
                        }
                        else
                        {
                                temp_IDs <- dataset_Entreztoprobe[[v]][t==1]
                                temp_Entrez <- apply(dataset_Entrez[temp_IDs,],2,FUN=mean)
                        }
                        Entrez_dataset_exprs <- t(as.data.frame(temp_Entrez))
                }
                
        }
        else
        {
                if (length(dataset_Entreztoprobe[[v]])==1)
                {
                        Entrez_dataset_exprs <- rbind(Entrez_dataset_exprs,exprs(dataset_Entrez)[dataset_Entreztoprobe[[v]],])
                }
                else
                {
                        t <- rep(0,length(dataset_Entreztoprobe[[v]]))
                        for (r in 1:length(dataset_Entreztoprobe[[v]]))
                        {
                                if (length(dataset_probetoEntrez[[dataset_Entreztoprobe[[v]][r]]])==1) 
                                {
                                        t[r] <- 1
                                }
                        }
                        if (identical(t,rep(0,length(dataset_Entreztoprobe[[v]]))))
                        {
                                temp_Entrez <- apply(dataset_Entrez[dataset_Entreztoprobe[[v]],],2,FUN=mean)
                        }
                        else
                        {
                                temp_IDs <- dataset_Entreztoprobe[[v]][t==1]
                                temp_Entrez <- apply(dataset_Entrez[temp_IDs,],2,FUN=mean)
                        }
                        Entrez_dataset_exprs <- rbind(Entrez_dataset_exprs,temp_Entrez)
                }
        }
        rownames(Entrez_dataset_exprs)[i] <- v
}
print(i)
print(length(rownames(Entrez_dataset_exprs)))

#Now we can create expression set object from obtained dataset
pheno_dataset <- phenoData(dataset_Entrez)
Entrez_dataset <- new("ExpressionSet", exprs=Entrez_dataset_exprs, phenoData=pheno_dataset)
dim(exprs(Entrez_dataset))

return(Entrez_dataset)
}