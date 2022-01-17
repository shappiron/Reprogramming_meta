######################################
#This script uses previously found Ensembl to Entrez dictionary and transforms
#Ensembl data into Entrez

Ensembl_to_entrez <- function(dataset,dictionary){
        Entrez_dataset <- NULL
        i <- 0
        for (v in ls(dictionary))
        {
                i <- i+1
                if (length(dictionary[[v]])==1)
                {
                        Entrez_dataset <- rbind(Entrez_dataset,exprs(dataset)[dictionary[[v]],])
                }
                else
                {
                        print(v)
                        temp_Entrez <- apply(exprs(dataset)[dictionary[[v]],],2,FUN=sum)
                        Entrez_dataset <- rbind(Entrez_dataset,temp_Entrez)
                }
                rownames(Entrez_dataset)[i] <- v
        }
        
        return(Entrez_dataset)

}
