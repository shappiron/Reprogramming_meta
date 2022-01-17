#########################################################
#This script creates dictionary connecting microarray IDs with Entrez IDs

Microarray_dictionary_create <- function(dataset){

        print("Initial number of probes:")
        print(nrow(dataset))
        index <- grepl("entrez",tolower(colnames(fData(dataset))))
        dataset_Entrez <- dataset[fData(dataset)[,index]!="" & !is.na(fData(dataset)[,index]),]

        print("Number of probes with Entrez ID:")
        print(nrow(dataset_Entrez))
        #1. Create library for every probe
        dataset_probetoEntrez <- new.env()
        for (i in 1:length(rownames(fData(dataset_Entrez))))
        {
                assign(rownames(fData(dataset_Entrez))[i],as.character(strsplit(as.character(fData(dataset_Entrez)[,index][i])," /// ")[[1]]),
                       envir=dataset_probetoEntrez)
        }
        print("Library 1 is created. ")

        #2. Create library for every Entrez ID
        dataset_Entreztoprobe <- new.env()
        for (i in 1:length(rownames(fData(dataset_Entrez))))
        {
                temp <- as.character(strsplit(as.character(fData(dataset_Entrez)[,index][i])," /// ")[[1]])
                for (j in 1:length(temp))
                {
                        if (exists(temp[j],dataset_Entreztoprobe))
                        {
                                assign(temp[j],c(get(temp[j],dataset_Entreztoprobe),rownames(fData(dataset_Entrez))[i]),envir=dataset_Entreztoprobe)
                        }
                        else
                        {
                                assign(temp[j],rownames(fData(dataset_Entrez))[i],envir=dataset_Entreztoprobe)
                        }
                }
        }
        print("Library 2 is created. ")
        print("Total number of Entrez IDs:")
        print(length(ls(dataset_Entreztoprobe)))

        return(list(dataset=dataset_Entrez,dataset_Entreztoprobe=dataset_Entreztoprobe,
                    dataset_probetoEntrez=dataset_probetoEntrez))
}



