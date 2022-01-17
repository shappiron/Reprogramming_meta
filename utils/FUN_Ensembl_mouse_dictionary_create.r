######################################################
#This function creates mice dictionary for Ensembl-Entrez mapping

Ensembl_mouse_dictionary_create <- function(dataset){
        #1) First let's create table where every Ensembl has corresponding Entrez ID; we will remove samples that do not map or map to several Entrez IDs
        library(org.Mm.eg.db)
        xx <- as.list(org.Mm.egENSEMBL2EG)
        ensembltoentrez <- data.frame(ID=rownames(dataset),ENTREZID=NA)
        rownames(ensembltoentrez) <- ensembltoentrez$ID
        l=0
        r=0
        for (i in 1:nrow(dataset))
        {
                if (length(xx[[rownames(dataset)[i]]])==1)
                {
                        ensembltoentrez[i,2] <- xx[[rownames(dataset)[i]]]
                }
                else if (is.null(xx[[rownames(dataset)[i]]]))
                {
                        r=r+1
                }
                else
                {
                        l=l+1 
                }
        }

        #Let's remove Ensembl IDs that do not map to any Entrez ID and that map to several Entrez IDs
        ensembltoentrez_noNA <- ensembltoentrez[!is.na(ensembltoentrez$ENTREZID),]
        dataset <- dataset[rownames(ensembltoentrez_noNA),]

        #2) Let's create library for every Entrez ID to see if any of them are mapped by more than 1 Ensembl ID
        EntreztoID_rnaseq <- new.env()
        for (i in 1:length(rownames(ensembltoentrez_noNA)))
        {
                temp <- strsplit(ensembltoentrez_noNA$ENTREZID[i]," /// ")[[1]]
                for (j in 1:length(temp))
                {
                        if (exists(temp[j],EntreztoID_rnaseq))
                        {
                                assign(temp[j],c(get(temp[j],EntreztoID_rnaseq),rownames(ensembltoentrez_noNA)[i]),envir=EntreztoID_rnaseq)
                        }
                        else
                        {
                                assign(temp[j],rownames(ensembltoentrez_noNA)[i],envir=EntreztoID_rnaseq)
                        }
                }
        }

        return(EntreztoID_rnaseq)
}