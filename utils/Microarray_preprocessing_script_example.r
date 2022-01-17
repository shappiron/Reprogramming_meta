#This is an example of scipt for preprocessing of individual gene expression datasets
#In this example liver samples of GSE55272 will be analyzed.
#The plan is to identify differentially expressed genes between wildtype mice and
#mice deficient for Myc gene (heterozygous deletion)

#0. Specify directory where functions are stored
function_path <- c("D:/Documents/Science/PhD/Functions/")

#1. Download data - GSE55272 (microarray)
library(lumi)
library(GEOquery)
current_gse <- "GSE55272"
data <- getGEO(current_gse, GSEMatrix =TRUE)
if (length(data) > 1) idx <- grep("GPL6103", attr(data, "names")) else idx <- 1
data <- data[[idx]]
pData(data)


#2. Update sample annotation
pheno_mod <- data.frame(ID=rownames(pData(data)),
                        Species=as.character(pData(data)$organism_ch1),
                        Group=as.character(pData(data)[,"genotype:ch1"]),
                        Age=as.character(pData(data)[,"age:ch1"]),
                        Sex=as.character(pData(data)[,"Sex:ch1"]),
                        Tissue=as.character(pData(data)[,"tissue:ch1"]),
                        Strain=as.character(pData(data)[,"strain:ch1"]))
rownames(pheno_mod) <- as.character(pheno_mod$ID)
#write.csv(pheno_mod,"Phenodata_mod.csv")
common_rows <- intersect(rownames(pheno_mod),colnames(data))
data <- data[,common_rows]
print(identical(colnames(data),rownames(pheno_mod)))
pData(data) <- pheno_mod
pData(data)
dim(data)

#Keep only liver samples
data <- data[,data$Tissue=="Liver"]
dim(data)


#3. Remove unexpressed genes and normalize data
#Density plot
plot(density(exprs(data[,1])))
for (i in 1:ncol(data)){
  lines(density(exprs(data[,i])),col=i)
}

#A) Remove unexpressed genes (for RNA-seq)
#library(edgeR)
#reads_thresh <- 10
#percentage_thresh <- 0.5
#r <- apply(exprs(data),1,function(x){sum(x>=reads_thresh)})>=ncol(data)*percentage_thresh
#table(r)
#density(data[r,])
#data_filtered <- data[r,]

#B) Log-transform data
#In our case data is already log-transformed
#exprs(data) <- log2(exprs(data)+1)

#C) Perform sample normalzation 
#scaling and quantile normalization - for microarray;
#RLE or TMM - for RNAseq
data_norm <- data
library(preprocessCore)
exprs(data_norm) <- scale(exprs(data_norm))
exprs(data_norm) <- log2(normalize.quantiles.robust(2^exprs(data_norm), remove.extreme="both", n.remove=2))
plot(density(exprs(data_norm)[,1]))
for (i in 2:ncol(data_norm)){
  lines(density(exprs(data_norm)[,i]),col=i)
}


#4. Clustering
#A) PCA
library(rgl)
#Specify names of columns corresponding to different interventions
colnames(pData(data_norm))
group_factors <- c("Age","Group")

pheno_groups <- apply(pData(data_norm)[,group_factors],1,function(x){factor(paste(as.character(x),collapse ="_"))})
temp <- exprs(data_norm)
pca_model <- prcomp(t(temp),scale=T)
variance_explained <- pca_model$sdev^2/sum(pca_model$sdev^2)*100

library(ggplot2)
ggplot(data=as.data.frame(pca_model$x[,1:2]),aes(PC1,PC2))+
  geom_point(aes(colour=pheno_groups),size=5)+
  theme_bw(base_size = 22)+
  geom_hline(yintercept = 0,lwd=1,lty=1)+
  geom_vline(xintercept = 0,lwd=1,lty=1)+
  scale_y_continuous(limits=c(-max(abs(pca_model$x[,2]))*1.1,max(abs(pca_model$x[,2]))*1.1))+
  scale_x_continuous(limits=c(-max(abs(pca_model$x[,1]))*1.1,max(abs(pca_model$x[,1]))*1.1))+
  labs(x=paste0("PC1 (",round(variance_explained[1],1),"%)"),
       y=paste0("PC2 (",round(variance_explained[2],1),"%)"))+
  theme(legend.position = c(0.8,0.55),
        legend.justification = c(0.73,0.15),
        axis.text.x = element_text(angle = 0, hjust = 0.5,size=14,colour="black"),
        axis.text.y = element_text(size=14,colour="black"),
        axis.title=element_text(size=20,face = "bold"),
        title=element_text(size = 20),
        legend.title=element_blank(),
        legend.text = element_text(size=22),
        legend.key.size = unit(1.5, 'lines'),
        plot.title = element_text(hjust=0.5))
#Nice separation

#B) MDS
dist.matrix_rnaseq <- as.dist((1-cor(exprs(data_norm), method="spearman"))/2)
mdsdiffage_rnaseq <- cmdscale(dist.matrix_rnaseq)
colnames(mdsdiffage_rnaseq) <- c("Coordinate1","Coordinate2")
mdsdiffage_rnaseq <- as.data.frame(mdsdiffage_rnaseq)

plot(mdsdiffage_rnaseq[,1],mdsdiffage_rnaseq[,2], pch=16, 
     col=rainbow(length(levels(pheno_groups)))[as.numeric(pheno_groups)],
     cex=1.5,xlab="Coordinate 1", ylab="Coordinate 2")
text(mdsdiffage_rnaseq[,1]+0.001,mdsdiffage_rnaseq[,2]+0.0005,
     rownames(pData(data_norm)),cex=0.9,
     col=rainbow(length(levels(pheno_groups)))[as.numeric(pheno_groups)])
legend('topright', levels(pheno_groups), pch=16, col=rainbow(length(levels(pheno_groups))))

ggplot(data=mdsdiffage_rnaseq,aes(Coordinate1,Coordinate2))+
  geom_point(aes(colour=pheno_groups),size=5)+
  theme_bw(base_size = 22)+
  geom_hline(yintercept = 0,lwd=1,lty=1)+
  geom_vline(xintercept = 0,lwd=1,lty=1)+
  scale_y_continuous(limits=c(-max(abs(mdsdiffage_rnaseq[,2]))*1.1,
                              max(abs(mdsdiffage_rnaseq[,2]))*1.1))+
  scale_x_continuous(limits=c(-max(abs(mdsdiffage_rnaseq[,1]))*1.1,
                              max(abs(mdsdiffage_rnaseq[,1]))*1.1))+
  labs(x="Coordinate 1",
       y="Coordinate 2")+
  theme(legend.position = c(0.55,0.55),
        legend.justification = c(0.73,0.15),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        axis.text.x = element_text(angle = 0, hjust = 0.5,size=14,colour="black"),
        axis.text.y = element_text(size=14,colour="black"),
        axis.title=element_text(size=20,face = "bold"),
        title=element_text(size = 20),
        legend.title=element_text(size=20,face="bold"),
        legend.text = element_text(size=22),
        legend.key.size = unit(1.5, 'lines'),
        plot.title = element_text(hjust=0.5))


#5. Translate gene IDs to Entrez IDs and run diff expression analysis
#A) First, we need to find annotation for microarray IDs translating them to Entrez IDs
head(fData(data_norm))
#There is no Entrez ID annotation in this data. Will add it manually using the platform annotation
map = fData(data_norm) # annotation, see '?fData'
map_new <- matrix(nrow=nrow(map),ncol=2)
map_new <- as.data.frame(map_new)
colnames(map_new) <- c("probe ID", "ENTREZ")
map_new[,1] <- rownames(map)
rownames(map_new) <- map_new[,1]
dim(map_new)

library(mogene10sttranscriptcluster.db)
for (i in 1:nrow(map_new))
{
  temp_entrez <- paste(mogene10sttranscriptclusterENTREZID[[rownames(map_new)[i]]],
                       collapse=" /// ")
  map_new[i,2] <- temp_entrez
}
identical(rownames(map_new),rownames(map))
fData(data_norm)$Entrez <- as.character(map_new$ENTREZ)
head(fData(data_norm)[fData(data_norm)$Entrez!="NA",])
sum(map_new$ENTREZ!="NA")
#For 22088 microarray IDs Entrez ID was identified


#B) Diff expression analysis

#We want to find MYC deficiency associated gene expression changes for both ages

#Specify factors corresponding to control group
print(levels(pheno_groups))
pheno_dataframe <- data.frame(ID=colnames(data_norm))
rownames(pheno_dataframe) <- pheno_dataframe$ID
pheno_dataframe$Group <- as.character(data_norm$Group)
pheno_dataframe$Age <- as.character(data_norm$Age)
control_group <- "Wild type"

source(paste0(function_path,"FUN.Update_data_for_association_test.R"))
new_data <- list()
limma_output <- Update_data_for_association_test(data_list = new_data,data_norm = data_norm,
                                         phenodata=pheno_dataframe,
                                         group_col="Group",control_group=control_group,
                                         batch_col="Age",type="microarray",
                                         current_gse = current_gse,
                                         function_path = function_path)

#C) Add gene symbols and sort by adjusted p-value
library(annotate)
library(org.Mm.eg.db)
for (i in names(limma_output)){
  limma_output[[i]]$Gene_symbol <- getSYMBOL(rownames(limma_output[[i]]),data="org.Mm.eg")[rownames(limma_output[[i]])]
  limma_output[[i]] <- limma_output[[i]][order(limma_output[[i]]$adj.P.Val,decreasing = F),]
}
head(limma_output[[1]])

FDR_thresh <- 0.05
sum(limma_output[[1]]$adj.P.Val<FDR_thresh)
#373 genes are differentially expressed in response to MYC deficiency 
#regardless of age (both in young and old mice)


##################################################
#Homework exercise:
#1) Identify differentially expressed genes between Wild type and MYC deficiency 
#separately for young and old mice

#2) Check if the overlap of diff expressed genes (adjusted p < 0.05) is 
#significant, using hypergeometric test

#3) Compare the overlap of diff expressed genes (adjusted p < 0.05) calculated 
#separately for young and old mice with diff expressed genes (adjusted p < 0.05)
#calculated for both ages simultaneously. Which set of genes is larger?