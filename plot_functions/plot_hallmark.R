#plot hallmark ontologies sorted by dendrogram
plot_hallmark <- function(gsea, orderRow=TRUE, col_order=c(), fontsize=18, return_table=FALSE){
    tmp <- gsea[[1]]
    #terms_list <- rownames(tmp[grepl("HALLMARK", rownames(tmp)),])
    terms_list <- c('HALLMARK_ANDROGEN_RESPONSE',
                    'HALLMARK_ANGIOGENESIS',
                    'HALLMARK_APICAL_JUNCTION',
                    'HALLMARK_APOPTOSIS',
                    'HALLMARK_BILE_ACID_METABOLISM',
                    'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
                    'HALLMARK_DNA_REPAIR',
                    'HALLMARK_E2F_TARGETS',
                    'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
                    #'HALLMARK_ESTROGEN_RESPONSE_EARLY',
                    'HALLMARK_FATTY_ACID_METABOLISM',
                    'HALLMARK_G2M_CHECKPOINT',
                    'HALLMARK_HEME_METABOLISM',
                    'HALLMARK_HYPOXIA',
                    'HALLMARK_IL6_JAK_STAT3_SIGNALING',
                    'HALLMARK_INFLAMMATORY_RESPONSE',
                    'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                    'HALLMARK_KRAS_SIGNALING_UP',
                    #'HALLMARK_MITOTIC_SPINDLE',
                    'HALLMARK_MYC_TARGETS_V1',
                    'HALLMARK_MYOGENESIS',
                    'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                    'HALLMARK_P53_PATHWAY',
                    'HALLMARK_PI3K_AKT_MTOR_SIGNALING',
                    'HALLMARK_PROTEIN_SECRETION',
                    'HALLMARK_SPERMATOGENESIS',
                    'HALLMARK_TGF_BETA_SIGNALING',
                    'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                    'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
                    'HALLMARK_UV_RESPONSE_DN',
                    'HALLMARK_XENOBIOTIC_METABOLISM',
                    'KEGG_COMPLEMENT_AND_COAGULATION_CASCADES',
                    'REACTOME_BASE_EXCISION_REPAIR',
                    'GOBP_ATP_METABOLIC_PROCESS',
                    'KEGG_CITRATE_CYCLE_TCA_CYCLE',
                    'GOBP_MITOCHONDRIAL_TRANSLATION',
                    'KEGG_RIBOSOME',
                    'REACTOME_MITOCHONDRIAL_BIOGENESIS'
                    )
    #within-group geometric pvalue
    SLR <- list()
    SLI <- list()
    SLA <- list()
    cols <- c('pathway', 'pval', 'NES')

    for (col in cols){
        SLR[[col]] <- data.frame(row.names=gsea[['Reprogramming:Mouse']]$pathway)
        SLI[[col]] <- data.frame(row.names=gsea[['Interventions:Median_lifespan']]$pathway)
        SLA[[col]] <- data.frame(row.names=gsea[['Aging:Brain']]$pathway)
        for (name in names(gsea)){
            tmp <- gsea[[name]][col]
            rownames(tmp) <- gsea[[name]]$pathway
            colnames(tmp) <- c(name)
            if (grepl("Aging", name, fixed=T)){    
                SLA[[col]] <- transform(merge(SLA[[col]], tmp, 
                                                    by=0, sort=F,), 
                                                    row.names=Row.names, Row.names=NULL)
            }else if (grepl("Interventions", name, fixed=T)){
                SLI[[col]] <- transform(merge(SLI[[col]], tmp, 
                                                    by=0, sort=F,), 
                                                    row.names=Row.names, Row.names=NULL)
            }else{
                SLR[[col]] <- transform(merge(SLR[[col]], tmp, 
                                        by=0, sort=F,), 
                                        row.names=Row.names, Row.names=NULL) 
            }
        }
    }
    #
    slr <- data.frame(apply(SLR[['pval']], 1, function(x) exp(mean(log(x)))))
    colnames(slr) <- c('P.Value')
    slr$FDR <- NaN
    #
    sli <- data.frame(apply(SLI[['pval']], 1, function(x) exp(mean(log(x)))))
    colnames(sli) <- c('P.Value')
    sli$FDR <- NaN
    #
    sla <- data.frame(apply(SLA[['pval']], 1, function(x) exp(mean(log(x)))))
    colnames(sla) <- c('P.Value')
    sla$FDR <- NaN
    
    #select by gmean criterion
    ff <- gsea
    df <- c()
    for (name in names(ff)){
        fgs <- as.data.frame(ff[[name]])
        fgs <- fgs[terms_list,]
        tmp <- cbind(fgs[c("pathway", "NES", "size", "pval", "padj")], rep(name, nrow(fgs)))
        rownames(tmp) <- NULL
        colnames(tmp) <- c("Function", "NES", "Size", "P.val", "P.adj", "Signature")
        df <- rbind(df, tmp)
    }

    slr <- slr[terms_list,]
    sli <- sli[terms_list,]
    sla <- sla[terms_list,]
    
    slr$FDR <- p.adjust(slr$P.Value, method="BH")
    sli$FDR <- p.adjust(sli$P.Value, method="BH")
    sla$FDR <- p.adjust(sla$P.Value, method="BH")

    df <- df[complete.cases(df),]
    gmean <- function(x) exp(mean(log(x)))
    passed <- c()
    for (t in df$Function){

        rp <- slr[t,]$FDR < 0.05
        ip <- sli[t,]$FDR < 0.05
        ap <- sla[t,]$FDR < 0.05
        
        if (rp  
            #& any(c(ip, ap))
            ){
            passed <- c(passed, t)
        }
    }
    passed <- unique(passed)

    df <- df[df$Function %in% passed,]

    #function to first letter up
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        return(x)
    }
    df$Function <- gsub("_", " ", df$Function)
    df$Function <- gsub("HALLMARK ", "", df$Function)
    df$Function <- firstup(tolower(df$Function))
    df$Function <- gsub("Dna repair", "DNA repair", df$Function)
    df$Function <- gsub("E2f targets", "E2F targets", df$Function)
    df$Function <- gsub('Estrogen response early', "Estrogen response", df$Function)
    df$Function <- gsub('G2m checkpoint', "G2M checkpoint", df$Function)
    df$Function <- gsub('Il6 jak stat3 signaling', 'IL6 JAK STAT3 signaling', df$Function)
    df$Function <- gsub('Kras signaling up', 'KRAS signaling up', df$Function)
    df$Function <- gsub('Myc targets v1', 'MYC targets', df$Function)
    df$Function <- gsub('Pi3k akt mtor signaling', 'PI3K AKT mTOR signaling', df$Function)
    df$Function <- gsub('Tgf beta signaling', 'TGF beta signaling', df$Function)
    df$Function <- gsub('Tnfa signaling via nfkb', 'TNFA signaling via NFkB', df$Function)
    df$Function <- gsub('Uv response dn', 'UV response down', df$Function)
    df$Function <- gsub('Kegg complement and coagulation cascades', 'Complement and coagulation (KEGG)', df$Function)
    df$Function <- gsub('Reactome base excision repair', 'Base excision repair (Reactome)', df$Function)
    df$Function <- gsub('Gobp atp metabolic process', 'ATP metabolic process (GO:BP)', df$Function)
    df$Function <- gsub('Gobp atp metabolic process', 'ATP metabolic process (GO:BP)', df$Function)
    df$Function <- gsub('Kegg citrate cycle tca cycle', 'Citrate cycle TCA cycle (KEGG)', df$Function)
    df$Function <- gsub('Gobp mitochondrial translation', 'Mitochondrial translation (GO:BP)', df$Function)
    df$Function <- gsub('Kegg ribosome', 'Ribosome (KEGG)', df$Function)
    
    if (return_table == TRUE){
        return(list('data'=df, 'slr'=slr, 'sli'=sli, 'sla'=sla))
    }

    predist <- reshape2::dcast(data = df, formula = Function~Signature, fun.aggregate = sum, value.var = "NES")
    predist <- as.matrix(predist)
    
    library(amap)
    # Cluster rows
    if (orderRow) {
        dd.row <- as.dendrogram(hclust(dist(predist, method = "euclidean"), method = "complete"))
        row.ord <- order.dendrogram(dd.row)
        ordered_row_names <- predist[row.ord, 'Function']
        df$Function <- factor(df$Function, levels=ordered_row_names)
        if (length(col_order) != 0){
            df$Signature <- factor(df$Signature, levels=col_order)
        }
    }

    # df$Signature <- gsub("\\$", ":", df$Signature)
    # df$Signature <- gsub("_", " ", df$Signature)
    # df$Signature <- gsub("\\+dox_mef", "", df$Signature)
    # df$Signature <- gsub(" liver", "", df$Signature)
    # df$Signature <- gsub("All", "Global", df$Signature)

    g <- ggplot(df, aes(x=Signature, y=Function, fill=NES)) + 
        geom_tile(colour = "black",size=0.1)+
        scale_fill_gradient2(high="red3",low="blue3",
                                #mid=background_color,
                                guide = "colorbar",
                                #limits=c(-max_abs_score, max_abs_score),
                                ) +
        labs(x = "",y = "", title="") +
        # scale_x_discrete(expand = c(0, 0)) + 
        # scale_y_discrete(expand = c(0, 0)) +
        scale_y_discrete(position='right') +
        #theme_minimal() + 
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=fontsize, colour="black"),
            axis.text.y = element_text(size=fontsize, colour="black"),
            panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
            legend.key.width = unit(1.5, "cm"),
            legend.title=element_text(size=fontsize, face = "bold"),
            legend.text=element_text(size=fontsize),
            legend.position="top",
            panel.grid.major = element_blank(), #element_line(colour="black"), 
            panel.grid.minor = element_blank(),# = element_line(colour="black"),
            plot.margin = unit(c(-0.6, 0.,0.,0.5), "cm"))
            
    return(g)

}
