rej_genes_enrich <- function(enr, fontsize=18){
    sel <- c(
    'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
    'REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS',
    'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION',
    'REACTOME_COLLAGEN_FORMATION',
    'GOBP_CRISTAE_FORMATION',
    'KEGG_FOCAL_ADHESION'
     )
    enr <- enr[enr$ID %in% sel,]

    enr$ID <- gsub('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', "EMT", enr$ID)
    enr$ID <- gsub('REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS', "Integrin cell\ninteractions", enr$ID)
    enr$ID <- gsub('REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION', "ECM\norganization", enr$ID)
    enr$ID <- gsub('REACTOME_COLLAGEN_FORMATION', "Collagen\nformation", enr$ID)
    enr$ID <- gsub('GOBP_CRISTAE_FORMATION', "Cristae\nformation", enr$ID)
    enr$ID <- gsub('KEGG_FOCAL_ADHESION', "Focal\nadhesion", enr$ID)

    enr <- enr[order(enr$p.adjust),]
    enr$plog10 <- -log10(enr$p.adjust)

    enr$ID <- factor(enr$ID, levels=rev(enr$ID))

    p <-ggplot(enr, aes(x=plog10, y=ID)) +   
        geom_bar(stat="identity", colour="black", fill='firebrick')+
        geom_vline(xintercept=-log10(0.1), col="black", size=0.7, linestyle='dashed') +
        theme_minimal()+
        theme(axis.text=element_text(size=fontsize, hjust=1, color='black'),
                axis.title=element_text(size=fontsize, face="bold"), 
                plot.title=element_text(size=fontsize, face = "bold"),
                legend.title=element_text(size=fontsize, face = "bold"),
                legend.text=element_text(size=fontsize),
                legend.position="top",
                panel.background = element_rect(fill='white', colour='black', size=1.5, linetype='solid'),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.margin = unit(c(1.0,0.3,0,0.), "cm")) +
        labs(x = "-log10(adj.P)", y = "")
    return(p)
}