#plot venn diagram for significant up/down regulated genes
plot_venn <- function(SM, SH){
    upm <- rownames(SM[(SM$FDR < 0.05) & (SM$logFC > 0.),])
    uph <- rownames(SH[(SH$FDR < 0.05) & (SH$logFC > 0.),])
    dwm <- rownames(SM[(SM$FDR < 0.05) & (SM$logFC < 0.),])
    dwh <- rownames(SH[(SH$FDR < 0.05) & (SH$logFC < 0.),])

    library("ggvenn")

    list_up <- list(Mouse = upm,
                    Human = uph)
    list_dw <- list(Mouse = dwm,
                    Human = dwh)                

    gup <- ggvenn(list_up, c("Mouse", "Human"), 
                    fill_color = c("yellow", "blue"),
                    set_name_size=10, text_size=5) + 
                ggtitle('Upregulated') + 
                theme(plot.title = element_text(size=24, hjust = 0.5, face='bold', colour='red'))

    gdw <- ggvenn(list_dw, c("Mouse", "Human"), 
                    fill_color = c("yellow", "blue"),
                    set_name_size=10, text_size=5) + 
                ggtitle('Downregulated') + 
                theme(plot.title = element_text(size=24, hjust = 0.5, face='bold', colour='blue'))   

    venn <- cowplot::plot_grid(gup, gdw, labels = c("D", ""), nrow=1, label_size = 32)  

    return(venn)
}