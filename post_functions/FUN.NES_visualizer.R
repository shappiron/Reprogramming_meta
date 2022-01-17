#####################################################################
#This function visualizes NES scores for certain signature

NES_visualizer <- function(test_output,signature_name,intervention_name){
#  library(ggpubr)
  if (length(test_output$NES_vector[[signature_name]][[intervention_name]])==2){
    temp_vector_up <- test_output$NES_vector[[signature_name]][[intervention_name]]$Up
    temp_vector_down <- test_output$NES_vector[[signature_name]][[intervention_name]]$Down
    NES_data <- list(Up=data.frame(Index=c(1:length(temp_vector_up)),
                                   NES=temp_vector_up,Direction="Upregulated"),
                     Down=data.frame(Index=c(1:length(temp_vector_down)),
                                     NES=temp_vector_down,Direction="Downregulated"))
    
    NES_data$Up$Gene_position = ifelse(NES_data$Up$Index %in% test_output$Gene_table[[signature_name]][[intervention_name]]$Up$Index,1,0)
    NES_data$Down$Gene_position = ifelse(NES_data$Down$Index %in% test_output$Gene_table[[signature_name]][[intervention_name]]$Down$Index,1,0)
    NES_data$Up$Index <- NES_data$Up$Index/length(temp_vector_up)
    NES_data$Down$Index <- NES_data$Down$Index/length(temp_vector_down)
    visual_NES <- rbind(NES_data$Up,NES_data$Down)

    visual_NES_max_up <- data.frame(x=visual_NES[visual_NES$Direction=="Upregulated",]$Index[which.max(abs(visual_NES[visual_NES$Direction=="Upregulated",]$NES))],
                                    y=visual_NES[visual_NES$Direction=="Upregulated",]$NES[which.max(abs(visual_NES[visual_NES$Direction=="Upregulated",]$NES))])
    visual_NES_max_down <- data.frame(x=visual_NES[visual_NES$Direction=="Downregulated",]$Index[which.max(abs(visual_NES[visual_NES$Direction=="Downregulated",]$NES))],
                                    y=visual_NES[visual_NES$Direction=="Downregulated",]$NES[which.max(abs(visual_NES[visual_NES$Direction=="Downregulated",]$NES))])

    segment_data_up = data.frame(x = visual_NES[visual_NES$Direction=="Upregulated" & visual_NES$Gene_position==1,]$Index,
                                   xend = visual_NES[visual_NES$Direction=="Upregulated" & visual_NES$Gene_position==1,]$Index,y = 0,
                                   yend =visual_NES[visual_NES$Direction=="Upregulated" & visual_NES$Gene_position==1,]$NES)

    segment_data_down = data.frame(x = visual_NES[visual_NES$Direction=="Downregulated" & visual_NES$Gene_position==1,]$Index,
                                   xend = visual_NES[visual_NES$Direction=="Downregulated" & visual_NES$Gene_position==1,]$Index,y = 0,
                                   yend =visual_NES[visual_NES$Direction=="Downregulated" & visual_NES$Gene_position==1,]$NES)
    visual_NES$Direction <- factor(visual_NES$Direction,levels=c("Upregulated","Downregulated"))

    p1 <- ggplot(visual_NES, aes(Index,NES,color=Direction))+
      geom_segment(data=visual_NES_max_up,aes(x=x,xend=x,
                                              y=y,yend=0),lwd=1.2,lty=2,col="red3",alpha=0.8)+
      geom_segment(data=visual_NES_max_down,aes(x=x,xend=x,
                                                y=y,yend=0),lwd=1.2,lty=2,col="blue3",alpha=0.8)+
      geom_hline(yintercept = 0,lwd=0.5,lty=1)+
      geom_vline(xintercept = 0,lwd=0.2,lty=1)+
      geom_vline(xintercept = 1,lwd=0.2,lty=1)+
      geom_line(stat="identity",lwd=1,alpha=0.8)+
      labs(title = paste("(",intervention_name,") vs (",gsub("_"," ",signature_name),")",sep=""),x="",y="NES")+
      scale_colour_manual(values = c("red", "blue"))+
      theme_bw(base_size = 16)+
      geom_hline(yintercept = visual_NES_max_up$y,lwd=1,lty=3,color="red")+
      geom_hline(yintercept = visual_NES_max_down$y,lwd=1,lty=3,color="blue")+
      coord_cartesian(xlim=c(0,1))+
      theme(legend.position = c(0.87,0.82),
            legend.direction = "vertical",
            legend.title = element_blank(),
            legend.key.height = unit(1,"cm"),
            #              legend.spacing.x = unit(0.3,"cm"),
            #              legend.box.spacing = unit(0.7,"cm"),
            legend.text = element_text(size=18,colour="black"),
            plot.title = element_text(size=26,colour="black",face="bold",hjust = 0.5,
                                      margin=margin(0,0,18,0)),
            axis.text.x = element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.length.y = unit(0.1,"cm"),
            axis.text.y = element_text(size=14,colour="black"),
            axis.title.y =element_text(size=21,color="black",face="bold"),
            #              legend.key.size = unit(0.2, "cm"),
            #              legend.key.width = unit(0.7,"cm"),
            #legend.background = element_blank(),
            legend.background = element_rect(colour = "black",size=0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())+
      guides(color = guide_legend(override.aes = list(alpha = 1)))
    
    #ggbld <- ggplot_build(p1)
    #ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$y.major_source
    #temp <- ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$y.major_source
    
    p2 <- ggplot(segment_data_up)+
      geom_hline(yintercept = 0,col="white")+
      geom_vline(xintercept=0.5,lwd=1.2,lty=2,col="black",alpha=0.8)+
      geom_vline(xintercept = 0,lwd=0.2,lty=1)+
      geom_vline(xintercept = 1,lwd=0.2,lty=1)+
      geom_vline(aes(xintercept=xend),lwd=0.5,lty=1,col="red",alpha=0.7)+
      geom_vline(data=segment_data_down,aes(xintercept = xend),lwd=0.5,lty=1,col="blue",alpha=0.7)+
      geom_vline(xintercept = visual_NES_max_up$x,lwd=1.2,lty=2,col="red3",alpha=1)+
      geom_vline(xintercept = visual_NES_max_down$x,lwd=1.2,lty=2,col="blue3",alpha=1)+
      labs(x="Normalized index",y="NES")+
      coord_cartesian(xlim=c(0,1),ylim=c(-100,100))+
      #scale_y_continuous(breaks=temp)+
      theme_bw()+
      theme(axis.text.x = element_text(hjust = 0.5,size=14,color="black",
                                       margin = margin(b=4,t=3)),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.length.y = unit(0.1,"cm"),
            axis.title.x=element_text(size=21,color="black",face="bold"),
            axis.title.y =element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
  }else{
    if (names(test_output$NES_vector[[signature_name]][[intervention_name]])=="Up"){
      temp_vector_up <- test_output$NES_vector[[signature_name]][[intervention_name]]$Up
      NES_data <- list(Up=data.frame(Index=c(1:length(temp_vector_up)),
                                     NES=temp_vector_up,Direction="Upregulated"))
      
      NES_data$Up$Gene_position = ifelse(NES_data$Up$Index %in% test_output$Gene_table[[signature_name]][[intervention_name]]$Up$Index,1,0)
      NES_data$Up$Index <- NES_data$Up$Index/length(temp_vector_up)
      visual_NES <- NES_data$Up
      
      visual_NES_max_up <- data.frame(x=visual_NES[visual_NES$Direction=="Upregulated",]$Index[which.max(abs(visual_NES[visual_NES$Direction=="Upregulated",]$NES))],
                                      y=visual_NES[visual_NES$Direction=="Upregulated",]$NES[which.max(abs(visual_NES[visual_NES$Direction=="Upregulated",]$NES))])
      
      segment_data_up = data.frame(x = visual_NES[visual_NES$Direction=="Upregulated" & visual_NES$Gene_position==1,]$Index,
                                   xend = visual_NES[visual_NES$Direction=="Upregulated" & visual_NES$Gene_position==1,]$Index,y = 0,
                                   yend =visual_NES[visual_NES$Direction=="Upregulated" & visual_NES$Gene_position==1,]$NES)

      p1 <- ggplot(visual_NES, aes(Index,NES,color=Direction))+
        geom_segment(data=visual_NES_max_up,aes(x=x,xend=x,
                                                y=y,yend=0),lwd=1.2,lty=2,col="red3",alpha=0.8)+
        geom_hline(yintercept = 0,lwd=0.5,lty=1)+
        geom_vline(xintercept = 0,lwd=0.2,lty=1)+
        geom_vline(xintercept = 1,lwd=0.2,lty=1)+
        geom_line(stat="identity",lwd=1,alpha=0.8)+
        labs(title = paste("(",intervention_name,") vs (",gsub("_"," ",signature_name),")",sep=""),x="",y="NES")+
        scale_colour_manual(values = c("red"))+
        theme_bw(base_size = 16)+
        geom_hline(yintercept = visual_NES[visual_NES$Direction=="Upregulated",]$NES[which.max(abs(visual_NES[visual_NES$Direction=="Upregulated",]$NES))],lwd=1,lty=3,color="red")+
        coord_cartesian(xlim=c(0,1))+
        theme(legend.position = c(0.87,0.84),
              legend.direction = "vertical",
              legend.title = element_blank(),
              #legend.key.height = unit(0.3,"cm"),
              #              legend.spacing.x = unit(0.3,"cm"),
              #              legend.box.spacing = unit(0.7,"cm"),
              legend.text = element_text(size=18,colour="black"),
              plot.title = element_text(size=26,colour="black",face="bold",hjust = 0.5,
                                        margin=margin(0,0,18,0)),
              axis.text.x = element_blank(),
              axis.ticks.x=element_blank(),
              axis.ticks.length.y = unit(0.1,"cm"),
              axis.text.y = element_text(size=14,colour="black"),
              axis.title.y =element_text(size=21,color="black",face="bold"),
              #              legend.key.size = unit(0.2, "cm"),
              #              legend.key.width = unit(0.7,"cm"),
              legend.background = element_blank(),
              legend.box.background = element_rect(colour = "black"))+
        guides(color = guide_legend(override.aes = list(alpha = 1)))
      
      ggbld <- ggplot_build(p1)
      temp <- ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$y.major_source
      
      p2 <- ggplot(segment_data_up)+
        geom_hline(yintercept = 0,col="white")+
        geom_vline(xintercept=0.5,lwd=1.2,lty=2,col="black",alpha=0.8)+
        geom_vline(xintercept = 0,lwd=0.2,lty=1)+
        geom_vline(xintercept = 1,lwd=0.2,lty=1)+
        geom_vline(aes(xintercept=xend),lwd=0.7,lty=1,col="red",alpha=0.7)+
        geom_vline(xintercept = visual_NES_max_up$x,lwd=1.2,lty=2,col="red3",alpha=1)+
        labs(x="Normalized index",y="NES")+
        coord_cartesian(xlim=c(0,1),ylim=c(-100,100))+
        scale_y_continuous(breaks=temp)+
        theme_bw()+
        theme(axis.text.x = element_text(hjust = 0.5,size=14,color="black",
                                         margin = margin(b=4,t=3)),
              axis.text.y = element_text(size=14,colour="white"),
              axis.ticks.y = element_line(colour="white"),
              axis.ticks.length.y = unit(0.1,"cm"),
              axis.title.x=element_text(size=21,color="black",face="bold"),
              axis.title.y =element_text(size=21,color="white",face="bold"))
      
    }else{
      temp_vector_down <- test_output$NES_vector[[signature_name]][[intervention_name]]$Down
      NES_data <- list(Down=data.frame(Index=c(1:length(temp_vector_down)),
                                       NES=temp_vector_down,Direction="Downregulated"))
      
      NES_data$Down$Gene_position = ifelse(NES_data$Down$Index %in% test_output$Gene_table[[signature_name]][[intervention_name]]$Down$Index,1,0)
      NES_data$Down$Index <- NES_data$Down$Index/length(temp_vector_down)
      visual_NES <- NES_data$Down
      visual_NES_max_down <- data.frame(x=visual_NES[visual_NES$Direction=="Downregulated",]$Index[which.max(abs(visual_NES[visual_NES$Direction=="Downregulated",]$NES))],
                                        y=visual_NES[visual_NES$Direction=="Downregulated",]$NES[which.max(abs(visual_NES[visual_NES$Direction=="Downregulated",]$NES))])
      
      segment_data_down = data.frame(x = visual_NES[visual_NES$Direction=="Downregulated" & visual_NES$Gene_position==1,]$Index,
                                     xend = visual_NES[visual_NES$Direction=="Downregulated" & visual_NES$Gene_position==1,]$Index,y = 0,
                                     yend =visual_NES[visual_NES$Direction=="Downregulated" & visual_NES$Gene_position==1,]$NES)

      p1 <- ggplot(visual_NES, aes(Index,NES,color=Direction))+
        geom_segment(data=visual_NES_max_down,aes(x=x,xend=x,
                                                  y=y,yend=0),lwd=1.2,lty=2,col="blue3",alpha=0.8)+
        geom_hline(yintercept = 0,lwd=0.5,lty=1)+
        geom_vline(xintercept = 0,lwd=0.2,lty=1)+
        geom_vline(xintercept = 1,lwd=0.2,lty=1)+
        geom_line(stat="identity",lwd=1,alpha=0.8)+
        labs(title = paste("(",intervention_name,") vs (",gsub("_"," ",signature_name),")",sep=""),x="",y="NES")+
        scale_colour_manual(values = c("blue"))+
        theme_bw(base_size = 16)+
        geom_hline(yintercept = visual_NES[visual_NES$Direction=="Downregulated",]$NES[which.max(abs(visual_NES[visual_NES$Direction=="Downregulated",]$NES))],lwd=1,lty=3,color="blue")+
        coord_cartesian(xlim=c(0,1))+
        theme(legend.position = c(0.87,0.84),
              legend.direction = "vertical",
              legend.title = element_blank(),
              #legend.key.height = unit(0.3,"cm"),
              #              legend.spacing.x = unit(0.3,"cm"),
              #              legend.box.spacing = unit(0.7,"cm"),
              legend.text = element_text(size=18,colour="black"),
              plot.title = element_text(size=26,colour="black",face="bold",hjust = 0.5,
                                        margin=margin(0,0,18,0)),
              axis.text.x = element_blank(),
              axis.ticks.x=element_blank(),
              axis.ticks.length.y = unit(0.1,"cm"),
              axis.text.y = element_text(size=14,colour="black"),
              axis.title.y =element_text(size=21,color="black",face="bold"),
              #              legend.key.size = unit(0.2, "cm"),
              #              legend.key.width = unit(0.7,"cm"),
              legend.background = element_blank(),
              legend.box.background = element_rect(colour = "black"))+
        guides(color = guide_legend(override.aes = list(alpha = 1)))
      
      ggbld <- ggplot_build(p1)
      temp <- ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$y.major_source
      
      p2 <- ggplot(segment_data_down)+
        geom_hline(yintercept = 0,col="white")+
        geom_vline(xintercept=0.5,lwd=1.2,lty=2,col="black",alpha=0.8)+
        geom_vline(xintercept = 0,lwd=0.2,lty=1)+
        geom_vline(xintercept = 1,lwd=0.2,lty=1)+
        geom_vline(aes(xintercept = xend),lwd=0.7,lty=1,col="blue",alpha=0.7)+
        geom_vline(xintercept = visual_NES_max_down$x,lwd=1.2,lty=2,col="blue3",alpha=1)+
        labs(x="Normalized index",y="NES")+
        coord_cartesian(xlim=c(0,1),ylim=c(-100,100))+
        scale_y_continuous(breaks=temp)+
        theme_bw()+
        theme(axis.text.x = element_text(hjust = 0.5,size=14,color="black",
                                         margin = margin(b=4,t=3)),
              axis.text.y = element_text(size=14,colour="white"),
              axis.ticks.y = element_line(colour="white"),
              axis.ticks.length.y = unit(0.1,"cm"),
              axis.title.x=element_text(size=21,color="black",face="bold"),
              axis.title.y =element_text(size=21,color="white",face="bold"))
      
    }
  }
#  p3 <- ggarrange(p1, p2, heights = c(2, 0.7),
#                  ncol = 1, nrow = 2, align = "v")
  
    return(list(p1=p1,p2=p2))
}