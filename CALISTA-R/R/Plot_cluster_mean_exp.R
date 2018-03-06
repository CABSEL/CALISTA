Plot_cluster_mean_exp<-function(Results,DATA){
  writeLines('\nPlotting mean gene expressions...\n')
  p1=sqrt(DATA$numGENES)
  D=list()
  MEAN=matrix(0,ncol(DATA$totDATA),Results$expected_clusters)
  p_ggplot=list()
  for(i in 1:Results$expected_clusters){
    D[i]=list(t(DATA$totDATA[which(Results$final_groups==i),]))
    MEAN[,i]=rowMeans(D[[i]])
  }
  
  #x11(title = 'Plotting  mean gene expression')
  #par(mfrow=c(p1,p1))
  for(j in 1:DATA$numGENES){
    gg_data=data.frame()
    label_color=c()
    label_names=c()
    for(i in 1:Results$expected_clusters){
        gg_data=rbind(gg_data,data.frame(x=Results$cluster_progression[i],
                                         y=MEAN[j,i],
                                         gg_names=paste('cluster',i),
                                         gg_color=Results$colorMARK_calista[i]))
          label_color=c(label_color,Results$colorMARK_calista[i])
          label_names=c(label_names,paste('cluster',i))
    }
    
    names(label_color)=label_names
    p<-ggplot(gg_data)
    p<-p+geom_point(mapping = aes(x=x,y=y,color=gg_names),data = gg_data)
    #p<-p+labs(x='PC1',y='PC2')
    p<-p+scale_color_manual(name='cluster',
                           #breaks=label_names,
                           values = label_color,
                           labels=label_names)
    hh=Results$TRANSITION$final_graph
    V(hh)$size=1
    V(hh)$color=Results$colorMARK_calista
    E(hh)$width=2.5
    h_after_layout=data.frame(Results$cluster_progression,MEAN[j,],label_names,stringsAsFactors = FALSE)
    colnames(h_after_layout)=c('x_center','y_center','label_name')
    p<-p+geom_point(mapping =aes(x=x_center,y=y_center,color=label_name),data = h_after_layout)
    h_edge=get.edgelist(hh)
    segment_data=data.frame(x=h_after_layout[h_edge[,1],1],
                            y=h_after_layout[h_edge[,1],2],
                            xend=h_after_layout[h_edge[,2],1],
                            yend=h_after_layout[h_edge[,2],2])
    p<-p+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend))
    #p_title=DATA$genes[j]
    #p<p+ggtitle(paste('gene',p_title))
    p_ggplot[[j]]=p
    #par(new=TRUE)
  }
  #grid.arrange(grobs=p_ggplot)
  #x11()
  ml=marrangeGrob(grobs = p_ggplot,ncol = 2,nrow = 2)
  ggsave("ploting  meean gene expression.pdf",ml)
  writeLines('\n finished! You can find it in "ploting  meean gene expression.pdf" in current working directory  \n')
  Results$singleCellClusterDATA=D
  return(Results)

  
}