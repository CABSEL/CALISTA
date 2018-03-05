visualization<-function(reduction,DATA,Results){
  writeLines('\nPlotting...\n')
  expected_clusters=Results$expected_clusters
  final_groups=Results$final_groups
  totDATA=DATA$totDATA^(1/2)
  switch (reduction,
          {
            initial_dims=DATA$numGENES
            perplexity=30
            no_dims=3
            dataANALYSING=totDATA
            score3=tsne(dataANALYSING,k=no_dims,initial_dims = initial_dims,perplexity = perplexity)
          },
          {score3_list=princomp(scale(totDATA))
           score3=score3_list$scores
           }, # this is just exactly opposite in  matlab 
          {
            no_dims=4
            t=10
            S=50
            writeLines("I am sorry, but I don't find any diffusionmap algorithm in R to reduce dimention")
            
          }
  )
  ###colormap for the original cell info
  #jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
  #                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  jet.colors <-colorRampPalette(c('#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                                  '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))

  colorMARK_time=jet.colors(DATA$nc)
  #colorMARK_time=terrain.colors(DATA$nc)
  
  #################################
  #if no time info
  #some information in find_progression2
  
  
  
  ################################
  #coloar map after clustering
  cluster_predicted=Results$cluster_predicted
  ClusterGroup_calista=final_groups
  #colorMARK_calista=topo.colors(length(cluster_predicted))
  jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  colorMARK_calista=jet.colors(length(cluster_predicted))
  Null_color='noncolor'
  c=data.frame(rep(Null_color,DATA$nvars),stringsAsFactors = FALSE)
  for (i in 1:expected_clusters) {
    idx_tmp=which(ClusterGroup_calista==Results$cluster_predicted[i])
    c[idx_tmp,]=rep(colorMARK_calista[Results$cluster_predicted[i]],length(idx_tmp))
  }
  #plot using par
  x11()
  #par(mfrow=c(1,2))
  xMinLab=min(score3[,1])
  xMaxLab=max(score3[,1])
  yMinLab=min(score3[,2])
  yMaxLab=max(score3[,2])
  zMinLab=min(score3[,3])
  zMaxLab=max(score3[,3])
  mtextlab=character()
  for(i in 1:DATA$num_time_points){
    scatterplot3d(score3[DATA$timeline==DATA$time[i],1],
                  score3[DATA$timeline==DATA$time[i],2],
                  score3[DATA$timeline==DATA$time[i],3],
                  xlim = c(xMinLab,xMaxLab),
                  ylim = c(yMinLab,yMaxLab),
                  zlim = c(zMinLab,zMaxLab),
                  xlab = "COMP1",
                  ylab = "COMP2",
                  zlab = "COMP3",
                  color = colorMARK_time[i],
                  pch=16)
    mtextlab=c(mtextlab,paste('time/stage',DATA$time[i]))
    par(new=TRUE)
  }
  title(main = 'original time/cell stage info(3D)')
  legend('topright',mtextlab, pch = 16,col = colorMARK_time[1:DATA$num_time_points])
  
#######################test######################
  x11(title = 'original time/cell stage info(2D)')
  gg_data=data.frame()
  label_color=c()
  label_names=c()
  for(k in 1:DATA$num_time_points){
    gg_data=rbind(gg_data,data.frame(x=score3[DATA$timeline==DATA$time[k],1],
                                     y=score3[DATA$timeline==DATA$time[k],2],
                                     z=score3[DATA$timeline==DATA$time[k],3],
                                     stage_colors=rep(colorMARK_time[k],length(score3[DATA$timeline==DATA$time[k],1]))
                                     #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1])),
                                     #gg_names=rep(paste('cluster',k),length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1]))
    ))
    label_color=c(label_color,colorMARK_time[k])
    label_names=c(label_names,paste('time/stage',DATA$time[k]))
  }
  names(label_color)=label_names
  p1<-ggplot(gg_data,mapping = aes(x=x,y=y))
  p1<-p1+geom_point(data = gg_data,color=gg_data[,"stage_colors"])
  #p1<-p1+geom_point(data = gg_data,aes(color=stage_colors))
  #p1<-p1+scale_colour_manual(values =label_color )
  p1<-p1+labs(x='COMP1',y='COMP2')

  p2<-ggplot(gg_data,mapping = aes(x=y,y=z))
  p2<-p2+geom_point(data = gg_data,color=gg_data[,"stage_colors"])
  p2<-p2+labs(x='COMP2',y='COMP3')
  
  p3<-ggplot(gg_data,mapping = aes(x=x,y=z))
  p3<-p3+geom_point(data = gg_data,color=gg_data[,"stage_colors"])
  p3<-p3+labs(x='COMP1',y='COMP3')
  
  grid.arrange(p1,p2,p3)
###########test###################  
  
  
  x11(title="1000")
  Results$figure1000=dev.cur()
  
  for (i in 1:length(cluster_predicted)) {
    scatterplot3d(score3[ClusterGroup_calista==cluster_predicted[i],1],
                  score3[ClusterGroup_calista==cluster_predicted[i],2],
                  score3[ClusterGroup_calista==cluster_predicted[i],3],
                  xlim = c(xMinLab,xMaxLab),
                  ylim = c(yMinLab,yMaxLab),
                  zlim = c(zMinLab,zMaxLab),
                  xlab = "PC1",
                  ylab = "PC2",
                  zlab = "PC3",
                  color = colorMARK_calista[i],
                  pch=16)
    par(new=TRUE)
  }
  title(main = 'Cell Clustering(3D)')
  legend('topright',unlist(Results$legendInfo_calista), pch = 16,col = colorMARK_calista[1:length(cluster_predicted)])
  
  #######################test######################
  x11(title = 'Cell Clustering(2D)')
  Results$figure1001=dev.cur()
  gg_data=data.frame()
  label_color=c()
  label_names=c()
  for(i in 1:length(cluster_predicted)){
    gg_data=rbind(gg_data,data.frame(x=score3[ClusterGroup_calista==cluster_predicted[i],1],
                                     y=score3[ClusterGroup_calista==cluster_predicted[i],2],
                                     z=score3[ClusterGroup_calista==cluster_predicted[i],3],
                                     stage_colors=rep(colorMARK_calista[i],length(score3[ClusterGroup_calista==cluster_predicted[i],1]))
                                     #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1])),
                                     #gg_names=rep(paste('cluster',k),length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1]))
    ))
  }
  
  #plot(gg_data[,1],gg_data[,2],type='p',
  #     col=gg_data[,4])
  p1<-ggplot(gg_data,mapping = aes(x=x,y=y))
  p1<-p1+geom_point(data = gg_data,color=gg_data[,"stage_colors"])
  p1<-p1+labs(x='PC1',y='PC2')
  
  p2<-ggplot(gg_data,mapping = aes(x=y,y=z))
  p2<-p2+geom_point(data = gg_data,color=gg_data[,"stage_colors"])
  p2<-p2+labs(x='PC2',y='PC3')
  
  p3<-ggplot(gg_data,mapping = aes(x=x,y=z))
  p3<-p3+geom_point(data = gg_data,color=gg_data[,"stage_colors"])
  p3<-p3+labs(x='PC1',y='PC3')
  
  grid.arrange(p1,p2,p3)
  ###########test###################  
  
  Results$colorMARK_calista=colorMARK_calista
  Results$colorMARK_time=colorMARK_time
  Results$score3=score3
  Results$c=c
  
  return(Results)
}





