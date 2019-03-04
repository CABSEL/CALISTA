#' A calista Function
#'
#' @param reduction,DATA,Results
#' @keywords calista
#' @export
#' @examples
#' visualization()



visualization<-function(reduction,DATA,Results){
  #writeLines('\nPlotting...\n')
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
  #display.brewer.all() #display all palettes
  colorMARK_time = brewer.pal(11, "Spectral")
  colorMARK_time = colorMARK_time[-6] #remove white colors
  if (DATA$nc<=length(colorMARK_time)){
    colorMARK_time=colorMARK_time[1:DATA$nc]
  }else{
    # I can add more tones to this palette :
    colorMARK_time = colorRampPalette(colorMARK_time)(DATA$nc)
    # Plot it
    pie(rep(1, length(colorMARK_time)), col = colorMARK_time , main="")
  }

  #################################
  #if no time info
  #some information in find_progression2



  ################################
  #colormap after clustering
  cluster_predicted=Results$cluster_predicted
  ClusterGroup_calista=final_groups

  colorMARK_calista = brewer.pal(8, "Set2")
  #colorMARK_calista = colorMARK_calista[-6] #remove white colors
  if (DATA$nc<=length(colorMARK_calista)){
    colorMARK_calista=colorMARK_calista[1:length(cluster_predicted)]
  }else{
    # I can add more tones to this palette :
    colorMARK_calista = colorRampPalette(colorMARK_calista)(length(cluster_predicted))
    # Plot it
    #pie(rep(1, length(colorMARK_calista)), col = colorMARK_calista , main="")
  }

  Null_color='noncolor'
  c=data.frame(rep(Null_color,DATA$nvars),stringsAsFactors = FALSE)
  for (i in 1:expected_clusters) {
    idx_tmp=which(ClusterGroup_calista==Results$cluster_predicted[i])
    c[idx_tmp,]=rep(colorMARK_calista[Results$cluster_predicted[i]],length(idx_tmp))
  }

  xMinLab=min(score3[,1])
  xMaxLab=max(score3[,1])
  yMinLab=min(score3[,2])
  yMaxLab=max(score3[,2])
  zMinLab=min(score3[,3])
  zMaxLab=max(score3[,3])

  ## Plot time info
  plotly_time=data.frame()
  label_color=c()
  label_names=c()
  for(k in 1:DATA$num_time_points){
    plotly_time=rbind(plotly_time,data.frame(x=score3[DATA$timeline==DATA$time[k],1],
                                     y=score3[DATA$timeline==DATA$time[k],2],
                                     z=score3[DATA$timeline==DATA$time[k],3],
                                     stage_colors=rep(colorMARK_time[k],length(score3[DATA$timeline==DATA$time[k],1])),
                                     #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1])),
                                     gg_names=rep(paste('time/stage',DATA$time[k]),length(score3[DATA$timeline==DATA$time[k],1]))
    ))
    label_color=c(label_color,colorMARK_time[k])
    label_names=c(label_names,paste('time/stage',DATA$time[k]))
  }
  names(label_color)=label_names

  p <- plot_ly(plotly_time, x = ~x, y = ~y, z = ~z, color = ~gg_names, colors = colorMARK_time,
               size=1, alpha=0.8) %>%
    add_markers() %>%
    layout(title = 'original time/cell stage info',
           scene = list(camera = list(eye = list(x = -0.1, y = -2.5, z = 0.9)),
                        xaxis = list(title = 'COMP1',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     range = c(xMinLab,xMaxLab),
                                     #type = 'log',
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwidth = 2),
                        yaxis = list(title = 'COMP2',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     range = c(yMinLab,yMaxLab),
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwith = 2),
                        zaxis = list(title = 'COMP3',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     #type = 'log',
                                     range = c(zMinLab,zMaxLab),
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwith = 2)),
           paper_bgcolor = 'rgb(243, 243, 243)',
           plot_bgcolor = 'rgb(243, 243, 243)')
  Results$p_time=p


  ## Plot cell clustering results
  Results$figure1000=dev.cur()
  plotly_cluster=data.frame()
  label_color=c()
  label_names=c()
  for(i in 1:length(cluster_predicted)){
    plotly_cluster=rbind(plotly_cluster,data.frame(x=score3[ClusterGroup_calista==cluster_predicted[i],1],
                                     y=score3[ClusterGroup_calista==cluster_predicted[i],2],
                                     z=score3[ClusterGroup_calista==cluster_predicted[i],3],
                                     stage_colors=rep(colorMARK_calista[i],length(score3[ClusterGroup_calista==cluster_predicted[i],1])),
                                     #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1])),
                                     gg_names=rep(paste('cluster',i),length(score3[Results$final_groups==Results$cluster_predicted[i],1]))
    ))
  }

  names(label_color)=label_names

  p <- plot_ly(plotly_cluster, x = ~x, y = ~y, z = ~z, color = ~gg_names, colors = colorMARK_calista,
               size=1, alpha=0.8) %>%
    add_markers() %>%
    layout(title = 'Cell Clustering',
      scene = list(camera = list(eye = list(x = -0.1, y = -2.5, z = 0.9)),
                       xaxis = list(title = 'COMP1',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     range = c(xMinLab,xMaxLab),
                                     #type = 'log',
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwidth = 2),
                        yaxis = list(title = 'COMP2',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     range = c(yMinLab,yMaxLab),
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwith = 2),
                        zaxis = list(title = 'COMP3',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     #type = 'log',
                                     range = c(zMinLab,zMaxLab),
                                     zerolinewidth = 1,
                                     ticklen = 5,
                                     gridwith = 2)),
           paper_bgcolor = 'rgb(243, 243, 243)',
           plot_bgcolor = 'rgb(243, 243, 243)')
  Results$p_cluster=p


  ###########

  Results$colorMARK_calista=colorMARK_calista
  Results$colorMARK_time=colorMARK_time
  Results$score3=score3
  Results$c=c
  Results$plotly_cluster=plotly_cluster
  Results$plotly_time=plotly_time

  return(Results)
}





