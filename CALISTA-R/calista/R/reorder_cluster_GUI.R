#' A calista Function
#'
#' @param DATA,INPUTS,Results
#' @keywords calista
#' @export
#' @examples
#' reorder_clusters_GUI()



reorder_clusters_GUI<-function(DATA,INPUTS,Results){
###########
###########if there's no time info, reorder cluster based on the connection
h=Results$TRANSITION$final_graph

  Results$TRANSITION$clusterTEXT=1
  hUG=h
  E(hUG)$weight=abs(E(hUG)$weight*2)
  dist=numeric()
  path=list()
  dist[1]=0

  for(i in 2:Results$expected_clusters){
    dist[i]=distances(hUG,v=1,to=i,algorithm = "dijkstra")
    if(length(unlist(shortest_paths(hUG,from = 1,to=i,output = "vpath",inbound.edges = FALSE)$vpath))>1){
      path[[i]]=unlist(shortest_paths(hUG,from = 1,to=i,output = "vpath",inbound.edges = FALSE)$vpath)
    }
  }
  cluster_progression=sort(dist)
  bb=order(dist)
  final_groups2=Results$final_groups
  c2=Results$c
  Order=numeric()
  Order[1]=1
  for(i in 2:Results$expected_clusters){
    Order[i]=path[[bb[i]]][length(path[[bb[i]]])]
    find_temp_idx=which(Results$final_groups==Order[i])
    final_groups2[find_temp_idx]=Results$cluster_predicted[i]
    c2[find_temp_idx,]=rep(Results$colorMARK_calista[i],length(find_temp_idx))
  }
  h=permute(h,Order)   ## check out the meaning in matlab code
  Results$final_groups=final_groups2
  Results$c=c2
  Results$cluster_progression=(cluster_progression-min(cluster_progression))/(max(cluster_progression)-min(cluster_progression))
  Results$cell_cluster_progression=matrix(0,length(Results$final_groups),1)
  for(i in 1:Results$expected_clusters){
    Results$cell_cluster_progression[which(Results$final_groups==Results$cluster_predicted[i])]=Results$cluster_progression[i]
  }
  # Results$cell_cluster_progression=(Results$cell_cluster_progression-min(Results$cell_cluster_progression))/max(Results$cell_cluster_progression-min(Results$cell_cluster_progression))
  x_center=matrix(0,Results$expected_clusters,1)
  y_center=x_center
  z_center=x_center
  for (i in 1:Results$expected_clusters) {
    x_center[i]=mean(Results$cell_cluster_progression[Results$final_groups==i])
    y_center[i]=mean(Results$score3[Results$final_groups==i,1])
    z_center[i]=mean(Results$score3[Results$final_groups==i,2])
  }
#####################################################################################################################
  ## Plot cell clustering results
  score3=Results$score3
  colorMARK_calista=Results$colorMARK_calista
  ClusterGroup_calista=Results$final_groups
  expected_clusters=Results$expected_clusters
  score3= Results$score3
  cluster_predicted=Results$cluster_predicted
  xMinLab=min(score3[,1])
  xMaxLab=max(score3[,1])
  yMinLab=min(score3[,2])
  yMaxLab=max(score3[,2])
  zMinLab=min(score3[,3])
  zMaxLab=max(score3[,3])

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
           scene = list(camera = list(eye = list(x = 2, y = 2, z = 0.5)),
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


  Results$c=c
  Results$plotly_cluster=plotly_cluster

  #####################################################################################################################
  ## Plot transition results
  normed_pseudotime=Results$cell_cluster_progression
  ClusterGroup2=ClusterGroup_calista
  colorMARK2=Results$colorMARK_calista
  legendInfo= Results$legendInfo_calista_transition
  gg_data=data.frame()
  label_color=c()
  label_names=c()
  for(k in 1:expected_clusters){
    gg_data=rbind(gg_data,data.frame(x=as.vector(normed_pseudotime[ClusterGroup2==k]),
                                     y=score3[ClusterGroup2==k,1],
                                     z=score3[ClusterGroup2==k,2],
                                     stage_colors=rep(colorMARK2[[k]],length(score3[ClusterGroup2==k,1])),
                                     gg_names=rep(paste('cluster',k),length(score3[ClusterGroup2==k,1]))
                                     #stages=rep(legendInfo[[k]],length(score3[ClusterGroup2==k,1]))

    ))
    label_color=c(label_color,colorMARK2[[k]])
    label_names=c(label_names,legendInfo[[k]])
  }
  h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
  colnames(h_after_layout)=c('x','y','z','stage_colors')
  txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
  # gg_data=rbind(gg_data,h_after_layout)

  #legend('topright',legend=legendInfo,pch = 16,col = colorMARK2)

  AddGraph_list=Results$TRANSITION$AddGraph_list
  AddGraph_list$h=h
  #AddGraph_list$i=NumberOfEdges
  #AddGraph_list$p=p
  AddGraph_list$gg_data=gg_data
  #AddGraph_list$segment_data=segment_data
  AddGraph_list$label_color=label_color
  AddGraph_list$label_names=label_names


  h_edge=get.edgelist(h)
  X=matrix(0,2,length(h_after_layout[h_edge[,1],1]))
  Y=X
  Z=X
  X[1,]=h_after_layout[h_edge[,1],1]
  X[2,]=h_after_layout[h_edge[,2],1]
  Y[1,]=h_after_layout[h_edge[,1],2]
  Y[2,]=h_after_layout[h_edge[,2],2]
  Z[1,]=h_after_layout[h_edge[,1],3]
  Z[2,]=h_after_layout[h_edge[,2],3]


  segment_data=data.frame(x=c(X),
                          y=c(Y),
                          z=c(Z))

  splitting=c(t(replicate(2, 1:(length(segment_data$x)/2))))
  ###########
  p <- plot_ly(gg_data) %>%
    add_trace(x = ~x, y = ~y, z = ~z, color = ~gg_names, colors = colorMARK2, type = 'scatter3d',
              mode = 'markers',size=1, alpha=0.8) %>%
    add_trace(x = segment_data$x, y = segment_data$y,
              z=segment_data$z, split=splitting,name='Graph',type = 'scatter3d',
              mode = 'lines', opacity = 1, line = list(color = 'rgb(22, 96, 167)',width = 5,alpha=0.8),
              showlegend = FALSE)  %>%
    layout(title = 'Lineage Inference',
           scene = list(camera = list(eye = list(x = -1.7, y = -1.7, z = 0.5)),
                        xaxis = list(title = 'Cluster progression',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     #range = c(xMinLim,xMaxLim),
                                     #type = 'log',
                                     #zerolinewidth = 1,
                                     #ticklen = 1,
                                     nticks=10,
                                     gridwidth = 3),

                        yaxis = list(title = 'COMP1',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     #range = c(yMinLim,yMaxLim),
                                     #zerolinewidth = 1,
                                     #ticklen = 1,
                                     nticks=10,
                                     gridwith = 3),
                        zaxis = list(title = 'COMP2',
                                     gridcolor = 'rgb(255, 255, 255)',
                                     #type = 'log',
                                     #range = c(zMinLim,zMaxLim),
                                     #zerolinewidth = 1,
                                     #ticklen = 1,
                                     nticks=10,
                                     gridwith = 3)),
           paper_bgcolor = 'rgb(243, 243, 243)',
           plot_bgcolor = 'rgb(243, 243, 243)')


  Results$TRANSITION$x_center=x_center
  Results$TRANSITION$y_center=y_center
  Results$TRANSITION$z_center=z_center
  Results$TRANSITION$final_graph=h
  Results$TRANSITION$p=p
  Results$TRANSITION$h_edge=h_edge
  Results$TRANSITION$AddGraph_list=AddGraph_list

  return(Results)
}
