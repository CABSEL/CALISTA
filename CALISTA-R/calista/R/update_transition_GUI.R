#' A calista Function
#'
#' @param DATA,INPUTS,Results,nodes,edges_chosen
#' @keywords calista
#' @export
#' @examples
#' update_transition_GUI()



update_transition_GUI<-function(DATA,INPUTS,Results,nodes,edges_chosen){

  h=NULL
  non=length(unique(as.vector(nodes[,1:2])))
  h=graph.empty(n=non,directed = FALSE)
  h_begin_layout=matrix(c(1:non,rep(1,non)),ncol = 2)
  l=layout.reingold.tilford(h)
  l=h_begin_layout
  for(i in 1:length(edges_chosen)){
    h=add_edges(h,c(nodes[i,1],nodes[i,2]),weight=nodes[i,3])
  }
  score3=Results$score3
  colorMARK2=Results$colorMARK_calista
  ClusterGroup2=Results$final_groups
  expected_clusters=Results$expected_clusters
  normed_pseudotime=Results$cell_cluster_progression
  legendInfo= Results$legendInfo_calista_transition
  gg_data=Results$TRANSITION$AddGraph_list$gg_data
  label_color=Results$TRANSITION$AddGraph_list$label_color
  label_names=Results$TRANSITION$AddGraph_list$label_names
  x_center=Results$TRANSITION$x_center
  y_center=Results$TRANSITION$y_center
  z_center=Results$TRANSITION$z_center
  h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
  colnames(h_after_layout)=c('x','y','z','stage_colors')
  txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))

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
           scene = list(camera = list(eye = list(x = -0.1, y = -2.5, z = 0.9)),
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

  #p
  Results$TRANSITION$x_center <- x_center
  Results$TRANSITION$y_center <- y_center
  Results$TRANSITION$z_center <- z_center
  Results$TRANSITION$final_graph <- h
  Results$TRANSITION$p <- p
  return(Results)

}
