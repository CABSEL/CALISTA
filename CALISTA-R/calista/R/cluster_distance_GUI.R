#'  cluster_distance_GUI infer lineage progression among cell clusters
#'
#'  CALISTA uses the cluster distances - a measure of cluster-cluster
#'  dissimilarity - to infer the lineage progression or cluster-cluster
#'  relationship. The lineage progression graph is built based on adding
#'  edges between clusters in increasing magnitude of cluster distance.
#'  CALISTA also provides a simple user-interface to add and remove edges
#'  based on the cluster distances.
#'
#'  Usage:
#'  Run CALISTA lineage inference using CALISTA single-clustering result
#'  Results=CALISTA_transition_main(DATA,Results);
#'
#'  Run CALISTA lineage Inference with user-defined cell clusters
#' A calista Function
#'
#' @param DATA,INPUTS,Results
#' @keywords calista
#' @export
#' @examples
#' cluster_distance_GUI()



cluster_distance_GUI<-function(DATA,INPUTS,Results){
  writeLines('\nCALISTA_transtion is running...\n\n')
  if(nargs()<2){
    stop('Not enough input variables')
  }

  ## CHECK INPUT ARGUMENTS
  if (!exists('data_type', where=INPUTS)){
    stop('Please specify the data type in INPUTS.data_type')
  }

  if  (!exists('format_data', where=INPUTS)){
    stop('Please specify the format of the data in INPUTS$format_data')
  }

  # if  (!exists('cluster_time', where=INPUTS)){
  # INPUTS$cluster_time=0
  # }

  if  (!exists('cells_2_cut', where=INPUTS)){
    INPUTS$cells_2_cut=0;
  }

  if  (!exists('cut_variable_genes', where=INPUTS)){
    INPUTS$cut_variable_genes=1000
  }

  if(!exists('data_selection', where=INPUTS)){
    INPUTS$data_selection=integer()
  }

  #   if (INPUTS$cluster_time~=0 && ~INPUTS$data_selection==0){
  #   writelines('\nCALISTA Time CLustering is active. INPUTS$data_selection set as 0. All time points are considered for the analysis\n')
  #   INPUTS$data_selection=0
  # }

  if  (!exists('perczeros_genes', where=INPUTS)){
    INPUTS$perczeros_genes=100
  }

  if  (!exists('perczeros_cells', where=INPUTS)){
    INPUTS$perczeros_cells=100

  }

  if  (!exists('cells_2_cut', where=INPUTS)){
    INPUTS$cells_2_cut=0

  }

  if  (!exists('optimize', where=INPUTS)){
    INPUTS$optimize=0

  }

  if  (!exists('parallel', where=INPUTS)){
    INPUTS$parallel=1

  }

  if  (!exists('runs', where=INPUTS)){
    INPUTS$runs=50

  }

  if  (!exists('max_iter', where=INPUTS)){
    INPUTS$max_iter=100

  }


  if  (!exists('Cluster', where=INPUTS)){
    INPUTS$Cluster='kmedoids'

  }

  if  (!exists('thr_transition_genes', where=INPUTS)){
    INPUTS$thr_transition_genes=50

  }



  if(length(unique(DATA$timeline))==1 && unique(DATA$timeline)==0){
    Results$legendInfo_calista_transition=Results$legendInfo_calista
  }
  else{
    stt1='Cluster:'
    stt2='Cluster pseudotime:'
    for (clust in 1:Results$expected_clusters) {
      #as.numeric(format(round(0.905, 2), nsmall = 2))
      Results$legendInfo_calista_transition[clust]=paste(stt1,clust,stt2,format(round(Results$cluster_progression[clust], 2), nsmall = 2))
    }
  }



  Results$TRANSITION$clusterTEXT=0 # for calista_path
  my_results_final=Results$clustering_struct
  my_results_final$all$all$distance[my_results_final$all$all$distance==0]=Inf
  neighbour=as.matrix(apply(my_results_final$all$all$distance,
                            1, which.min))
  all_in_one=cbind(neighbour,t(Results$final_groups))
  all_in_one_sorted=all_in_one[order(all_in_one[,2]),]
  aa=matrix(0,Results$expected_clusters,Results$expected_clusters)
  for (i in 1:Results$expected_clusters) {
    for (j in 1:Results$expected_clusters) {
      aa[i,j]=sum(all_in_one_sorted[all_in_one_sorted[,2]==i,1]==j)
    }
  }
  my_mean=matrix(0,Results$expected_clusters,Results$expected_clusters)
  xxx=my_results_final$all$all$cell_prob/matrix(rep(rowSums(my_results_final$all$all$cell_prob),Results$expected_clusters),ncol = Results$expected_clusters)
  mean_prob_in_out_cluster=matrix(0,Results$expected_clusters,2)
  for (i in 1:Results$expected_clusters) {
    aaa=which(Results$final_groups==i)
    my_mean[i,]=colMeans(my_results_final$all$all$cell_prob[aaa,])
    mean_prob_in_out_cluster[i,1]=my_mean[i,i]
    temp_my_mean=my_mean[i,]
    temp_my_meanp=temp_my_mean[-i]
    mean_prob_in_out_cluster[i,2]=mean(temp_my_mean)
    # ask why hold on here matlab
  }

  ############################################
  distance=(my_mean-diag(my_mean)%*%matrix(1,1,Results$expected_clusters))/DATA$numGENES
  distance_vector=as.vector(distance)
  sorted_dist=sort(distance_vector,decreasing = TRUE)
  idx_dist=order(distance_vector,decreasing = TRUE)
  I=matrix(0,Results$expected_clusters*Results$expected_clusters,3)
  I[,2]=((idx_dist-1)%%Results$expected_clusters)+1
  I[,1]=floor((idx_dist-1)/Results$expected_clusters)+1
  I[,3]=sorted_dist
  nodes_n=I
  nodes_n=nodes_n[-(1:Results$expected_clusters),]
  nodes_all=nodes_n
  nodes_all[,1]=apply(nodes_n[,1:2], 1, min)
  nodes_all[,2]=apply(nodes_n[,1:2], 1, max)
  unique_nodes_all=unique(nodes_all[,1:2])
  unique_nodes_all=unique_nodes_all[order(unique_nodes_all[,1]),]
  ind_dir=apply(unique_nodes_all, 1, function(x){
    which(apply(nodes_all[,1:2], 1, function(z) all(x==z)))[1]
  })
  nodes=nodes_all[ind_dir,]
  nodes=nodes[order(nodes[,3],decreasing = TRUE),]
  cluster_distance=nodes;
  cluster_distance[,3]=abs(cluster_distance[,3])
  cluster_distance_tot=cluster_distance

  x_center=matrix(0,Results$expected_clusters,1)
  y_center=x_center
  z_center=x_center
  for (i in 1:Results$expected_clusters) {
    x_center[i]=mean(Results$cell_cluster_progression[Results$final_groups==i])
    y_center[i]=mean(Results$score3[Results$final_groups==i,1])
    z_center[i]=mean(Results$score3[Results$final_groups==i,2])
  }

  az=-37.5000
  el=30
  score3=Results$score3
  colorMARK2=Results$colorMARK_calista
  ClusterGroup2=Results$final_groups
  expected_clusters=Results$expected_clusters
  normed_pseudotime=Results$cell_cluster_progression
  legendInfo= Results$legendInfo_calista_transition
  AddGraph_list=AddGraph(nodes,colorMARK2,x_center,y_center,z_center,Results$final_groups,Results$expected_clusters,
                         Results$cell_cluster_progression,Results$score3,Results$legendInfo_calista_transition,az,el,Results$cluster_progression,DATA)

  h=AddGraph_list$h
  i=AddGraph_list$i
  nodes=AddGraph_list$nodes
  gg_data=AddGraph_list$gg_data
  label_color=AddGraph_list$label_color
  label_names=AddGraph_list$label_names


  MaxNumberOfEdges=nrow(nodes)
  connected=FALSE
  while(connected==FALSE){
    NumberOfConnectedNodes=sum(!is.na(dfs(h,1,unreachable = FALSE)$order))
    if(NumberOfConnectedNodes==length(V(h))){
      connected=TRUE
    }else{
      MaxEdgesAdded=CheckNumberOfEdges(i+1,MaxNumberOfEdges)
      if(MaxEdgesAdded==FALSE){
        i=i+1
        h=add_edges(h,c(nodes[i,1],nodes[i,2]),weight=nodes[i,3])
      }
    }


    h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
    colnames(h_after_layout)=c('x','y','z','stage_colors')
    txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
    #gg_data=rbind(gg_data,h_after_layout)
  }


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

  p
  Results$TRANSITION$x_center=x_center
  Results$TRANSITION$y_center=y_center
  Results$TRANSITION$z_center=z_center
  Results$TRANSITION$final_graph=h
  Results$TRANSITION$p=p
  Results$TRANSITION$h_edge=h_edge
  Results$TRANSITION$AddGraph_list=AddGraph_list
  Results$TRANSITION$cluster_distance
  Results$TRANSITION$nodes=nodes

  ## Select top cluster distances to show

  idx_edges_of_h=match(data.frame(t(h_edge)), data.frame(t(cluster_distance_tot[,1:2])))
  #idx_edges_of_h[is.na(idx_edges_of_h)] <- Inf
  idx_edges_of_h=sort(idx_edges_of_h)
  edges_cutoff_2_display=min(200,nrow(cluster_distance_tot))
  cluster_distance_temp=cluster_distance_tot[-idx_edges_of_h,]
  cluster_distance_selected=rbind(cluster_distance_tot[idx_edges_of_h,],cluster_distance_temp[1:(edges_cutoff_2_display-length(idx_edges_of_h)),])
  ww=paste('**GUI displaying top',nrow(cluster_distance_selected),'cluster distances **')
  writeLines(ww)
  cluster_distance_selected[,3]=round(cluster_distance_selected[,3],2)

  Results$TRANSITION$check_boxes=paste(cluster_distance_selected[,1],'-',cluster_distance_selected[,2],'  (',cluster_distance_selected[,3],') ')
  Results$TRANSITION$cluster_distance_selected=cluster_distance_selected

  return(Results)


}
