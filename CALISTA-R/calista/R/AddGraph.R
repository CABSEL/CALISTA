#' A calista Function
#'
#' @param nodes,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,az,el,cluster_progression,DATA,...
#' @keywords calista
#' @export
#' @examples
#' AddGraph()


AddGraph<-function(nodes,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,
                   normed_pseudotime,score3,legendInfo,az,el, cluster_progression,DATA,...){
  ###test
  #AddGraph(nodes,Results$colorMARK_calista,x_center,y_center,z_center,Results$final_groups,Results$expected_clusters,
  #         Results$cell_cluster_progression,Results$score3,Results$legendInfo_calista,az,el)
  # colorMARK2=Results$colorMARK_calista
  # ClusterGroup2=Results$final_groups
  # expected_clusters=Results$expected_clusters
  # normed_pseudotime=Results$cell_cluster_progression
  # score3=Results$score3
  # legendInfo=Results$legendInfo_calista
  #first form of calling
  #AddGraph(nodes,Results$colorMARK_caession,Results$score3,Results$legendInfo_calista,az,el)
  varargin=list(...)
  if (length(varargin)>0){   #check out lista,x_center,y_center,z_center,Results$final_groups,Results$expected_clusters,
    #         Results$cell_cluster_progrhere, when coming into further code
    h=varargin[[1]]
    NumberOfEdges=nrow(get.edgelist(h))
  }else{
    h=NULL
    non=length(unique(as.vector(nodes[,1:2])))
    h=graph.empty(n=non,directed = FALSE)
    h_begin_layout=matrix(c(1:non,rep(1,non)),ncol = 2)
    l=layout.reingold.tilford(h)
    l=h_begin_layout

    if(length(unique(DATA$timeline))==1 && unique(DATA$timeline)==0){
      NumberOfEdges=1
      h=add_edges(h,c(nodes[1,1],nodes[1,2]),weight=nodes[1,3])
    } else {
      idx_outliers=which(nodes[,3] %in% boxplot(nodes[,3],plot=FALSE)$out)
      nodes_temp=nodes
      nodes_added=NULL
      NumberOfEdges=0
      unique_cluster_progression=unique(cluster_progression)
      jj_tot=NULL
      normalized_timeline=(DATA$timeline-min(DATA$timeline))/(max(DATA$timeline)-min(DATA$timeline))
      bin_tot_time_edge=unique(normalized_timeline)
      occurrence_tot_time=matrix(0, expected_clusters, length(bin_tot_time_edge))
      for (i in 1:expected_clusters){
        idx_cells_in_k=which(ClusterGroup2 %in% i)
        labels_in_k=normalized_timeline[idx_cells_in_k]
        occurrence_counts = as.data.frame(table(labels_in_k))
        idx_occurrence_counts=which(bin_tot_time_edge %in% occurrence_counts$labels_in_k)
        occurrence_tot_time[i,idx_occurrence_counts]=occurrence_counts$Freq
      }
      for (j in 2:length(unique_cluster_progression)){
        idx_actual_clusters=which(cluster_progression %in% unique_cluster_progression[j])

        for (k in 1:length(idx_actual_clusters)){
          cells_in_cluster=which(ClusterGroup2 %in% idx_actual_clusters[k])
          actual_times_in_cluster=normalized_timeline[cells_in_cluster]
          bin_edge=unique(actual_times_in_cluster)
          #%                 occurrence=hist(actual_times_in_cluster,bin_edge);
          #%                 % find time info for cells that occurr at least >5% of the total cells in the cluster
          #%                 selected_actual_times_in_cluster=bin_edge(find((occurrence/length(cells_in_cluster))>.10));
          #% remove clusters at actual/future cluster pseudotime
          idx_excluded_clusters=which(cluster_progression>=unique_cluster_progression[j])
          #%                 % find the cluster that has more cells at the previous time
          #%                 % point
          index_closest_time_point = which.min(abs(quantile(actual_times_in_cluster,.05)-bin_edge))    # time point at 5 percentile
          exit_while=0
          temp_occurrence_tot_time=occurrence_tot_time[,max(c(1,which(bin_tot_time_edge %in% bin_edge[index_closest_time_point])-1))]
          temp_occurrence_tot_time[idx_excluded_clusters]=0
          while(!exit_while)  {
            idx_thr=which.max(temp_occurrence_tot_time)
            exit_while2=0
            counting=1
            while(!exit_while2){
              if (sum(nodes[counting,1:2]-c(idx_thr, idx_actual_clusters[k]))==0){
                exit_while2=1
                idx_edge_to_check=counting
              }else{
                counting=counting+1
              }
            }

            if ( length(which(idx_outliers %in% idx_edge_to_check))==0){
              exit_while=1
            }else{
              temp_occurrence_tot_time[idx_thr]=0 #find the next max
            }
          }
          thr=which(unique_cluster_progression %in% cluster_progression[idx_thr])
          idx_previous_clusters=which(cluster_progression<unique_cluster_progression[j] & cluster_progression>=unique_cluster_progression[max(1,thr)])
          target=which(nodes[,2] %in% idx_actual_clusters[k])
          sources=which(nodes[,1] %in% idx_previous_clusters)
          idx_intersected_nodes=intersect(sources,target)
          idx_edge_2_add=which.max(nodes[idx_intersected_nodes,3])
          jj=idx_intersected_nodes[idx_edge_2_add]
           h=add_edges(h,c(nodes[jj,1],nodes[jj,2]),weight=nodes[jj,3])
           jj_tot=c(jj_tot, jj)
          NumberOfEdges=NumberOfEdges+1
          nodes_added=rbind(nodes_added, nodes[jj,])
        }
      }

      nodes_temp=nodes_temp[-jj_tot,]
      nodes=rbind(nodes_added, nodes_temp)
    }
    #plot.igraph(h,layout=l)
  }
  xMinLim=min(normed_pseudotime)#-0.1
  xMaxLim=max(normed_pseudotime)#+0.1
  yMinLim=min(score3[,1])-1
  yMaxLim=max(score3[,1])+1
  zMinLim=min(score3[,2])-1
  zMaxLim=max(score3[,2])+1
  #x11()
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
  h_edge=get.edgelist(h)
  V(h)$size=1
  V(h)$color=colorMARK2
  E(h)$width=2.5

  # X=matrix(0,2,length(h_after_layout[h_edge[,1],1]))
  # Y=X
  # Z=X
  # X[1,]=h_after_layout[h_edge[,1],1]
  # X[2,]=h_after_layout[h_edge[,2],1]
  # Y[1,]=h_after_layout[h_edge[,1],2]
  # Y[2,]=h_after_layout[h_edge[,2],2]
  # Z[1,]=h_after_layout[h_edge[,1],3]
  # Z[2,]=h_after_layout[h_edge[,2],3]
  #
  #
  # segment_data=data.frame(x=c(X),
  #                         y=c(Y),
  #                         z=c(Z))
  #
  # splitting=c(t(replicate(2, 1:(length(segment_data$x)/2))))
  # ###########
  # p <- plot_ly(gg_data) %>%
  #   add_trace(x = ~x, y = ~y, z = ~z, color = ~gg_names, colors = colorMARK2, type = 'scatter3d',
  #             mode = 'markers',size=1, alpha=0.8) %>%
  #   add_trace(x = segment_data$x, y = segment_data$y,
  #             z=segment_data$z, split=splitting,name='Graph',type = 'scatter3d',
  #             mode = 'lines', opacity = 1, line = list(color = 'rgb(22, 96, 167)',width = 5,alpha=0.8),
  #             showlegend = FALSE)  %>%
  #   layout(title = 'Lineage Inference',
  #          scene = list(camera = list(eye = list(x = -1.7, y = -1.7, z = 0.5)),
  #                       xaxis = list(title = 'Cluster progression',
  #                                    gridcolor = 'rgb(255, 255, 255)',
  #                                    #range = c(xMinLim,xMaxLim),
  #                                    #type = 'log',
  #                                    #zerolinewidth = 1,
  #                                    #ticklen = 1,
  #                                    nticks=10,
  #                                    gridwidth = 3),
  #
  #                       yaxis = list(title = 'COMP1',
  #                                    gridcolor = 'rgb(255, 255, 255)',
  #                                    #range = c(yMinLim,yMaxLim),
  #                                    #zerolinewidth = 1,
  #                                    #ticklen = 1,
  #                                    nticks=10,
  #                                    gridwith = 3),
  #                       zaxis = list(title = 'COMP2',
  #                                    gridcolor = 'rgb(255, 255, 255)',
  #                                    #type = 'log',
  #                                    #range = c(zMinLim,zMaxLim),
  #                                    #zerolinewidth = 1,
  #                                    #ticklen = 1,
  #                                    nticks=10,
  #                                    gridwith = 3)))
  #          # paper_bgcolor = 'rgb(243, 243, 243)',
  #          # plot_bgcolor = 'rgb(243, 243, 243)')
  # p
  #
  #############################

  AddGraph_list=list()
  AddGraph_list$h=h
  AddGraph_list$i=NumberOfEdges
  #AddGraph_list$p=p
  AddGraph_list$gg_data=gg_data
  #AddGraph_list$segment_data=segment_data
  AddGraph_list$label_color=label_color
  AddGraph_list$label_names=label_names
  AddGraph_list$nodes=nodes
  return(AddGraph_list)

}

