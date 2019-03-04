#' A calista Function
#'
#' @param DATA,Results
#' @keywords calista
#' @export
#' @examples
#' jump_transition()


jump_transition <- function(DATA,Results){

  ClusterGroup2=Results$final_groups
  expected_clusters=Results$expected_clusters
  normed_pseudotime=Results$cell_cluster_progression
  score3=Results$score3
  legendInfo=Results$legendInfo_calista
  colorMARK2=Results$colorMARK_calista

  x_center=matrix(0,Results$expected_clusters,1)
  y_center=x_center
  z_center=x_center
  for (i in 1:Results$expected_clusters) {
    x_center[i]=mean(Results$cell_cluster_progression[Results$final_groups==i])
    y_center[i]=mean(Results$score3[Results$final_groups==i,1])
    z_center[i]=mean(Results$score3[Results$final_groups==i,2])
  }


  x11(title = '101')
  figure101=dev.cur()
  az=-37.5000
  el=30

  xMinLim=min(normed_pseudotime)-1
  xMaxLim=max(normed_pseudotime)+1
  yMinLim=min(score3[,1])-1
  yMaxLim=max(score3[,1])+1
  zMinLim=min(score3[,2])-1
  zMaxLim=max(score3[,2])+1

  gg_data=data.frame()
  label_color=c()
  label_names=c()
  for(k in 1:Results$expected_clusters){
    gg_data=rbind(gg_data,data.frame(x=Results$cell_cluster_progression[Results$final_groups==k],
                                     y=Results$score3[Results$final_groups==k,1],
                                     z=Results$score3[Results$final_groups==k,2],
                                     #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[Results$final_groups==k,1])),
                                     stage_colors=rep(Results$colorMARK_calista[k],length(Results$score3[Results$final_groups==k,1]))))
    label_color=c(label_color,Results$colorMARK_calista[k])
    label_names=c(label_names,Results$legendInfo_calista[[k]])
  }
  names(label_color)=label_names
  h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
  colnames(h_after_layout)=c('x','y','z','stage_colors')
  txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
  gg_data=rbind(gg_data,h_after_layout)
  p1<-ggplot(gg_data,mapping = aes(x=y,y=z))
  p1<-p1+geom_point(data = gg_data,color=gg_data[,"stage_colors"])+geom_text(aes(label=txt_lab),hjust=0,vjust=0,size=20)
  p1<-p1+labs(x='PC1',y='PC2')
  # p1<-p1+scale_color_manual(name='cluster',
  #                           breaks=label_names,
  #                           values = label_color,
  #                           labels=label_names)
  p1<-p1+geom_point(mapping =aes(x=y_center,y=z_center),data = h_after_layout,size=4,color=h_after_layout[,"stage_colors"])
  p1<-p1+ggtitle('Lineage inference')
  plot(p1)

  # Select paths maunally

  add_paths=TRUE
  count=1
  CELL_path=list()
  tot_nodes=vector()
  while(add_paths){
    writeLines(paste("Path num: ",count))
    writeLines("\n*********************************\n")
    writeLines("Please check the figure 'Lineage Inference.R' for reference and  ")
    writeLines('Key the clusters in the path (e.g. 1 2 3 4):')
    temp_path=strsplit(readLines(n=1), " ")
    temp_path=c(as.integer(temp_path[[1]]))
    CELL_path[[count]]=temp_path
    path_len=length(CELL_path[[count]])
    nodes=cbind(CELL_path[[count]][1:path_len-1],CELL_path[[count]][2:path_len])
    tot_nodes=rbind(tot_nodes,nodes)
    writeLines('Please press 1 to add another path, 0 otherwise:')
    continue_add=as.integer(readLines(n=1))
    if(continue_add){
      count=count+1
    }else{
      add_paths=FALSE
    }
  }


  tot_nodes=unique(tot_nodes)# Construct the graph
  h=NULL
  non=length(unique(as.vector(tot_nodes[,1:2])))
  h=graph.empty(n=non,directed = FALSE)
  h_begin_layout=matrix(c(1:non,rep(1,non)),ncol = 2)
  l=layout.reingold.tilford(h)
  l=h_begin_layout
  NumberOfEdges=non-1
  h=add_edges(h,as.vector(t(tot_nodes[,1:2])),weight=1)

  # Check to have a connected graph

  NumberOfConnectedNodes=sum(!is.na(dfs(h,1,unreachable = FALSE)$order))
  if(NumberOfConnectedNodes!=length(V(h))){
    stop('The graph must be CONNECTED. Please run CALISTA again')
  }

  h_edge=get.edgelist(h)
  V(h)$size=1
  V(h)$color=colorMARK2
  E(h)$width=2.5
  p1<-p1+geom_point(mapping =aes(x=y_center,y=z_center),data = h_after_layout,size=4,color=h_after_layout[,"stage_colors"])
  segment_data=data.frame(x=h_after_layout[h_edge[,1],2],
                          y=h_after_layout[h_edge[,1],3],
                          xend=h_after_layout[h_edge[,2],2],
                          yend=h_after_layout[h_edge[,2],3])
  p1<-p1+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                    size=1)
  plot(p1)

  Results$PATH$CELL_path=CELL_path
  Results$TRANSITION$x_center=x_center
  Results$TRANSITION$y_center=y_center
  Results$TRANSITION$z_center=z_center
  Results$TRANSITION$final_graph=h

  #plot each predicted cluster separately
  p1=ceiling(sqrt(length(unique(Results$cluster_progression))))
  stages_name=unique(Results$cluster_progression)
  x11(title = 'plot each predicted cluster separately')
  par(mfrow=c(p1,p1))
  for(K in 1:length(stages_name)){
    actual_stage=which(Results$cluster_progression==stages_name[K])
    #print(actual_stage)
    # print('actural stage')
    # print(actual_stage)
    # print(length(actual_stage))
    xMinLim=min(Results$score3[,1])
    xMaxLim=max(Results$score3[,1])
    yMinLim=min(Results$score3[,2])
    yMaxLim=max(Results$score3[,2])
    zMinLim=min(Results$score3[,3])
    zMaxLim=max(Results$score3[,3])
    for(ii in 1:length(actual_stage)){  #  remeber in one veteor
      #print(actual_stage[ii])
      #print(Results$colorMARK_calista[actual_stage[ii]])
      temp_plot=which(Results$final_groups==actual_stage[ii])
      scatterplot3d(x=Results$score3[temp_plot,1],
                    y=Results$score3[temp_plot,2],
                    z=Results$score3[temp_plot,3],
                    xlim = c(xMinLim,xMaxLim),
                    ylim = c(yMinLim,yMaxLim),
                    zlim = c(zMinLim,zMaxLim),
                    xlab = "",
                    ylab = "",
                    zlab = "",
                    pch=16,
                    color=Results$colorMARK_calista[actual_stage[ii]])
      ## ?????????֣??????ڶ????????Σ?ÿ??Ӧ????��????ɫ??????????
      par(new=TRUE)
    }
    idx_remove=which(Results$cluster_predicted %in% actual_stage)

    cluster_plot=Results$cluster_predicted
    cluster_plot=cluster_plot[-idx_remove]
    for(iii in 1:length(cluster_plot)){
      temp_plot=which(Results$final_groups==cluster_plot[iii])
      scatterplot3d(Results$score3[temp_plot,1],
                    Results$score3[temp_plot,2],
                    Results$score3[temp_plot,3],
                    xlim = c(xMinLim,xMaxLim),
                    ylim = c(yMinLim,yMaxLim),
                    zlim = c(zMinLim,zMaxLim),
                    xlab = "PC1",
                    ylab = "PC2",
                    zlab = "PC3",
                    main="CALISTA cluster pseudotime",
                    pch=16,
                    color='grey')
      par(new=TRUE)
    }
    par(new=FALSE)

  }
  Results$stages_name=stages_name
  hh=Results$TRANSITION$final_graph
  nodes_connection3=get.edgelist(hh)
  Results$TRANSITION$nodes_connection=nodes_connection3


  return(Results)
}
