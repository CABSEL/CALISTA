CALISTA_transition<-function(Results,DATA){
  writeLines('\nCALISTA_transtion is running...\n\n')
  if(nargs()<2){
    stop('Not enough input variables')
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
  
  x_center=matrix(0,Results$expected_clusters,1)
  y_center=x_center
  z_center=x_center
  for (i in 1:Results$expected_clusters) {
    x_center[i]=mean(Results$cell_cluster_progression[Results$final_groups==i])
    y_center[i]=mean(Results$score3[Results$final_groups==i,1])
    z_center[i]=mean(Results$score3[Results$final_groups==i,2])
  }
  
################################################
  x11(title = '101')
  figure101=dev.cur()
  az=-37.5000
  el=30
  
  AddGraph_list=AddGraph(nodes,Results$colorMARK_calista,x_center,y_center,z_center,Results$final_groups,Results$expected_clusters,
                           Results$cell_cluster_progression,Results$score3,Results$legendInfo_calista_transition,az,el)
  h=AddGraph_list$h
  i=AddGraph_list$i
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
      ########plot ############
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
        label_names=c(label_names,Results$legendInfo_calista_transition[[k]])
        
      }
      names(label_color)=label_names
      h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
      colnames(h_after_layout)=c('x','y','z','stage_colors')
      txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
      gg_data=rbind(gg_data,h_after_layout)
      p1<-ggplot(gg_data,mapping = aes(x=y,y=z))
      p1<-p1+geom_point(data = gg_data,color=gg_data[,"stage_colors"])+geom_text(aes(label=txt_lab),hjust=0,vjust=0,size=50)
      p1<-p1+labs(x='COMP1',y='COMP2')
      p1<-p1+geom_point(mapping =aes(x=y_center,y=z_center),data = h_after_layout,size=4,color=h_after_layout[,"stage_colors"])
      h_edge=get.edgelist(h)
      segment_data=data.frame(x=h_after_layout[h_edge[,1],2],
                              y=h_after_layout[h_edge[,1],3],
                              xend=h_after_layout[h_edge[,2],2],
                              yend=h_after_layout[h_edge[,2],3])
      p1<-p1+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                          size=1)
      p1<-p1+ggtitle('Lineage Progression')
      plot(p1)
    }
  
           stop_adding=FALSE
           add_edge=TRUE
           MaxNumberOfEdges=nrow(nodes)
           while(stop_adding==FALSE)
           {
             NumberOfConnectedNodes=sum(!is.na(dfs(h,1,unreachable = FALSE)$order))
             writeLines(paste(i, 'edge(s) have been added and the graph is connected. \nIf you want to add another edge press "1" (then Enter)  \nIf you want to remove an edge press "2" (then Enter) \nIf you want to continue with the next step press "0" (then Enter) \n '))
             t=readinteger()
             true_key=FALSE
             #0:continue
             while(true_key==FALSE){
               if(t==0){
                 if(NumberOfConnectedNodes!=length(V(h))){
                   res <- tkmessageBox(title = "!! Warning !!",
                                       message = "Each cluster has to be connected by AT LEAST one edge!", icon = "info", type = "ok")
                   break
                 }
                 dev.off(figure101)
                 stop_adding=TRUE
                 x11(title="102")
                 figure102=dev.cur()
                 AddGraph_list=AddGraph(nodes,Results$colorMARK_calista,x_center,y_center,z_center,Results$final_groups,Results$expected_clusters,
                                          Results$cell_cluster_progression,Results$score3,Results$legendInfo_calista_transition,az,el,h)
                 
                 h=AddGraph_list$h
                 i=AddGraph_list$i
                 p=AddGraph_list$p
                 #p<-p+ggtitle('Cluster ordering') Lineage inference?????
                 p<-p+ggtitle('Lineage Progression')
                 plot(p)
               }
               #1: +
               if(t==1){
                 add_edge=TRUE
               }else if(t==2){
                 add_edge=FALSE
               }
               if(t %in% c(1,2,0)){
                 true_key=TRUE
               }else{
                 res1<- tkmessageBox(title = "!! Warning !!",
                                     message = "Please press '1','2','0'!", icon = "info", type = "ok")
                 t=readinteger()
               }
             }
           
             if(add_edge==TRUE && stop_adding==FALSE && true_key==TRUE){
               dev.off(figure101)
               x11(title = "101")
               figure101=dev.cur()
               MaxEdgesAdded=CheckNumberOfEdges(i+1,MaxNumberOfEdges)
               if(MaxEdgesAdded==FALSE){
                 i=i+1
                 h=AddNewEdge(h,nodes[i,],Results$colorMARK_calista,
                                            x_center,y_center,z_center,Results$final_groups,
                                            Results$expected_clusters,Results$cell_cluster_progression,
                                            Results$score3,Results$legendInfo_calista_transition,az,el);
                 #AddNewEdge(h,nodes(i,:),Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista,az,el);
               }else{
                 h=AddNewEdge(h,nodes[i,],Results$colorMARK_calista,
                              x_center,y_center,z_center,Results$final_groups,
                              Results$expected_clusters,Results$cell_cluster_progression,
                              Results$score3,Results$legendInfo_calista_transition,az,el);
                 writeLines('You have already added the maximum number of edges \n')
               }
             }
             if(add_edge==FALSE && stop_adding==FALSE && true_key==TRUE){
               h_temp=temp_RmEdge(h,nodes[i,],Results$colorMARK_calista,x_center,y_center,z_center,
                               Results$final_groups,Results$expected_clusters,
                               Results$cell_cluster_progression,Results$score3,
                               Results$legendInfo_calista_transition,az,el);
               NumberOfConnectedNodes2=sum(!is.na(dfs(h_temp,1,unreachable = FALSE)$order))
               if(NumberOfConnectedNodes2==length(V(h))){
                 dev.off(figure101)
                 x11(title = "101")
                 figure101=dev.cur()
                 h=RmEdge(h,nodes[i,],Results$colorMARK_calista,x_center,y_center,z_center,
                          Results$final_groups,Results$expected_clusters,
                          Results$cell_cluster_progression,Results$score3,
                          Results$legendInfo_calista_transition,az,el);
                 i=i-1
               }else{
                 res2<- tkmessageBox(title = "!! Warning !!",
                                     message = "The edge can not be removed. The graph must be CONNECTED!", icon = "info", type = "ok")
               }
             }
             
             
           }

  ###########if there's no time info, reorder cluster based on the connection
  hh=h  
  if(length(unique(DATA$timeline))==1 && unique(DATA$timeline)==0){
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
    for (i in 1:Results$expected_clusters) {
      x_center[i]=mean(Results$cell_cluster_progression[Results$final_groups==i])
      y_center[i]=mean(Results$score3[Results$final_groups==i,1])
      z_center[i]=mean(Results$score3[Results$final_groups==i,2])
    }
    dev.off(figure102)
    x11(title = '102')
    figure102=dev.cur()
    gg_data=data.frame()
    label_color=c()
    label_names=c()
    stt1='Cluster:'
    stt2='Cluster pseudotime:'
    for (clust in 1:Results$expected_clusters) {
      Results$legendInfo_calista_transition[clust]=paste(stt1,clust,stt2,Results$cluster_progression[clust])
    }
    for(k in 1:Results$expected_clusters){
      gg_data=rbind(gg_data,data.frame(x=Results$cell_cluster_progression[Results$final_groups==k],
                                       y=Results$score3[Results$final_groups==k,1],
                                       z=Results$score3[Results$final_groups==k,2],
                                       #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[Results$final_groups==k,1])),
                                       stage_colors=rep(Results$colorMARK_calista[k],length(Results$score3[Results$final_groups==k,1]))))
      label_color=c(label_color,Results$colorMARK_calista[k])
      label_names=c(label_names,Results$legendInfo_calista_transition[[k]])
    }
    names(label_color)=label_names
    h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
    colnames(h_after_layout)=c('x','y','z','stage_colors')
    txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
    gg_data=rbind(gg_data,h_after_layout)
    p1<-ggplot(gg_data,mapping = aes(x=y,y=z))
    p1<-p1+geom_point(data = gg_data,color=gg_data[,"stage_colors"])+geom_text(aes(label=txt_lab),hjust=0,vjust=0,size=20)
    p1<-p1+labs(x='COMP1',y='COMP2')
    # p1<-p1+scale_color_discrete(name='cluster',
    #                           #breaks=label_names,
    #                           values = label_color, 
    #                           labels=label_names
    #                         )
    p1<-p1+geom_point(mapping =aes(x=y_center,y=z_center),data = h_after_layout,size=4,color=h_after_layout[,"stage_colors"])
    h_edge=get.edgelist(h)
    segment_data=data.frame(x=h_after_layout[h_edge[,1],2],
                            y=h_after_layout[h_edge[,1],3],
                            xend=h_after_layout[h_edge[,2],2],
                            yend=h_after_layout[h_edge[,2],3])
    p1<-p1+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                        size=1)
    p1<-p1+ggtitle('plot after cluster relabeling')
    plot(p1)
    
    dev.off(Results$figure1000)
    dev.off(Results$figure1001)
    x11(title="1000")
    Results$figure1000=dev.cur()
    xMinLab=min(Results$score3[,1])
    xMaxLab=max(Results$score3[,1])
    yMinLab=min(Results$score3[,2])
    yMaxLab=max(Results$score3[,2])
    zMinLab=min(Results$score3[,3])
    zMaxLab=max(Results$score3[,3])
    for (i in 1:Results$expected_clusters) {
      scatterplot3d(Results$score3[Results$final_groups==i,1],
                    Results$score3[Results$final_groups==i,2],
                    Results$score3[Results$final_groups==i,3],
                    xlim = c(xMinLab,xMaxLab),
                    ylim = c(yMinLab,yMaxLab),
                    zlim = c(zMinLab,zMaxLab),
                    xlab = "PC1",
                    ylab = "PC2",
                    zlab = "PC3",
                    color = Results$colorMARK_calista[i],
                    pch=16)
      par(new=TRUE)
    }
    title(main = 'Cell Clustering(3D) after relabelling')
    legend('topright',unlist(Results$legendInfo_calista), pch = 16,col = Results$colorMARK_calista[1:Results$expected_clusters])
    
  }
  ###remove unwanted edges before additional analysis
  hh=h
  writeLines('Press 1 if you want to remove edges, 0 otherwise:')
  proceed=as.integer(readLines(n=1))
  # if(1){
  #   proceed=readinteger()
  # }
  
  if(proceed){
    Results$TRANSITION$clusterTEXT=1
    #dev.off(figure102)
    #x11(title = '102')
    #figure102=dev.cur()
    nodes_connection=get.edgelist(h)
    add_paths=TRUE
    count=0
    remove_edge=list()
    while(add_paths){
      writeLines('\n ******************************************************** \n')
      Edge=read2integer()
      remove_edge=list(Edge)
      remove_edge=matrix(unlist(remove_edge),byrow = TRUE,ncol=2)
      hh_temp=delete.edges(hh,get.edge.ids(h,remove_edge))
      mininal_graph=sum(!is.na(dfs(hh_temp,1,unreachable = FALSE)$order))
      if(mininal_graph==length(V(h))){
        count=count+1
        hh=hh_temp
      }else{
        writeLines("WARNING:the choosen edge can not be removed.All nodes should be connected by at least one edge\n")
      }
      writeLines(' Press 1 to remove another edge, 0 otherwise:')
      continue_add=as.integer(readinteger())
      if(continue_add!=1){
        add_paths=FALSE
      }
    }
    writeLines(paste(count,'edges to remove'))
    #h<-delete.edges(h,get.edge.ids(h,c(node[1],node[2])))
    dev.off(figure102)
    x11(title = '')
    figure102=dev.cur()
    gg_data=data.frame()
    label_color=c()
    label_names=c()
    
    stt1='Cluster:'
    stt2='Cluster pseudotime:'
    for (clust in 1:Results$expected_clusters) {
      Results$legendInfo_calista_transition[clust]=paste(stt1,clust,stt2,Results$cluster_progression[clust])
    }
    for(k in 1:Results$expected_clusters){
      gg_data=rbind(gg_data,data.frame(x=Results$cell_cluster_progression[Results$final_groups==k],
                                       y=Results$score3[Results$final_groups==k,1],
                                       z=Results$score3[Results$final_groups==k,2],
                                       #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[Results$final_groups==k,1])),
                                       stage_colors=rep(Results$colorMARK_calista[k],length(Results$score3[Results$final_groups==k,1]))))
      label_color=c(label_color,Results$colorMARK_calista[k])
      label_names=c(label_names,Results$legendInfo_calista_transition[[k]])
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
    h_edge=get.edgelist(hh)
    segment_data=data.frame(x=h_after_layout[h_edge[,1],2],
                            y=h_after_layout[h_edge[,1],3],
                            xend=h_after_layout[h_edge[,2],2],
                            yend=h_after_layout[h_edge[,2],3])
    p1<-p1+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                        size=1)
    p1<-p1+ggtitle('Lineage Progression')
    plot(p1)
  }
  ##########
  
  
  
  
  
  Results$TRANSITION$x_center=x_center
  Results$TRANSITION$y_center=y_center
  Results$TRANSITION$z_center=z_center
  Results$TRANSITION$final_graph=hh
  if(exists("Order")){
    #print()
    Results$mean_prob_in_out_cluster=mean_prob_in_out_cluster[Order,]
  }
  
  Results$mean_prob_in_out_cluster=mean_prob_in_out_cluster
  
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
    #unique_cluster_progression=unique(Results$cluster_progression)
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
                    main=paste("CALISTA cluster pseudotime"),
                    pch=16,
                    color='grey')
      
      #legend('topright',paste("CALISTA cluster pseudotime",unique_cluster_progression[iii]), pch = 16)
      par(new=TRUE)
    }
    par(new=FALSE)
    
  }
  Results$stages_name=stages_name
  hh=Results$TRANSITION$final_graph 
  nodes_connection3=get.edgelist(hh)
  Results$TRANSITION$nodes_connection=nodes_connection3
  Results$cluster_distance=cluster_distance
  return(Results)
}