CALISTA_ordering2<-function(Results,Parameters,DATA,n_points_interp){
  writeLines('\nCALISTA_ordering is running...\n')
  if(nargs()<4){
    stop('Not enough input variables.')
  }
  h=Results$TRANSITION$final_graph
  cell_prob_3D=get_3D_cell_prob(Results,Parameters,DATA)
  nodes_connection3=get.edgelist(h)
  n_edges=nrow(nodes_connection3)
  cell_ordering=matrix(0,DATA$nvars,2)
  for(CELL in 1:DATA$nvars){
    actual_cluster=Results$final_groups[CELL]
    edges_2_check=which(nodes_connection3[,1]==actual_cluster | nodes_connection3[,2]==actual_cluster)
    max_log_P_each_edge=numeric()
    time_max_log_P_each_edge=numeric()
    for(i in 1:length(edges_2_check)){
      clust1=which(Results$cluster_predicted==nodes_connection3[edges_2_check[i],1])
      clust2=which(Results$cluster_predicted==nodes_connection3[edges_2_check[i],2])
      t1=Results$cluster_progression[clust1]
      t2=Results$cluster_progression[clust2]
      same_t=0
      if(t1==t2){
        same_t=1
        tequi=t1
        t1=0
        t2=1
        t_interp=seq(t1,t2,length=n_points_interp)
      }else{
        t_interp=seq(t1,t2,length=n_points_interp)
      }
      sum_log_P_genes=matrix(0,1,length(t_interp))
      for(j in 1:DATA$numGENES){
        int1=cell_prob_3D[CELL,clust1,j]
        int2=cell_prob_3D[CELL,clust2,j]
        Xnew=approx(c(t1,t2),c(int1,int2),xout = t_interp)
        sum_log_P_genes=sum_log_P_genes+log2(Xnew[[2]])
      }
      max_log_P_each_edge[i]=max(sum_log_P_genes)
      idx_max_log_P_each_edge=which(sum_log_P_genes==max_log_P_each_edge[i])
      if (same_t){
        time_max_log_P_each_edge[i]=tequi
      }else{
        time_max_log_P_each_edge[i]=t_interp[idx_max_log_P_each_edge]
      }
    }
    
    idx_max_log_P_cell=which(max_log_P_each_edge==max(max_log_P_each_edge))[1]
    if(same_t){
      cell_ordering[CELL,1]=tequi
    }else{
      cell_ordering[CELL,1]=time_max_log_P_each_edge[idx_max_log_P_cell]
    }
    cell_ordering[CELL,2]=edges_2_check[idx_max_log_P_cell]
  }
  normed_cell_ordering=(cell_ordering[,1]-min(cell_ordering[,1]))/max(cell_ordering[,1]-min(cell_ordering[,1]))
  
  
  x_center=matrix(0,Results$expected_clusters,1)
  y_center=Results$TRANSITION$y_center
  z_center=Results$TRANSITION$z_center
  for( i in 1:Results$expected_clusters){
    x_center[i]=mean(normed_cell_ordering[Results$final_groups==i])
  }
  
  
  x11(title = 'CALISTA_ordering2.R')
  #######################dong dong dong, plot again 
  gg_data=data.frame()
  label_color=c()
  label_names=c()
  for(k in 1:Results$expected_clusters){
    gg_data=rbind(gg_data,data.frame(x=normed_cell_ordering[Results$final_groups==Results$cluster_predicted[k]],
                                     y=Results$score3[Results$final_groups==Results$cluster_predicted[k],1],
                                     z=Results$score3[Results$final_groups==Results$cluster_predicted[k],2],
                                     gg_names=rep(paste('cluster',k),length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1])),
                                     gg_color=rep(Results$colorMARK_calista[k],length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1]))))
    label_color=c(label_color,Results$colorMARK_calista[k])
    label_names=c(label_names,paste('cluster',k))
  }
  names(label_color)=label_names
  p1<-ggplot(gg_data)
  p1<-p1+geom_point(mapping = aes(x=y,y=z,color=gg_names),data = gg_data)
  p1<-p1+labs(x='PC1',y='PC2')
  p1<-p1+scale_color_manual(name='cluster',
                          breaks=label_names,
                          values = label_color,
                          labels=label_names)
  h_after_layout=data.frame(x_center,y_center,z_center,label_names,stringsAsFactors = FALSE)
  colnames(h_after_layout)=c('x_center','y_center','z_center','label_name')
  p1<-p1+geom_point(mapping =aes(x=y_center,y=z_center,color=label_name),data = h_after_layout,size=10)
  h_edge=get.edgelist(h)
  segment_data=data.frame(x=h_after_layout[h_edge[,1],2],
                          y=h_after_layout[h_edge[,1],3],
                          xend=h_after_layout[h_edge[,2],2],
                          yend=h_after_layout[h_edge[,2],3])
  p1<-p1+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                    size=1)
  ################################
  p2<-ggplot(gg_data)
  p2<-p2+geom_point(mapping = aes(x=x,y=y,color=gg_names),data = gg_data)
  p2<-p2+labs(x='cell ordering',y='PC1')
  p2<-p2+scale_color_manual(name='cluster',
                            breaks=label_names,
                            values = label_color,
                            labels=label_names)
  p2<-p2+geom_point(mapping =aes(x=x_center,y=y_center,color=label_name),data = h_after_layout,size=10)
  h_edge=get.edgelist(h)
  segment_data=data.frame(x=h_after_layout[h_edge[,1],1],
                          y=h_after_layout[h_edge[,1],2],
                          xend=h_after_layout[h_edge[,2],1],
                          yend=h_after_layout[h_edge[,2],2])
  p2<-p2+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                      size=1)
  ##############################
  #######################
  p3<-ggplot(gg_data)
  p3<-p3+geom_point(mapping = aes(x=x,y=z,color=gg_names),data = gg_data)
  p3<-p3+labs(x='cell ordering',y='PC2')
  p3<-p3+scale_color_manual(name='cluster',
                            breaks=label_names,
                            values = label_color,
                            labels=label_names)
  p3<-p3+geom_point(mapping =aes(x=x_center,y=z_center,color=label_name),data = h_after_layout,size=10)
  h_edge=get.edgelist(h)
  segment_data=data.frame(x=h_after_layout[h_edge[,1],1],
                          y=h_after_layout[h_edge[,1],3],
                          xend=h_after_layout[h_edge[,2],1],
                          yend=h_after_layout[h_edge[,2],3])
  p3<-p3+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                      size=1)
  ##########################
  grid.arrange(p1,p2,p3)
  Results$ORDERING$cell_ordering=cell_ordering
  Results$ORDERING$normed_cell_ordering=normed_cell_ordering
  Results$TRANSITION$nodes_connection=nodes_connection3
  idx_actual_edge=list()
  cells_assigned_to_edge=numeric()
  for(i in 1:n_edges){
    idx_actual_edge[[i]]=which(Results$ORDERING$cell_ordering[,2]==i)
    bb=order(Results$ORDERING$cell_ordering[idx_actual_edge[[i]],1])
    idx_actual_edge[[i]]=idx_actual_edge[[i]][bb]
    cells_assigned_to_edge[i]=length(idx_actual_edge[[i]])
  }
  Results$ORDERING$cells_assigned_to_edge=cells_assigned_to_edge
  Results$ORDERING$idx_actural_edge=idx_actual_edge
  Results$ORDERING$idx_sorted_cells=order(Results$ORDERING$cell_ordering[,1])
  return(Results)
  
  
}