#   %CALISTA_ORDERING_MAIN perform pseudotemporal ordering of cells
#   % For pseudotemporal ordering of cells, CALISTA performs maximum
#   % likelihood optimization for each cell using a linear interpolation of the
#   % cell likelihoods between any two connected clusters. 
#   % 
#   % Usage:
#   % 
#   % 1- Run CALISTA cell ordering using results of CALISTA clustering and lineage inference
#   %    Results=CALISTA_ordering_main(DATA,Results);
#   %  
#   % 2- Run CALISTA cell ordering using user-defined cell assignments
#   %    Results=list()
#   %    Results=CALISTA_ordering_main(DATA,Results,cell_assignments)
#   %
#   % CALISTA will ask users to specify sequences of connected clusters, i.e. 
#   % paths, in the lineage graph. The list of edges in the graph is the union 
#   % of all edges in the user-specified paths.
#   %    
#   % Inputs:   
#   % DATA - a structure containing preprocessed single-cell expression data
#   % Use 'import_data' to upload and preprocess single-cell expression values.
#   % 
#   % Results - a structure of CALISTA clustering and lineage inference results
#   % Run 'CALISTA_clustering_main' and/or 'CALISTA_transition_main'
#   %
#   % cell_assignments - 1xN vector of INTEGERS with N = number of cells. 
#   % The n-th element of cell_assignments contains the cluster assignment of
#   % the n-th cell of the expression data uploaded. Cluster names must
#   % be assigned in sequence (e.g. 1,2,3,4 and not 1,2,4).
#   %
#   % Outputs:
#   % Results - a structure containing the results of CALISTA analysis. 
#   % The most relevant fields containing the cell ordering are:
#   %
#   % Results$ORDERING$normed_cell_ordering - a vector contaning the pseudotime
#   % of the cells.
#   % 
#   % Results$ORDERING$idx_sorted_cells - a vector containing the ordered cell
#   % indices in increasing order of the cell pseudotime.
#   % 
#   % Results$ORDERING$idx_actual_edge: a 1xE cell array with the i-th cell 
#   % containing the vector of cell indices assigned in the i-th edge. The
#   % edges are defined in Results.TRANSITION.nodes_connection.
#   % 
#   % Created by Nan Papili Gao (R version implemented by Tao Fang)
#   %            Institute for Chemical and Bioengineering 
#   %            ETH Zurich
#   %            E-mail:  nanp@ethz.ch
#   %
#   % Copyright. June 1, 2017.

CALISTA_ordering_main<-function(DATA,Results,cell_assignments){
  
   if(nargs()<2){
    stop('Not enough input variables.')
  }
  
  if (length(Results)==0){
    if (nargs()<3){
      stop("Not enough input arguments")
    }
    Results=jump_clustering(DATA,cell_assignments)
    Results=jump_transition(DATA,Results)
    
    ###3.b- plot mean expression for each cluster
    Results=Plot_cluster_mean_exp(Results,DATA)
    
    ###3.c- cell-celll variability analysis
    Results=cell_variability(Results,DATA)
  }
  
  
  writeLines('\nCALISTA_ordering is running...\n')
  Parameters=DATA$Parameters
  n_points_interp=100
  hh=Results$TRANSITION$final_graph 
  cell_prob_3D=get_3D_cell_prob(Results,Parameters,DATA)
  nodes_connection3=get.edgelist(hh)
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
  
  
  x11(title = 'Cell Ordering.R')
  #######################dong dong dong, plot again 
  gg_data=data.frame()
  label_color=c()
  label_names=c()
  for(k in 1:Results$expected_clusters){
    gg_data=rbind(gg_data,data.frame(x=normed_cell_ordering[Results$final_groups==Results$cluster_predicted[k]],
                                     y=Results$score3[Results$final_groups==Results$cluster_predicted[k],1],
                                     z=Results$score3[Results$final_groups==Results$cluster_predicted[k],2],
                                     stage_colors=rep(Results$colorMARK_calista[k],length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1]))
                                     #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1])),
                                     #gg_names=rep(paste('cluster',k),length(Results$score3[Results$final_groups==Results$cluster_predicted[k],1]))
                                     ))
    label_color=c(label_color,Results$colorMARK_calista[k])
    label_names=c(label_names,Results$legendInfo_calista[[k]])
  }
  names(label_color)=label_names
  #p1<-ggplot(gg_data)
  #p1<-p1+geom_point(mapping = aes(x=y,y=z),data = gg_data,color=gg_data[,"stage_colors"])
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
  ################################
  #p2<-ggplot(gg_data)
  #p2<-p2+geom_point(mapping = aes(x=x,y=y),data = gg_data,color=gg_data[,"stage_colors"])
  #h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
  #colnames(h_after_layout)=c('x','y','z','stage_colors')
  #txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
  #gg_data=rbind(gg_data,h_after_layout)
  p2<-ggplot(gg_data,mapping = aes(x=x,y=y))
  p2<-p2+geom_point(data = gg_data,color=gg_data[,"stage_colors"])+geom_text(aes(label=txt_lab),hjust=0,vjust=0,size=20)
  
  p2<-p2+labs(x='cell ordering',y='PC1')
  # p2<-p2+scale_color_manual(name='cluster',
  #                           breaks=label_names,
  #                           values = label_color,
  #                           labels=label_names)
  p2<-p2+geom_point(mapping =aes(x=x_center,y=y_center),data = h_after_layout,size=4,color=h_after_layout[,"stage_colors"])
  h_edge=get.edgelist(hh)
  segment_data=data.frame(x=h_after_layout[h_edge[,1],1],
                          y=h_after_layout[h_edge[,1],2],
                          xend=h_after_layout[h_edge[,2],1],
                          yend=h_after_layout[h_edge[,2],2])
  p2<-p2+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                      size=1)
  ##############################
  #######################
  #p3<-ggplot(gg_data)
  #p3<-p3+geom_point(mapping = aes(x=x,y=z),data = gg_data,color=gg_data[,"stage_colors"])
  #h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
  #colnames(h_after_layout)=c('x','y','z','stage_colors')
  #txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
  #gg_data=rbind(gg_data,h_after_layout)
  p3<-ggplot(gg_data,mapping = aes(x=x,y=z))
  p3<-p3+geom_point(data = gg_data,color=gg_data[,"stage_colors"])+geom_text(aes(label=txt_lab),hjust=0,vjust=0,size=20)
  
  p3<-p3+labs(x='cell ordering',y='PC2')
  # p3<-p3+scale_color_manual(name='cluster',
  #                           breaks=label_names,
  #                           values = label_color,
  #                           labels=label_names)
  p3<-p3+geom_point(mapping =aes(x=x_center,y=z_center),data = h_after_layout,size=4,color=h_after_layout[,"stage_colors"])
  h_edge=get.edgelist(hh)
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
  
  Results$orderP1=p1
  Results$orderP2=p2
  Results$orderP3=p3
  return(Results)
}