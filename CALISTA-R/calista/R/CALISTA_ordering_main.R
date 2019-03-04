#' A calista Function
#'
#' @param DATA,Results,cell_assignments
#' @keywords calista
#' @export
#' @examples
#'
#'
#' CALISTA_ordering_main()



CALISTA_ordering_main<-function(DATA,Results,cell_assignments){

  if(nargs()<2){
    stop('Not enough input variables.')
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
  Z_temp=Parameters$Parameters[[2]]

  cell_ordering=matrix(0,DATA$nvars,2)
  placement=matrix(0,DATA$nvars,1)
  edges_of_cell=matrix(0,DATA$nvars,2)
  n_points_interp=100

  h=Results$TRANSITION$final_graph
  nodes_connection3=get.edgelist(h)
  n_edges=nrow(nodes_connection3)
  max_log_P_cell_all=matrix(0,DATA$nvars,1)
  # colorMARK2=Results$colorMARK_calista
  # cell_prob_3D=get_3D_cell_prob(Results,Parameters,DATA)
  # h_edge=nodes_connection3

  for (clust in 1:Results$expected_clusters){

    actual_cluster=clust
    edges_2_check=which(nodes_connection3[,1]==actual_cluster | nodes_connection3[,2]==actual_cluster)
    idx_cells_actual_clust=which(Results$final_groups %in% clust)
    max_log_P_each_edge=matrix(0,length(edges_2_check),length(idx_cells_actual_clust))
    idx_max_log_P_each_edge=matrix(0,length(edges_2_check),length(idx_cells_actual_clust))
    time_max_log_P_each_edge=matrix(0,length(edges_2_check),length(idx_cells_actual_clust))
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
      sum_log_P_cells=matrix(0,n_points_interp,length(idx_cells_actual_clust))
      for(j in 1:DATA$numGENES){
        int1=Z_temp[, Results$clustering_struct$all$all$param_idx[clust1,j]]
        int2=Z_temp[, Results$clustering_struct$all$all$param_idx[clust2,j]]
        M=cbind(int1,int2)
        Xnew=apply(M,1, function(x) approx(c(t1,t2),x,xout=t_interp))
        exp_cells_actual_clust_actual_gene=DATA$totDATA[idx_cells_actual_clust,j]+1
        xnew=NULL
        for (ii in 1:length(Xnew)){
          xnew=cbind(xnew,Xnew[[ii]][[2]])
        }
        sum_log_P_cells=sum_log_P_cells+log2(xnew[,exp_cells_actual_clust_actual_gene])
      }
      max_log_P_each_edge[i,]=apply(sum_log_P_cells, 2, function(x) max(x))
      idx_max_log_P_each_edge[i,]=apply(sum_log_P_cells, 2, function(x) which.max(x))
      if (same_t){
        time_max_log_P_each_edge[i,]=tequi
      }else{
        time_max_log_P_each_edge[i,]=t_interp[idx_max_log_P_each_edge[i,]]
      }

    }
    if (length(edges_2_check)==1){
      max_log_P_cell=max_log_P_each_edge
      idx_max_log_P_cell=matrix(1,1,length(max_log_P_cell))
      cell_ordering[idx_cells_actual_clust,1]=time_max_log_P_each_edge
      placement[idx_cells_actual_clust]=idx_max_log_P_each_edge
    }else{
      max_log_P_cell=apply(max_log_P_each_edge, 2, function(x) max(x))
      idx_max_log_P_cell=apply(max_log_P_each_edge, 2, function(x) which.max(x))

      A=t(time_max_log_P_each_edge)
      A2=t(idx_max_log_P_each_edge)
      J=idx_max_log_P_cell
      I = 1 : nrow(A)
      M=cbind(I,J)
      C=apply(M,1, function(x) A[x[1],x[2]])
      C2=apply(M,1, function(x) A2[x[1],x[2]])
      cell_ordering[idx_cells_actual_clust,1] = C
      placement[idx_cells_actual_clust] = C2
    }

    max_log_P_cell_all[idx_cells_actual_clust]=max_log_P_cell
    cell_ordering[idx_cells_actual_clust,2]=edges_2_check[idx_max_log_P_cell]
    edges_of_cell[idx_cells_actual_clust,]=nodes_connection3[cell_ordering[idx_cells_actual_clust,2],]
  }
  normed_cell_ordering=(cell_ordering[,1]-min(cell_ordering[,1]))/(max(cell_ordering[,1]-min(cell_ordering[,1])))

  x_center=matrix(0,Results$expected_clusters,1)
  y_center=Results$TRANSITION$y_center
  z_center=Results$TRANSITION$z_center
  for( i in 1:Results$expected_clusters){
    x_center[i]=mean(normed_cell_ordering[Results$final_groups==i])
  }

  XData=x_center
  YData=y_center
  ZData=z_center

# Place cells on the progression graph

  #Ycoordinate
  int1=YData[edges_of_cell[,1]]
  int2=YData[edges_of_cell[,2]]
  M=cbind(int1,int2)
  t1=1
  t2=100
  t_interp=t1:t2
  Xnew=apply(M,1, function(x) approx(c(t1,t2),x,xout=t_interp))
  xnew=NULL
  for (ii in 1:length(Xnew)){
    xnew=cbind(xnew,Xnew[[ii]][[2]])
  }
  A=t(xnew)
  J=placement
  I = 1 : nrow(A)
  M=cbind(I,J)
  Ycoordinate=apply(M,1, function(x) A[x[1],x[2]])

  #Xcoordinate
  int1=XData[edges_of_cell[,1]]
  int2=XData[edges_of_cell[,2]]
  M=cbind(int1,int2)
  Xnew=apply(M,1, function(x) approx(c(t1,t2),x,xout=t_interp))
  xnew=NULL
  for (ii in 1:length(Xnew)){
    xnew=cbind(xnew,Xnew[[ii]][[2]])
  }
  A=t(xnew)
  M=cbind(I,J)
  Xcoordinate=apply(M,1, function(x) A[x[1],x[2]])

  #Zcoordinate
  int1=ZData[edges_of_cell[,1]]
  int2=ZData[edges_of_cell[,2]]
  M=cbind(int1,int2)
  Xnew=apply(M,1, function(x) approx(c(t1,t2),x,xout=t_interp))
  xnew=NULL
  for (ii in 1:length(Xnew)){
    xnew=cbind(xnew,Xnew[[ii]][[2]])
  }
  A=t(xnew)
  M=cbind(I,J)
  Zcoordinate=apply(M,1, function(x) A[x[1],x[2]])



  phi = runif(DATA$nvars, min=0, max=2*pi)
  costheta=runif(DATA$nvars, min=-1, max=1)
  theta = acos(costheta)
  radii = 0.3*(runif(DATA$nvars)^(1/3))
  Xnoise=radii * sin( theta) * cos( phi )
  Ynoise=radii * sin( theta) * sin( phi )
  Znoise=radii * cos( theta)

  Xcoordinate=Xcoordinate+Xnoise/8;
  Ycoordinate=Ycoordinate+Ynoise
  Zcoordinate=Zcoordinate+Znoise
  Results$ORDERING$Xcoordinate=Xcoordinate
  Results$ORDERING$Ycoordinate=Ycoordinate
  Results$ORDERING$Zcoordinate=Zcoordinate
  Results$ORDERING$x_center=x_center


  groups=Results$final_groups
  colorMARK=Results$colorMARK_calista
  group_names='cluster'
  Results$ORDERING$p_cluster = plot_4_calista_ordering(groups,colorMARK, group_names, Results)



  groups=DATA$timeline
  colorMARK=Results$colorMARK_time
  group_names='time'
  Results$ORDERING$p_time = plot_4_calista_ordering(groups,colorMARK, group_names, Results)



  groups=normed_cell_ordering
  colorMARK=normed_cell_ordering
  group_names='pseudotime'
  Results$ORDERING$p_pseudotime = plot_4_calista_ordering(groups,colorMARK, group_names, Results)


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



  return(Results)
}
