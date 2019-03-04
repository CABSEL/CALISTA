
#' A calista Function
#'
#' @param mRNA_all,log_p_mat_ma,K_new,expected_clusters,...
#' @keywords calista
#' @export
#' @examples
#' CALISTA_clustering()






CALISTA_clustering <- function(mRNA_all,log_p_mat_ma,K_new,expected_clusters,...){

  #my_results_c=CALISTA_clustering(DATA$totDATA,Parameters$Parameters[[3]],
  #Parameters$Parameters[[1]],Results$expected_clusters,
  #'parallel',parallel,'cluster',Cluster,'algorithm',algorithm)
  #test
   # mRNA_all=DATA$totDATA
   # log_p_mat_ma=Parameters$Parameters[[3]]
   # K_new=Parameters$Parameters[[1]]
   # expected_clusters=Results$expected_clusters
   # varagin=list('parallel',parallel,'cluster',Cluster,'algorithm',algorithm)

  varagin=list(...)
  if(nargs()<4){
    stop('Not enough input variables.')
  }
  if (length(expected_clusters)>1){
    stop('Please input only one value for the expected number of clusters')
  }
  ##default values
  #fun optimizization
  optimize=TRUE
  #initial population
  in_pop=FALSE
  #number of runs
  loops=50 # change to 50 later !!!!!!!!!
  #run in parallel
  run_parallel=FALSE
  #cluter algorithm
  algorithm='greedy_cabsel'
  #num of maximum iterations
  max_iter=100
  #number of cells and number of genes
  nvars=nrow(mRNA_all)
  n_genes=ncol(mRNA_all)
  nVarargs=length(varagin)
  #plot
  my_plot=TRUE
  #clustering
  my_cluster='kmedoids'
  as_all=list()
  as_get_k=NULL
  for (k in seq(from=1,to=nVarargs,by=2)) {
    if(varagin[[k]]=='get_k'){
      as_get_k=varagin[[k+1]]
      optimize=FALSE
      loops=1
      if(nrow(as_get_k)>1){
        stop('Please provide only 1 initial population in order to calculate the rate constants wihout optimization')
      }
    }
    if(varagin[[k]]=='set_initial_pop'){
      in_pop=TRUE
      as_in=varagin[[k+1]]
    }
    if(varagin[[k]]=='pop_size'){
      loops=varagin[[k+1]]
    }
    if(varagin[[k]]=='max_iteration'){
      max_iter=varagin[[k+1]]
    }
    if(varagin[[k]]=='parallel'){
      run_parallel=varagin[[k+1]]
    }
    if(varagin[[k]]=='algorithm'){
      algorithm=varagin[[k+1]]
    }
    if(varagin[[k]]=='cluster'){
      my_cluster=varagin[[k+1]]
    }
  }
  if(optimize==FALSE){
    as_all=as_get_k
    if(is.null(expected_clusters)){
      expected_clusters=length(unique(as_all))
    }
  }else if(length(expected_clusters==1)){
    set.seed(round(nvars/20*expected_clusters))
    as_all=matrix(sample(1:expected_clusters,loops*nvars,replace = TRUE),nrow = loops)
  }

  ### if user provides initial poplulation and at the same time want to
  ### calculate the optimal constants, throw warning
  if (in_pop==TRUE){
    if(optimize==FALSE){
      stop('No optimization is carried out, since option "get_k" was set')
    }else{
      if(nrow(as_in)<loops){
        loops=nrow(as_in)
        warning('Number of iteratios was set according to size of the initial population')
        as_all=as_in
      }
    }
  }
  #define variable used in the function greedy_cabsel
  in_population=matrix(0,loops,nvars)
  sum_prob_tot=matrix(0,1,loops)
  opt_idx=list(run=NULL)
  opt_idx_a=list()
  opt_idx_a[1:loops]=opt_idx
  p=NULL
  if(run_parallel==FALSE){
    p=1
  }else{
    p=detectCores()-1
    # if(is.null(p)){
    #   cl=makeCluster(6)
    #   registerDoParallel(cl)
    # }
    # p=getDoParWorkers()
  }
  # run greedy algorithm
  display=1
  # in_population=matrix(0,loops,nvars)
  # sum_prob_tot=matrix(0,1,loops)
  # opt_idx=list(run=NULL)
  # opt_idx_a=list()
  # opt_idx_a[1:loops]=opt_idx
  # p=NULL
  # if(run_parallel==FALSE){
  #   p=1
  # }else{
  #   p=detectCores()-1
  #   # if(is.null(p)){
  #   #   cl=makeCluster(6)
  #   #   registerDoParallel(cl)
  #   # }
  #   # p=getDoParWorkers()
  # }
  #num_clusters=length(expected_clusters)

  #write greedy_cabsel to get my_results here
  #here notice as_all
  #source('greedy_cabsel.R')

  my_results=greedy_cabsel_ff(as_all,log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,optimize,opt_idx_a,p,sum_prob_tot,in_population,loops,expected_clusters,algorithm,display)
  my_results$cluster=expected_clusters
  population=my_results$population
  #consensus
  #source('get_consensus.R')
  consensus=get_consensus(population) # implement get_consensus here
  all_my_results=list()
  all_my_results$all=my_results
  i=0
  final_results=list()
  #source('greedy_cabsel.R')

  if(optimize==TRUE){
    switch (my_cluster,
            # 'kmedoids' = {idx=pam(consensus,expected_clusters)$clustering},
            'kmedoids' = {dissimilarity_matrix=1-consensus/loops
            idx=as.integer(Cluster_Medoids(dissimilarity_matrix,clusters = expected_clusters,
                                           distance_metric = "euclidean",seed = 1)$clusters)},
            #'kmedoids' =
            # print(dim(consensus))
            #  idx=Cluster_Medoids(consensus,clusters = expected_clusters, distance_metric = "euclidean",seed = 10)$clusters
            #  idx=as.integer(idx)
            #  print(idx)
            #  str(idx)},
            'hierarchical'={Z=hclust(dist(scale(consensus)),method = 'complete')
            idx=cutree(Z,k=expected_clusters)}
    )
    }else{idx=as_all}
    display=0


    my_results_1=greedy_cabsel_f(t(idx),log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,FALSE,p,
                                 sum_prob_tot,population,1,expected_clusters,algorithm,display,opt_idx_a) #here,loops can only be one

    all_my_results$all$clusterprobabilities=my_results_1$clusterprobabilities
    all_my_results$all$idx=idx
    all_my_results$all$cell_prob = my_results_1$cell_prob
    # my_distance=my_results_1$cell_prob#matrix(0,nvars,expected_clusters)
    # for (m in 1:expected_clusters) {
    #   aaa=which(idx==m)
    #   my_distance[aaa,]=abs(rep( my_results_1$cell_prob[aaa,m],1*expected_clusters)- my_results_1$cell_prob[aaa,])   #check heres
    # }
    all_my_results$all$distance =  my_results_1$distance
    all_my_results$all$param_idx = my_results_1$param_idx
    all_my_results$all$consensus = consensus



  final_results$all=all_my_results
  return(final_results)
  # return(idx)
}


















