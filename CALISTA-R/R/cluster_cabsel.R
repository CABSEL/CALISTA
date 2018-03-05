cluster_cabsel <- function(mRNA_all,log_p_mat_ma,K_new,range,...){
  #test here
  #cluster_cabsel(mRNA_all,log_p_mat_ma,K_new,expected_clusters,'parallel',1,
  #               'pop_size',loops,'get_k',t(final_groups))
  #my_results=cluster_cabsel(mRNA_all,log_p_mat_ma,K_new,max_clust,'parallel',1,
  # pop_size',loops,'cluster','kmedoids')
  #range=max_clust
  #varagin=list('parallel',1,'pop_size',loops,'cluster','kmedoids')
  #my_results=cluster_cabsel(mRNA_all,log_p_mat_ma,K_new,max_clust,'parallel',1,'pop_size',loops,'cluster','kmedoids')
  
  
  
  #mRNA_all=as.matrix(mRNA_all)
  varagin=list(...)
  if(nargs()<4){
  stop('Not enough input variables.')
  }
  source('greedy_cabsel.R')
  min_C=min(range)
  max_C=max(range)
  expected_clusters=min_C:max_C
  as_all=list()
  as_get_k=NULL
  optimize=TRUE
  in_pop=FALSE
  idx=NULL
  loops=10
  run_parallel=FALSE
  algorithm='greedy_cabsel'
  max_iter=100
  nvars=nrow(mRNA_all)
  n_genes=ncol(mRNA_all)
  nVarargs=length(varagin)
  my_plot=TRUE
  my_cluster='kmedoids'
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
    if(varagin[[k]]=='parallel'){
      run_parallel=varagin[[k+1]]
    }
    if(varagin[[k]]=='plot'){
      my_plot=varagin[[k+1]]
    }
    if(varagin[[k]]=='cluster'){
      my_cluster=varagin[[k+1]]
    }
    if(varagin[[k]]=='algorithm'){
      algorithm=varagin[[k+1]]
    }
  }
  
 if(optimize==FALSE){
   as_all=as_get_k
   if(is.null(expected_clusters)){
     expected_clusters=length(unique(as_all))
   }
 }else if(length(expected_clusters==1)){
     set.seed(round(nvars/2*expected_clusters))
     as_all=matrix(sample(expected_clusters,loops*nvars,replace = TRUE),nrow = loops)
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
    p=0
  }else{
    p=6
    # if(is.null(p)){
    #   cl=makeCluster(6)
    #   registerDoParallel(cl)
    # }
    # p=getDoParWorkers()
  }
  
  source('greedy_cabsel.R')
  # run greedy algorithm
  display=1
  my_results=greedy_cabsel(as_all,log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,
                           optimize,opt_idx_a,p,sum_prob_tot,in_population,loops,expected_clusters,algorithm,display)
  my_results$cluster=expected_clusters[length(expected_clusters)]
  population=my_results$population
  source('get_consensus.R')
  consensus=get_consensus(population); # implement get_consensus here
  all_my_results=list()
  all_my_results$all=my_results
  i=0
  final_results=list()
  if(optimize ==TRUE &&length(expected_clusters)>1){
    if(my_cluster=='hierarchical'){
      Z=hclust(scale(consensus),method = 'complete')
    }
    idx_total=list()
    silh=list()
    for (sens in minC:max_C) {
      i=i+1
      switch (my_cluster,
        #'kmedoids' = {idx=kmedoids(consensus,sens,'Replicates',5,'Alogorithm','pam')},
        'kmedoids' = {idx=pam(consensus,sens)},  #notice here, replicates time
        'hierarchical'={idx=cutree(Z,k=sens)}   # notich here, might not right
      )
      idx_total[[i]]=idx
      all_my_results[[i]]$all$idx=idx
      all_my_results[[i]]$all$cluster=sens
      # omit silhouette here
      s=silhouette(idx,scale(consensus))
      silh[i]=mean(s)
      all_my_results[[i]]$all$silh=silh[is]
    }
    idx_total=t(data.frame(idx_total))
    # omit plot silhouette here
    if(my_plot==TRUE){
      window()
      plot(expected_clusters,silh,type = '*',col='blue',xlab = 'cluster',ylab = 'Silhouette')
      title('Silhouette')
    }
    #omit final_results here 
    idx_silh_max=which(silh=max(silh))
    final_results$opt_number_of_clusters_silh=expected_clusters[idx_silh_max]
    final_results$final_groups_silh=idx_total[idx_silh_max,]
  }
  
  source('greedy_cabsel.R')
  if(length(expected_clusters==1)&& optimize==TRUE){
    switch (my_cluster,
      'kmedoids' = {idx=pam(consensus,expected_clusters)$clustering},
      'hierarchical'={Z=hclust(scale(consensus),method = 'complete')
                      idx=cutree(Z,k=expected_clusters)}      
    )
  
    my_results_1=greedy_cabsel(t(idx),log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,
                               FALSE,opt_idx_a,p,sum_prob_tot,in_population,1,expected_clusters,algorithm) #here,loops can only be one
    all_my_results$all$clusterprobabilities=my_results_1$clusterprobabilities
    all_my_results$all$idx=idx
  }

  final_results$all=all_my_results
  return(final_results)
}