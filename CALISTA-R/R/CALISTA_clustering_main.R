#   %CALISTA_CLUSTERING_MAIN single-cell clustering in CALISTA
#   % CALISTA implements two-step clustering algorithm. The first step involves
#   % a likelihood-based clustering based on the stochastic two-state model of
#   % gene transcriptional process. A greedy algorithm is implemented to obtain
#   % cell clustering with maximizes the overall likelihood, and repeatedly run
#   % to produce a consensus matrix (i.e. the matrix of the number of times two
#   % cells are clustered together). The second and final step involves
#   % k-medoids or hierachical clustering using the consensus matrix. 
#   %
#   % Usage:
#   % Results=list()
#   % CALISTA_clustering_main_results=CALISTA_clustering_main(DATA,INPUTS)
#   % Results=CALISTA_clustering_main_results$Results
#   % DATA=CALISTA_clustering_main_results$DATA
#   % INPUTS=CALISTA_clustering_main_results$INPUTS
#   % Perform single-cell clustering using user-defined input values
#   % 
#   % Inputs:
#   % INPUTS - a structure containing settings for the single-cell clustering 
#   %
#   % ** INPUTS$optimize **
#   % 1- select the number of clusters by eigengap plot
#   % 0- define the number of clusters 
#   %
#   % ** INPUTS$parallel **
#   % 1- Use parallel processing (number of cores available - 1)  
#   % 0- Do not use processing
#   %
#   % ** INPUTS$max_iter **
#   % Maximum number of iterations in the greedy algorithm 
#   % INPUTS$max_iter = 100; 
#   % 
#   % ** INPUTS$runs **
#   % Number of clustering runs 
#   % INPUTS$runs = 50; 
#   %
#   % ** INPUTS$cluster **
#   % 'hierarchical'- hierachical clustering of consensus matrix
#   % 'kmedoids'-  kmedoids for the clustering of consensus matrix
#   % 
#   % DATA - a structure containing preprocessed single-cell expression data
#   % Use 'import_data' to upload and preprocess single-cell expression values
#   %
#   % Outputs:
#   % Results - a structure containing the results of CALISTA analysis.
#   % The most relevant fields containing the clustering analysis are:
#   %
#   % Results$final_groups - vector containing cluster assignment for the cells
#   %
#   % Results$singleCELLclusterDATA - a 1xK cell array where K is the number of 
#   % clusters. Each cell contains a GxNk (normalised) expression matrix of G 
#   % genes and Nk cells in the k-th cluster.
#   %
#   % Created by Nan Papili Gao (R version implemented by Tao Fang)
#   %            Institute for Chemical and Bioengineering 
#   %            ETH Zurich
#   %            E-mail:  nanp@ethz.ch
#   %
#   % Copyright. June 1, 2017.
#  

CALISTA_clustering_main<-function(DATA,INPUTS){
  if(nargs()<2){
    stop("Not enough input arguments")
  }
  INPUTS$algorithm='greedy_cabsel'
  Parameters=DATA$Parameters
  optimize=INPUTS$optimize
  algorithm=INPUTS$algorithm
  Cluster=INPUTS$Cluster
  parallel=INPUTS$parallel
  loops=INPUTS$runs;
  max_iter=INPUTS$max_iter;
  
  Results=list()
  
  if(optimize){
    #print(optimize)
    max_clust=12
    my_results=CALISTA_clustering(DATA$totDATA,Parameters$Parameters[[3]],Parameters$Parameters[[1]],max_clust,'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops)
    consensus=get_consensus(my_results$all$all$population)
    Results$expected_clusters=eigengap(consensus,max_clust)
  }else{
    writeLines('Number of clusters expected:')
    Results$expected_clusters=readLines(n=1)
    Results$expected_clusters=as.integer(Results$expected_clusters)
  }
  
  ###2.b-cell clustering
  my_results_c=CALISTA_clustering(DATA$totDATA,Parameters$Parameters[[3]],Parameters$Parameters[[1]],Results$expected_clusters,
                                  'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops)
  Results$final_groups=t(my_results_c$all$all$idx)
  
  ###2.c- relabeling based on time/cell stage info
  method=3
  find_progression2_list=find_progression2(Results,DATA,method)
  Results=find_progression2_list$Results
  DATA=find_progression2_list$DATA
  
  ###2.d- optimal parameter estimation for the final cluster assignment
  my_results_final=CALISTA_clustering(DATA$totDATA,Parameters$Parameters[[3]],Parameters$Parameters[[1]],Results$expected_clusters,
                                      'parallel',parallel,'get_k',Results$final_groups,'algorithm',algorithm)
 
  ###2.e- cluster visualization
  reduction=2;
  Results=visualization(reduction,DATA,Results)
  Results$clustering_struct=my_results_final
 
  CALISTA_clustering_main_results=list()
  CALISTA_clustering_main_results$Results=Results
  CALISTA_clustering_main_results$DATA=DATA
  CALISTA_clustering_main_results$INPUTS=INPUTS
  return(CALISTA_clustering_main_results)
}
  
  
