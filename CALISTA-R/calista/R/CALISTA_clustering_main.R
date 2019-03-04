#' A calista Function
#'
#' @param DATA,INPUTS
#' @keywords calista
#' @export
#' @examples
#'
#'
#' CALISTA_clustering_main()

CALISTA_clustering_main<-function(DATA,INPUTS){
  if(nargs()<2){
    stop("Not enough input arguments")
  }

  ##CHECK INPUT ARGUMENTS
  if (!exists('data_type', where=INPUTS)){
    stop('Please specify the data type in INPUTS.data_type')
  }

  if  (!exists('format_data', where=INPUTS)){
    stop('Please specify the format of the data in INPUTS$format_data')
  }

  #if  (!exists('cluster_time', where=INPUTS)){
  #INPUTS$cluster_time=0
  #}

  if  (!exists('cells_2_cut', where=INPUTS)){
    INPUTS$cells_2_cut=0;
  }

  if  (!exists('cut_variable_genes', where=INPUTS)){
    INPUTS$cut_variable_genes=1000
  }

  if(!exists('data_selection', where=INPUTS)){
    INPUTS$data_selection=integer()
  }

  #  if (INPUTS$cluster_time~=0 && ~INPUTS$data_selection==0){
  #  writelines('\nCALISTA Time CLustering is active. INPUTS$data_selection set as 0. All time points are considered for the analysis\n')
  #  INPUTS$data_selection=0
  #}

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






  INPUTS$algorithm='greedy_cabsel'
  Parameters=DATA$Parameters
  optimize=INPUTS$optimize
  algorithm=INPUTS$algorithm
  Cluster=INPUTS$Cluster
  parallel=INPUTS$parallel
  loops=INPUTS$runs;
  max_iter=INPUTS$max_iter;

  Results=list()

  if (exists('cell_assignments', where=INPUTS)){
    writeLines('Using user defined cell assignments..Jumping CALISTA clustering')
    cell_assignments=INPUTS$cell_assignments
    if (nrow(cell_assignments)>ncol(cell_assignments)) {
      cell_assignments=t(cell_assignments)
    }

    # % Create new Results structure
    Results=list()
    Results$cluster_predicted=unique(t(cell_assignments));
    Results$expected_clusters=length(Results$cluster_predicted);
    Results$final_groups=cell_assignments;
  }else{


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
    time1=Sys.time()
    ###2.b-cell clustering
    my_results_c=CALISTA_clustering(DATA$totDATA,Parameters$Parameters[[3]],Parameters$Parameters[[1]],Results$expected_clusters,
                                    'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops)
    Results$final_groups=t(my_results_c$all$all$idx)
    time2=Sys.time()
    Results$runtime=time2-time1
  }
  ###2.c- relabeling based on time/cell stage info
  method=3
  find_progression2_list=find_progression2(Results,DATA,method)
  Results=find_progression2_list$Results
  DATA=find_progression2_list$DATA

  ###2.d- optimal parameter estimation for the final cluster assignment
  my_results_final=CALISTA_clustering(DATA$totDATA,Parameters$Parameters[[3]],Parameters$Parameters[[1]],Results$expected_clusters,
                                      'parallel',parallel,'get_k',Results$final_groups,'algorithm',algorithm)


  Results$clustering_struct=my_results_final
  if  (exists('my_results_c')){
    Results$clustering_struct$all$all$population=my_results_c$all$all$population
    if  (exists('consensus', where=my_results_c$all$all)){
      Results$clustering_struct$all$all$consensus=my_results_c$all$all$consensus}}

  #Results$clustering_struct$all$all$parameter=my_results_final$all$all$parameter
  #Results$clustering_struct$all$all$clusterprobabilities=my_results_final$all$all$clusterprobabilities
  #Results$clustering_struct$all$all$distance=my_results_final$all$all$distance
  #Results$clustering_struct$all$all$cell_prob=my_results_final$all$all$cell_prob
  #Results$clustering_struct$all$all$param_idx=my_results_final$all$all$param_idx

  ###2.e- cluster visualization
  reduction=2;
  Results=visualization(reduction,DATA,Results)


  CALISTA_clustering_main_results=list()
  CALISTA_clustering_main_results$Results=Results
  CALISTA_clustering_main_results$DATA=DATA
  CALISTA_clustering_main_results$INPUTS=INPUTS

  return(CALISTA_clustering_main_results)
}


