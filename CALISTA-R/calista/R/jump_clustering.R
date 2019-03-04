#' A calista Function
#'
#' @param DATA,cell_assignments
#' @keywords calista
#' @export
#' @examples
#' jump_clustering()


jump_clustering <- function(DATA,cell_assignments){

  algorithm='greedy_cabsel'
  Parameters=DATA$Parameters
  parallel=0

  if (nrow(cell_assignments)>ncol(cell_assignments)) {
    cell_assignments=t(cell_assignments)
  }

  # % Create new Results structure
  Results=list()
  Results$cluster_predicted=unique(t(cell_assignments));
  Results$expected_clusters=length(Results$cluster_predicted);
  Results$final_groups=cell_assignments;

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

  return(Results)
}
