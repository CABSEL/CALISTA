#   %CALISTA_TRANSITION_MAIN infer lineage progression among cell clusters
#   % CALISTA uses the cluster distances - a measure of cluster-cluster
#   % dissimilarity - to infer the lineage progression or cluster-cluster
#   % relationship. The lineage progression graph is built based on adding
#   % edges between clusters in increasing magnitude of cluster distance.
#   % CALISTA also provides a simple user-interface to add and remove edges
#   % based on the cluster distances. 
#   % 
#   % Usage:
#   % Run CALISTA lineage inference using CALISTA single-clustering result
#   % Results=CALISTA_transition_main(DATA,Results);
#   % 
#   % Run CALISTA lineage Inference with user-defined cell clusters 
#   % Results=list()
#   % Results=CALISTA_transition_main(DATA,Results,cell_assignments)
#   %    
#   % Inputs:   
#   % DATA - a structure containing preprocessed single-cell expression data
#   % Use 'import_data' to upload and preprocess single-cell expression values
#   %
#   % Results - a structure of CALISTA clustering results
#   % Run 'CALISTA_clustering_main'
#   %
#   % cell_assignments - 1xN vector of INTEGERS with N = number of cells. 
#   % The n-th element of cell_assignments contains the cluster assignment of
#   % the n-th cell of the expression data uploaded. Cluster names must
#   % be assigned in sequence (e.g. 1,2,3,4 and not 1,2,4).
#   %    
#   % Outputs:
#   % Results - a structure containing the results of CALISTA analysis. 
#   % The most relevant field containing the clustering analysis is:
#   %
#   % Results$TRANSITION$nodes_connection - an Ex2 matrix containing E edges 
#   % of the lineage progression graph. The first and second columns give the
#   % clusters incident to the edges. If time/stage information is available,
#   % then the edges are directed where the first column indicates the source
#   % cluster and the second column indicates the target cluster. The user can
#   % also specify the starting cell or marker gene to determine the direction
#   % of the edges. 
#   %
#   % Results$TRANSITION$cluster_distance - a ranked list of cluster distances. 
#   % The first and second columns are cluster indices and the thrid column 
#   % gives the cluster distance.
#   %
#   % Created by Nan Papili Gao (R version implemented by Tao Fang)
#   %            Institute for Chemical and Bioengineering 
#   %            ETH Zurich
#   %            E-mail:  nanp@ethz.ch
#   %
#   % Copyright. June 1, 2017.

CALISTA_transition_main<-function(DATA,Results,cell_assignments){
  if (nargs()<2){
    stop("Not enough input arguments")
  }
  
  if (length(Results)==0){
    if (nargs()<3){
      stop("Not enough input arguments")
    }
    Results=jump_clustering(DATA,cell_assignments)
  }
  
  Results=CALISTA_transition(Results,DATA)
  
  ###3.b- plot mean expression for each cluster
  Results=Plot_cluster_mean_exp(Results,DATA)
  
  ###3.c- cell-celll variability analysis
  Results=cell_variability(Results,DATA)
  return(Results)
}