#   %CALISTA_PATH_MAIN perform post-analysis along developmental path(s) 
#   % This function gives transition genes, performs hiearchical clustering,
#   % smoothen gene expression and constructs gene co-expression network for
#   % each path defined by the users.
#   %
#   % Usage:
#   %
#   % Results=CALISTA_path_main(INPUTS,Results)
#   %
#   % CALISTA will ask the user to specify the paths for the post-analysis.
#   % Each path consists of a sequence of clusters in the lineage progression
#   % graph. 
#   %
#   % Inputs:
#   % INPUTS - a structure containing the settings for the post-analysis with
#   % following fields:
#   %
#   % INPUTS$plot_fig
#   % 1- plot figure of smoothened gene expression along path
#   % 0- otherwise
#   %
#   % INPUTS$hclustering
#   % 1- perform hierarchical clustering of gene expression for each path 
#   % 0- otherwise
#   %
#   % INPUTS$method
#   % 1- use partial correlation for the gene co-expression network (value_cutoff=0.4, pvalue_cutoff=0.05)
#   % 2- use pairwise correlation for the gene co-expression network (value_cutoff=0.8, pvalue_cutoff=0.01)
#   %
#   % INPUTS$moving_average_window - set the size of window (percent of cells in each path)
#   % used for the moving averaging. 
#   %
#   % Results - a structure of CALISTA clustering,lineage inference, transition genes detection and cell ordering results
#   % Run 'CALISTA_clustering_main',
#   %     'CALISTA_transition_main',
#   %     'CALISTA_transition_genes_main', and
#   %     'CALISTA_cell_ordering_main'
#   % 
#   % Outputs:
#   % Results$PATH$path_transition_genes - a 1xP cell array of P paths. The 
#   % i-th cell contains the transition genes of the i-th path.
#   % 
#   % Results$PATH$smoothExpr - a 1xP cell array of P paths. The i-th cell
#   % gives the moving-averaged single-cell expression matrix where each row
#   % corresponds to a gene.
#   %
#   % Created by Nan Papili Gao (R version implemented by Tao Fang)
#   %            Institute for Chemical and Bioengineering
#   %            ETH Zurich
#   %            E-mail:  nanp@ethz.ch
#   %
#   % Copyright. June 1, 2017.


CALISTA_path_main<-function(INPUTS,Results){
  if(nargs()<2){
    stop("Not enough input arguments")
  }
  # % ** CV_plot **
  #   % 1- calculate and plot the coefficients of variation
  # % 0- otherwise (DEFAULT)
  INPUTS$CV_plot=0
  Results=CALISTA_path(Results,INPUTS)
  Results=CALISTA_net_path(Results,INPUTS)
  
  return(Results)
}