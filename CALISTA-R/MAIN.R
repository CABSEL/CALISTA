##############################################################
##                         MAIN SCRIPT                      ##
##############################################################
## Before running CALISTA:
## - please set CALISTA main folder as working directory (with setwd command) and install the required packages
## - For any questions please check the README file or write to nanp@ethz.ch
rm(list=ls())
graphics.off()
options(warn=-1)

source("./R/initialization.R")

# %% *** 1-Data Import and Preprocessing ***
# %
# % Please check comments in 'import_data' for more information 
# %
# % Specify data types and settings for pre-processing
INPUTS=list()
INPUTS$data_type=1 # Single-cell RT-qPCR CT data
INPUTS$format_data=1 # Rows= cells and Columns= genes with time/stage info in the last column
INPUTS$data_selection=integer() # Include data from all time points
INPUTS$perczeros_genes=100 # Remove genes with > 100% of zeros 
INPUTS$perczeros_cells=100 # Remove cells with 100% of zeros
INPUTS$cells_2_cut=0 # No manual removal of cells
INPUTS$perc_top_genes=100 # Retain only top X the most variable genes with X=min(200, INPUTS.perc_top_genes * num of cells/100, num of genes)

# % Specify single-cell clustering settings
INPUTS$optimize=0 # The number of cluster is known a priori
INPUTS$parallel=1 # Use parallelization option
INPUTS$runs=50 # Perform 50 independent runs of greedy algorithm 
INPUTS$max_iter=100 # Limit the number of iterations in greedy algorithm to 100
INPUTS$Cluster='kmedoids' # Use k-medoids in consensus clustering

# % Specify transition genes settings
INPUTS$thr_transition_genes=50 # Set threshold for transition genes determination to 50%

# % Specify path analysis settings
INPUTS$plot_fig=1; # Plot figure of smoothened gene expression along path
INPUTS$hclustering=1; # Perform hierarchical clustering of gene expression for each path 
INPUTS$method=2; # Use pairwise correlation for the gene co-expression network (value_cutoff=0.8, pvalue_cutoff=0.01)
INPUTS$moving_average_window=10; # Set the size of window (percent of cells in each path) used for the moving averaging

# % Upload and pre-process data 
DATA=import_data(INPUTS)

# %% *** 2-SINGLE-CELL CLUSTERING ***
# %
# % Please check comments in 'CALISTA_clustering_main' for more information. 
# %

Results=list()
CALISTA_clustering_main_results=CALISTA_clustering_main(DATA,INPUTS)
Results=CALISTA_clustering_main_results$Results
DATA=CALISTA_clustering_main_results$DATA
INPUTS=CALISTA_clustering_main_results$INPUTS

###cluster cutting 
writeLines('Press 1 if you want to remove one cell cluster, 0 otherwise:')
cluster_cut=readLines(n=1)
cluster_cut=as.integer(cluster_cut)

if(cluster_cut==1){
  writeLines('key cluster number(e.g 1 or 5 3)')
  which_cut=strsplit(readLines(n=1), " ")
  which_cut=as.integer(unlist(which_cut))
  cells_2_cut2=which(Results$final_groups %in% which_cut)
  cells_2_cut2=matrix(cells_2_cut2,nrow = 1)
  write.table(cells_2_cut2,"Cells_2_remove.csv",row.names = FALSE,col.names = FALSE)
  writeLines("cell's indices to remove are saved in 'cells 2 remove.csv',\n
             please rerun MAIN script again \n")
  stopQuietly()
}

writeLines(' Press 1 if you want to perform additional analysis (e.g. lineage inference, cell ordering) , 0 otherwise: ')
Proceed=readLines(n=1)
if(as.integer(Proceed)!=1){
  stopQuietly()
}
  
  # %% *** 3-RECONSTRUCTION OF LINEAGE PROGRESSION ***
  # %
  # % Please check comments in 'CALISTA_transition_main' for more information. 
  Results=CALISTA_transition_main(DATA,Results)
  
  # %% *** 4-DETERMINATION OF TRANSITION GENES ***
  # %
  # % Please check comments in 'CALISTA_transition_genes_main' for more information.
  Results=CALISTA_transition_genes_main(DATA,INPUTS,Results)
  
  # %% *** 5-PSEUDOTEMPORAL ORDERING OF CELLS ***
  # %
  # % Please check comments in 'CALISTA_ordering_main' for more information. 
  # 
  Results=CALISTA_ordering_main(DATA,Results)
  
  # %% *** 6-PATH ANALYSIS ***
  # % Please check comments in 'CALISTA_path_main' for more information.   
  writeLines('Press 1 if you want to select transition paths and perform additional analysis, 0 otherwise:')
  proceed=as.integer(readLines(n=1))
  if(!proceed){
    stopQuietly()
  }
  Results=CALISTA_path_main(INPUTS,Results)
