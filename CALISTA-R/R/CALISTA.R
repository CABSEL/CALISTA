CALISTA<-function(INPUTS){
  
  cur_path=getSrcDirectory(function(x) {x})
  
  ###parameters loading
  Parameters=readMat('./Two-state model parameters/Parameters.mat')
  
  # %% *** 1-INITIALIZATION ***
  #   
  #   % INPUTS:
  #   % 
  # % ** data_type **
  #   % 1- for scRT-qPCR  (CT values)  (DEFAULT) 
  # % 2- for scRT-qPCR  (Expression values - Num. mRNA molecules)
  # % 3- for scRNA-seq  (Expression values - e.g log(TPM+1) or log(RPKM+1))
  # % 4- for scRNA-seq  (Expression values - e.g TPM+1 or RPKM)
  # %
  # %
  # % ** format_data **%
  #   % Five data formats are accepted:
  #   %
  # % 1- Rows= cells and Columns= genes with time info (last column of doubles) (DEFAULT)  
  # % ------------------------------------------
  #   % Gene1  Gene2  Gene3  ... Genej  Time/Stage
  # %  27     80     56    ...  69        0
  # %  73     20     90    ...  45        0
  # %   .      .      .    ...   .        .
  # %   .      .      .    ...   .        .
  # %   .      .      .    ...   .        .
  # % ------------------------------------------
  #   % 
  # % 
  # %
  # % 2- Rows= genes and Columns= cells with time info (first row of doubles) 
  # % ----------------------------------------------
  #   % Time/Stages   0      0      1   ...  16     16
  # %   Gene1      27     80     56   ...  69      0
  # %   Gene2      73     20     90   ...  45      0
  # %   Gene3       .      .      .   ...   .      .
  # %   Gene4       .      .      .   ...   .      .
  # %   Genej       .      .      .   ...   .      .
  # % ----------------------------------------------
  #   % 
  # % 3- Rows= cells and Columns= genes without time info
  # % ------------------------------
  #   % Gene1  Gene2  Gene3  ... Genej  
  # %  27     80     56    ...  69       
  # %  73     20     90    ...  45       
  # %   .      .      .    ...   .       
  # %   .      .      .    ...   .       
  # %   .      .      .    ...   .       
  # % ------------------------------
  #   % 
  # % 4- Rows= genes and Columns= cells without time info
  # % --------------------------------------------
  #   % Gene1      27     80     56   ...  69      0
  # % Gene2      73     20     90   ...  45      0
  # % Gene3       .      .      .   ...   .      .
  # % Gene4       .      .      .   ...   .      .
  # % Genej       .      .      .   ...   .      .
  # % --------------------------------------------
  #   % 
  # % 
  # %
  # % 5- Manual data selection (Please see README file)
  # %
  # %
  # % ** data_selection **
  #   % CALISTA groups single cells in “snapshots” based on the their time 
  # % or cell stage information (whenever available). 
  # % Hence, the user can select specific time snapshots for further analysis. 
  # % For instance, considering a datset with cells whose expresisons 
  # % were measured at 4 time points ( 0, 24, 48 and 96 h), 
  # % snapshots can be selected as follows:
  #   % data_selection =[ ] or (1:4)  for all time points/cells (DEFAULT)
  # % data_selection =1               for cells at 0 h
  # % data_selection =[1 3 4]       for cells at 0, 48 and 96 h
  # % data_selection =(2:3)          for cells at 24 and 48 h
  # % 
  # % 
  # % ** perczeros_genes **
  #   % Genes with more than perczeros_genes% of zero values (between 0 and 100) are
  # % removed. perczeros_genes=90 (DEFAULT)
  # % 
  # % 
  # % ** perczeros_cells **
  #   % Cells and cells with more than perczeros_cells% of zero values (between 0 and 100) are
  # % removed. perczeros_cells=100 (DEFAULT)
  # % 
  # % 
  # % ** cells_2_cut **
  #   % 1- Remove cells based on their indices on the expression matrix. Indices need to be uploaded as a csv file. (see README file)
  # % 0- No cells delection (DEFAULT)
  # % 
  # % 
  # % ** perc_top_genes **
  #   % This input is only required for RNA-seq datasets.
  # % The number of most variable genes selected for further analysis is min(200, perc_top_genes% of total num of genes)
  # % perc_top_genes=50 (DEFAULT)
  
  
  ### 1.a: Data loading
  
  
  DATA=import_data(INPUTS)
  
  # %% *** 2-CLUSTERING ***
  #   
  #   % INPUTS:
  #   %
  # % ** optimize **
  #   % 1- define the optimal number of clusters based on eigengap values
  # % 0- define the optimal number of clusters a priori (DEFAULT)
  # % 
  # % 
  # % ** algorithm **
  #   % 'greedy_cabsel' % swap cells based on maximum likelihood (DEFAULT)
  # % 'greedy_maxdiff' % swap cells based on maximum difference between clusters    
  # % 'cabsel_sabec' % swap cells based on maximum difference and simulated annealing (SABEC)    
  # % 
  # % 
  # % ** cluster **
  #   % 'hierarchical' % use hierachical clustering (DEFAULT)
  # % 'kmedoids' % use kmedoids for the clustering of consensus matrix
  # % 
  # % 
  # % ** parallel **
  #   % 1- Use parallel processing. The number of CPUs used for parallel execution minus 1 for system stability. (DEFAULT)
  # % 0- Do not use processing
  
  
  Results=list()
  CALISTA_clustering_main_results=CALISTA_clustering_main(DATA,Parameters,INPUTS)
  Results=CALISTA_clustering_main_results$Results
  DATA=CALISTA_clustering_main_results$DATA
  INPUTS=CALISTA_clustering_main_results$INPUTS
  if(Results$go==1 && Results$cluster_cut!=1){
    ###3-Lineage inference
    CALISTA_transition_main_results=CALISTA_transition_main(Results ,DATA)
    Results=CALISTA_transition_main_results$Results
    DATA=CALISTA_transition_main_results$DATA
    ###4-Cell ordering
    Results=CALISTA_ordering_main(Results,Parameters,DATA)
    ###5-Transition genes
    INPUTS$thr_transition_genes=75
    Results=CALISTA_transition_genes_main(Results,DATA,Parameters,INPUTS)
    # %% *** 6-PATH ANALYSIS ***
    #   
    #   % INPUTS:
    #   % 
    # % 
    # % ** plot_fig **
    #   % 1- plot figure (DEFAULT)
    # % 0- otherwise
    # %
    # %
    # % ** hclustering **
    #   % 1- run hierarchical clustering for each path (DEFAULT)
    # % 0- otherwise
    # % 
    # % 
    # % ** method **
    #   % 1- for partial correlation (value_cutoff=0.4, pvalue_cutoff=0.05)
    # % 2- for correlation (value_cutoff=0.8, pvalue_cutoff=0.01) (DEAFULT)
    # % 
    # % 
    
    
    Results=CALISTA_path_main(Results,INPUTS)
  }
}
