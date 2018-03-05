function [Results, DATA, INPUTS]=CALISTA_clustering_main(DATA,INPUTS)
%CALISTA_CLUSTERING_MAIN single-cell clustering in CALISTA
% CALISTA implements two-step clustering algorithm. The first step involves
% a likelihood-based clustering based on the stochastic two-state model of
% gene transcriptional process. A greedy algorithm is implemented to obtain
% cell clustering with maximizes the overall likelihood, and repeatedly run
% to produce a consensus matrix (i.e. the matrix of the number of times two
% cells are clustered together). The second and final step involves
% k-medoids or hierachical clustering using the consensus matrix. 
%
% Usage:
% [Results, DATA, INPUTS]=CALISTA_clustering_main(DATA,INPUTS)
% Perform single-cell clustering using user-defined input values
% 
% Inputs:
% INPUTS - a structure containing settings for the single-cell clustering 
%
% ** INPUTS.optimize **
% 1- select the number of clusters by eigengap plot
% 0- define the number of clusters 
%
% ** INPUTS.parallel **
% 1- Use parallel processing (number of cores available - 1)  
% 0- Do not use processing
%
% ** INPUTS.max_iter **
% Maximum number of iterations in the greedy algorithm 
% INPUTS.max_iter = 100; 
% 
% ** INPUTS.runs **
% Number of clustering runs 
% INPUTS.runs = 50; 
%
% ** INPUTS.cluster **
% 'hierarchical'- hierachical clustering of consensus matrix
% 'kmedoids'-  kmedoids for the clustering of consensus matrix
% 
% DATA - a structure containing preprocessed single-cell expression data
% Use 'import_data' to upload and preprocess single-cell expression values
%
% Outputs:
% Results - a structure containing the results of CALISTA analysis.
% The most relevant fields containing the clustering analysis are:
%
% Results.final_groups - vector containing cluster assignment for the cells
%
% Results.singleCELLclusterDATA - a 1xK cell array where K is the number of 
% clusters. Each cell contains a GxNk (normalised) expression matrix of G 
% genes and Nk cells in the k-th cluster.
%
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. June 1, 2017.

if nargin <2 
    error('Not enough input arguments')
end
INPUTS.algorithm='greedy_cabsel';  
Parameters=DATA.Parameters;
optimize=INPUTS.optimize;
algorithm=INPUTS.algorithm;  
Cluster=INPUTS.cluster;
parallel=INPUTS.parallel;
loops=INPUTS.runs;
max_iter=INPUTS.max_iter;

if optimize
    max_clust=12;
    [ my_results] = CALISTA_clustering(DATA.totDATA,Parameters.logP,Parameters.sets,max_clust,'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops);
    
    % Eigengap heuristic
    consensus=get_consensus(my_results.all.all.population);
    [Results.expected_clusters]=eigengap(consensus,max_clust);
    
else
    Results.expected_clusters=input(' Number of clusters expected: ');
end

% *** 2.b-Cell clustering ***
[ my_results_c] = CALISTA_clustering(DATA.totDATA,Parameters.logP,Parameters.sets,Results.expected_clusters,'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops);
Results.final_groups=my_results_c.all.all.idx';

% *** 2.c-Relabelling based on time/cell stage info ***
method=3;
[Results,DATA]= find_progression2(Results,DATA,method);

% *** 2.c-Optimal parameter estimation for the final cluster assignment ***
[ my_results_final] = CALISTA_clustering(DATA.totDATA,Parameters.logP,Parameters.sets,Results.expected_clusters,'parallel',parallel,'get_k',Results.final_groups,'algorithm',algorithm);

% *** 2.e-Cluster visualization ***
reduction=2;
[Results]=visualization(reduction,DATA,Results);

Results.clustering_struct=my_results_final;

