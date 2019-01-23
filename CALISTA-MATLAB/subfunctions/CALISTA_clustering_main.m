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
% ** INPUTS.cluster_time **
% 1- cluster data from each time point separately
% 0- otherwise (by default)
%
% ** INPUTS.optimize **
% 1- select the number of clusters by eigengap plot
% 2- select the number of clusters by BIC
% 0- define the number of clusters
% By default:
% if INPUTS.cluster_time==0 (see import_data.m) -> INPUTS.optimize=0;
% otherwise INPUTS.optimize=1;
%
% ** INPUTS.parallel **
% 1- Use parallel processing (number of cores available - 1) (by default)
% 0- Do not use processing
%
% ** INPUTS.max_iter **
% Maximum number of iterations in the greedy algorithm
% INPUTS.max_iter = 200; (by default)
%
% ** INPUTS.runs **
% Number of clustering runs
% If INPUTS.data_type<5 -> INPUTS.runs=50; 
% Otherwise INPUTS.runs=20; 
%
% ** INPUTS.cluster **
% 'hierarchical'- hierachical clustering of consensus matrix
% 'kmedoids'-  kmedoids for the clustering of consensus matrix (by default)
%
% ** INPUTS.plot_tsne **
% 1- for t-sne plot after clustering (default) 
% 0- otherwise 
% For more info visit https://github.com/KlugerLab/FIt-SNE and http://www.fftw.org
%
% ** INPUTS.tsne_opts.perplexity **
% Perplexity parameter for tnse (usually a value between 1 and 30 (default))
% 
% ** INPUTS.tsne_opts.input_data 
% 0- Use the expression data as input
% otherwise a values between 1 and number of genes to use top 'INPUTS.tsne_opts.input_data'
% principal components as input
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
% Copyright. October 1, 2018.
%% CHECK INPUT ARGUMENTS

if nargin <2
    error('Not enough input arguments')
end


if ~isfield(INPUTS,'cluster_time')
     INPUTS.cluster_time=0; 
end

if ~isfield(INPUTS,'optimize')
    if INPUTS.cluster_time==0
        INPUTS.optimize=0;
    else
        INPUTS.optimize=1;
    end
end

if ~isfield(INPUTS,'parallel')
     INPUTS.parallel=1; 
end

if ~isfield(INPUTS,'runs')
    if INPUTS.data_type<5
        INPUTS.runs=50; 
    else
        INPUTS.runs=20; 
    end
     
end

if ~isfield(INPUTS,'max_iter')
     INPUTS.max_iter=200;
end

if ~isfield(INPUTS,'algorithm')
     INPUTS.algorithm='greedy_cabsel';
end

if ~isfield(INPUTS,'cluster')
     INPUTS.cluster='kmedoids';
end

if ~isfield(INPUTS,'use_drop_prob_in_clustering')
    if INPUTS.data_type<5
        INPUTS.use_drop_prob_in_clustering=0;
    else
        INPUTS.use_drop_prob_in_clustering=1;
    end
end

if ~isfield(INPUTS,'plot_tsne')
    INPUTS.plot_tsne=0;
end

if  INPUTS.plot_tsne==1
    if ~isfield(INPUTS,'opts')
    INPUTS.tsne_opts.perplexity = 30;
    INPUTS.tsne_opts.input_data = 0;
    else
        if ~isfield(INPUTS.tsne_opts,'perplexity')
           INPUTS.tsne_opts.perplexity=30;
        end
         if ~isfield(INPUTS.tsne_opts,'input_data')
           INPUTS.tsne_opts.input_data=0;
         end
    end
end

tic
if INPUTS.cluster_time==0
    [Results, DATA, INPUTS]=CALISTA_clustering_classic(DATA,INPUTS);
else
    [Results, DATA, INPUTS]=CALISTA_clustering_time(DATA,INPUTS);
end
toc
fprintf('\nDONE.\n')



