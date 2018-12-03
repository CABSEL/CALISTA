function [Results]=CALISTA_transition_main(DATA,INPUTS,Results,cell_assignments)
%CALISTA_TRANSITION_MAIN infer lineage progression among cell clusters
% CALISTA uses the cluster distances - a measure of cluster-cluster
% dissimilarity - to infer the lineage progression or cluster-cluster
% relationship. 
% CALISTA also provides a simple user-interface to add and remove edges
% based on the cluster distances. 
% 
% Usage:
% Run CALISTA lineage inference using CALISTA single-clustering result
% [Results]=CALISTA_transition_main(DATA,Results)
% 
% Run CALISTA lineage Inference with user-defined cell clusters 
% [Results]=CALISTA_transition_main(DATA,[],cell_assignments)
%    
% Inputs:   
% DATA - a structure containing preprocessed single-cell expression data
% Use 'import_data' to upload and preprocess single-cell expression values
%
% In addition to the specification in 'import_data', users need to specify:
% ** INPUTS.transition_new **
% 1- CALISTA automatically detect transition edges based on time/stage info(required for num of clusters >15)
% 0- The lineage progression graph is built based on adding
% edges between clusters in increasing magnitude of cluster distance.
% By default:
% If INPUTS.data_type>=5 or Results.expected_clusters>15 -> INPUTS.transition_new=1;
% Otherwise INPUTS.transition_new=0;
%
% Results - a structure of CALISTA clustering results
% Run 'CALISTA_clustering_main'
%
% cell_assignments - 1xN vector of INTEGERS with N = number of cells. 
% The n-th element of cell_assignments contains the cluster assignment of
% the n-th cell of the expression data uploaded. Cluster names must
% be assigned in sequence (e.g. 1,2,3,4 and not 1,2,4).
% 
%
% Outputs:
% Results - a structure containing the results of CALISTA analysis. 
% The most relevant field containing the clustering analysis is:
%
% Results.TRANSITION.nodes_connection - an Ex2 matrix containing E edges 
% of the lineage progression graph. The first and second columns give the
% clusters incident to the edges. If time/stage information is available,
% then the edges are directed where the first column indicates the source
% cluster and the second column indicates the target cluster. The user can
% also specify the starting cell or marker gene to determine the direction
% of the edges. 
%
% Results.TRANSITION.cluster_distance - a ranked list of cluster distances. 
% The first and second columns are cluster indices and the thrid column 
% gives the cluster distance.
%
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright.  October 1, 2018.


% Error check
if nargin <3
    error('Not enough input arguments')
end 
 
if isempty(Results)
   if nargin <4
    error('Not enough input arguments')
   end 
   if ~isvector(cell_assignments)
       error('Please upload the cell assignments as a vector')
   end
   Results=jump_clustering(DATA,INPUTS,cell_assignments);
end

if ~isfield(INPUTS,'transition_new')
    INPUTS.transition_new=0;
end

if Results.expected_clusters>15
    if length(unique(DATA.timeline))==1 % No time info or only one time
         if INPUTS.transition_new==0
            fprintf('\nNum of clusters>15. Please await or re-rum CALISTA with time info and INPUTS.transition_new=1\n')
         end
    else
        if INPUTS.transition_new==0
            fprintf('\nNum of clusters>15. CALISTA automatically detect transition edges based on time/stage info\n')
            INPUTS.transition_new=1;
        end
    end
end

if INPUTS.transition_new==1
    [Results]=CALISTA_transition_new(DATA,Results);
else
    [Results]=CALISTA_transition(DATA,Results);
    % *** 3.b-Plot mean expressions for each cluster ***
    Results=plot_cluster_mean_exp(DATA,Results);
    
    % *** 3.c-Cell-Cell variability analysis ***
    [Results]=cell_variability(DATA,Results);
end

fprintf('\nDONE.\n')
