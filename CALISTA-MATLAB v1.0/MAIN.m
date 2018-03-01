%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             MAIN SCRIPT                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

% Add path for CALISTA folder
mfilename=matlab.desktop.editor.getActiveFilename;
curr_path=fileparts(which(mfilename));
cd(curr_path);

addpath(genpath(curr_path))
warning ('off','all');

%% *** 1-Data Import and Preprocessing ***
%
% Type 'help import_data' for more information 
%
% Specify data types and settings for pre-processing
INPUTS.data_type=1; % Single-cell RT-qPCR CT data
INPUTS.format_data=3; % Rows= cells and Columns= genes with time/stage info in the last column
INPUTS.data_selection=[]; % Include data from all time points
INPUTS.perczeros_genes=100; % Remove genes with > 100% of zeros 
INPUTS.perczeros_cells=100; % Remove cells with 100% of zeros
INPUTS.cells_2_cut=0; % No manual removal of cells
INPUTS.perc_top_genes=100; % Retain only top X the most variable genes with X=min(200, INPUTS.perc_top_genes * num of cells/100, num of genes)


% Specify single-cell clustering settings
INPUTS.optimize=0; % The number of cluster is known a priori
INPUTS.parallel=1; % Use parallelization option
INPUTS.runs=50; % Perform 50 independent runs of greedy algorithm 
INPUTS.max_iter=100; % Limit the number of iterations in greedy algorithm to 100
INPUTS.cluster='kmedoids'; % Use k-medoids in consensus clustering

% Specify transition genes settings
INPUTS.thr_transition_genes=50; % Set threshold for transition genes determination to 50%

% Specify path analysis settings
INPUTS.plot_fig=1; % Plot figure of smoothened gene expression along path
INPUTS.hclustering=1; % Perform hierarchical clustering of gene expression for each path 
INPUTS.method=2; % Use pairwise correlation for the gene co-expression network (value_cutoff=0.8, pvalue_cutoff=0.01)
INPUTS.moving_average_window=10; % Set the size of window (percent of cells in each path) used for the moving averaging
% 
% Upload and pre-process data 
DATA = import_data(INPUTS);

%% *** 2-SINGLE-CELL CLUSTERING ***
%
% Type 'help CALISTA_clustering_main' for more information. 
[Results, DATA, INPUTS]=CALISTA_clustering_main(DATA,INPUTS);

%  Cluster removal (if desired)

cluster_cut=input('Press 1 if you want to remove cell cluster(s). Press 0 otherwise: ');

if cluster_cut==1
    
   which_cut=input('Enter cluster index for removal (e.g 1 or [5 3]): ');
   cells_2_cut2=find(ismember(Results.final_groups,which_cut)==1);
   csvwrite('Cells_2_remove',cells_2_cut2)   
   RESULTS.cluster_cut = 1;
   fprintf('Cell''s indices to remove are saved in "Cells_2_remove.csv". Please run MAIN script again \n')
   return
end

Proceed=input(' Press 1 if you want to perform additional analysis (e.g. lineage inference, cell ordering) , 0 otherwise: ');
if ~Proceed
    return
end

%% *** 3-RECONSTRUCTION OF LINEAGE PROGRESSION ***
%
% Type 'help CALISTA_transition_main' for more information. 

[Results]=CALISTA_transition_main(DATA,Results);

%% *** 4-DETERMINATION OF TRANSITION GENES ***
%
% Type 'help CALISTA_transition_genes_main' for more information. 

[Results]=CALISTA_transition_genes_main(DATA,INPUTS,Results);

%% *** 5-PSEUDOTEMPORAL ORDERING OF CELLS ***
%
% Type 'help CALISTA_ordering_main' for more information. 

[Results]=CALISTA_ordering_main(DATA,Results);

%% *** 6-PATH ANALYSIS ***
%
% Type 'help CALISTA_path_main' for more information. 

Proceed2=input(' Press 1 if you want to select transition paths and perform additional analysis, 0 otherwise: ');

if ~Proceed2
    return
end
Results=CALISTA_path_main(INPUTS,Results);