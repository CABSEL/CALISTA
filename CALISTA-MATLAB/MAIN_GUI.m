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
INPUTS.data_type=1; % Single-cell RT-qPCR CT data
INPUTS.format_data=1; % Rows= cells and Columns= genes with time/stage info in the last column
INPUTS.data_selection=[]; % Include data from all time points
INPUTS.zeros_genes=1; % Remove genes with 100% of zeros 
INPUTS.zeros_cells=1; % Remove cells with 100% of zeros
INPUTS.cut_variable_genes=1000; % Select 1000 most variable genes
INPUTS.top_genes=1; % Retain only top Y the top 100% genes based on CALISTA gene selection step (all genes)
INPUTS.cells_2_cut=0; % No manual removal of cells
INPUTS.cluster_time=0; % Run CALISTA clustering on all datasets
INPUTS.plot=0; %Do not show additional plots
%
% Specify single-cell clustering settings
INPUTS.optimize=1; % Find the optimal number of clusters
INPUTS.parallel=1; % Use parallelization option
INPUTS.runs=50; % Perform 50 independent runs of greedy algorithm 
INPUTS.max_iter=200; % Limit the number of iterations in greedy algorithm to 200
INPUTS.cluster='kmedoids'; % Use k-medoids in consensus clustering
INPUTS.use_drop_prob_in_clustering=0; % Do not consider the dropout probability in CALISTA
%
% Specify lineage inference settings
INPUTS.transition_new=0; % Users can manually add and prune trandistion edges
%
% Specify transition genes settings
INPUTS.thr_transition_genes=50; % Set threshold for transition genes determination to 50%
INPUTS.plot_trans_genes=0; % Do not show additional plots

% Specify landscape plotting settings
INPUTS.ngrid=40; % Grid size for 2D interpolation
INPUTS.gridshift=3; % Shift for cell representation on the landscape
%
% Specify path analysis settings
INPUTS.plot_fig=1; % Plot figure of smoothened gene expression along path
INPUTS.hclustering=1; % Perform hierarchical clustering of gene expression for each path 
INPUTS.method=2; % Use pairwise correlation for the gene co-expression network (value_cutoff=0.8, pvalue_cutoff=0.01)
INPUTS.moving_average_window=10; % Set the size of window (percent of cells in each path) used for the moving averaging
INPUTS.path_auto=0; % Manually define paths for further analysis

%% Analysis with CALISTA:
% -SINGLE-CELL CLUSTERING
% -RECONSTRUCTION OF LINEAGE PROGRESSION
% -DETERMINATION OF TRANSITION GENES
% -PSEUDOTEMPORAL ORDERING OF CELLS
% -LANDSCAPE PLOTTING
% -PATH ANALYSIS

CALISTA_LI_GUI;
