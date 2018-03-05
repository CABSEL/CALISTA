%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           MAIN_GUI SCRIPT                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
closereq

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
INPUTS.format_data=1; % Rows= cells and Columns= genes with time/stage info in the last column
INPUTS.data_selection=[]; % Include data from all time points
INPUTS.perczeros_genes=100; % Remove genes with > 100% of zeros 
INPUTS.perczeros_cells=100; % Remove cells with 100% of zeros
INPUTS.cells_2_cut=0; % No manual removal of cells
INPUTS.perc_top_genes=100; % Retain only top X the most variable genes with X=min(200, INPUTS.perc_top_genes * num of cells/100, num of genes)


% Specify single-cell clustering settings
INPUTS.optimize=1; % The number of cluster is known a priori
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

%% Analysis with CALISTA:
% -SINGLE-CELL CLUSTERING
% -RECONSTRUCTION OF LINEAGE PROGRESSION
% -DETERMINATION OF TRANSITION GENES
% -PSEUDOTEMPORAL ORDERING OF CELLS
% -PATH ANALYSIS

CALISTA_LI_GUI;
