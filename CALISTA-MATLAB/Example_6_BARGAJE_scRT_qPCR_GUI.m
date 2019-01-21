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
% Specify pre-processing settings
INPUTS.data_type=1; % Single-cell RT-qPCR CT data
INPUTS.format_data=1; % Rows= cells and Columns= genes with time/stage info in the last column
%
% Specify single-cell clustering settings
INPUTS.optimize=1; % Find the optimal number of clusters
%
% Specify landscape plotting settings
INPUTS.ngrid=40; % Grid size for 2D interpolation
%
% Specify path analysis settings
INPUTS.plot_fig=1; % Plot figure of smoothened gene expression along path
INPUTS.hclustering=1; % Perform hierarchical clustering of gene expression for each path 
INPUTS.method=2; % Use pairwise correlation for the gene co-expression network (value_cutoff=0.8, pvalue_cutoff=0.01)
INPUTS.path_auto=0; % Manually define paths for further analysis

%% Analysis with CALISTA:
% -UPLOAD AND PRE_PROCESS DATA
% -SINGLE-CELL CLUSTERING
% -RECONSTRUCTION OF LINEAGE PROGRESSION
% -DETERMINATION OF TRANSITION GENES
% -PSEUDOTEMPORAL ORDERING OF CELLS
% -LANDSCAPE PLOTTING
% -PATH ANALYSIS

CALISTA_LI_GUI;
