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

%% *** 0-Data Import and Preprocessing ***
%
% Type 'help import_data' for more information 
%
INPUTS.data_type=3; % scRNA-seq data  (Expression values - e.g log(TPM+1) or log(RPKM+1))
INPUTS.format_data=5; % Manual selection from data table
INPUTS.zeros_genes=0.9;%0.9;%3; % Remove genes with > 100% of zeros 
INPUTS.top_genes=41;% Retain only top 41 genes (10% of num of cells)
INPUTS.plot=0; % Do not display additional plots
% Specify single-cell clustering settings
INPUTS.optimize=1; % Find the optimal number of clusters
%
% Upload and pre-process data 
[DATA,INPUTS] = import_data(INPUTS);

%% *** 2-SINGLE-CELL CLUSTERING ***
%
% Type 'help CALISTA_clustering_main' for more information. 
[Results, DATA, INPUTS]=CALISTA_clustering_main(DATA,INPUTS);

%  Cluster removal (if desired)
cluster_cut=input('Press 1 if you want to remove cell cluster(s). Press 0 otherwise: ');

if cluster_cut==1
    
   which_cut=input('Enter cluster index for removal (e.g 1 or [5 3]): ');
   cells_2_cut2=find(ismember(Results.final_groups,which_cut)==1);
   cells_2_cut2=DATA.cut_sort.abs_indices(cells_2_cut2);
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

[Results]=CALISTA_transition_main(DATA,INPUTS,Results);

%% *** 4-DETERMINATION OF TRANSITION GENES ***
%
% Type 'help CALISTA_transition_genes_main' for more information. 

[Results]=CALISTA_transition_genes_main(DATA,INPUTS,Results);

%% *** 5-PSEUDOTEMPORAL ORDERING OF CELLS ***
%
% Type 'help CALISTA_ordering_main' for more information. 

[Results]=CALISTA_ordering_main_2(DATA,INPUTS,Results);

%% *** 6-LANDSCAPE PLOTTING ***
%
% Type 'help CALISTA_landscape_plotting_main' for more information.

CALISTA_landscape_plotting_main(INPUTS,Results);
%% *** 7-PATH ANALYSIS ***
%
% Type 'help CALISTA_path_main' for more information. 

Proceed2=input(' Press 1 if you want to select transition paths and perform additional analysis, 0 otherwise: ');

if ~Proceed2
    return
end

Results=CALISTA_path_main(DATA,INPUTS,Results);


