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
INPUTS.data_type=5; % snRNA-seq data
INPUTS.format_data=7; % (RECOMMENDED FOR Droplet-based datasets) File .mat containing the following variables:
%                       'expMatrix': NxG Log-expression matrix with rows=cells and columns=genes;
%                       'cellID': cell array/character vector containing cell IDs (not required)
%                       'geneNAMES': cell array containing gene names (not required)
%                       'timeline': 1xN numerical vector containing the time information (not required)INPUTS.zeros_genes=0.9; % Remove genes with > 90% of zeros 
INPUTS.zeros_genes=0.9; % Remove genes with > 90% of zeros 
INPUTS.top_genes=100; % Retain only top 100 genes
INPUTS.cluster_time=1; % Run CALISTA clustering for each time point separately
INPUTS.plot=0; % Do not show additional plots
% 
% Specify single-cell clustering settings
INPUTS.optimize=1; % Find the optimal number of clusters
INPUTS.use_drop_prob_in_clustering=1; % Consider the dropout probability in CALISTA
%
% Specify lineage inference settings
INPUTS.transition_new=1; % CALISTA automatically detect transition edges(required for Droplet-based datasets)
%
% Specify transition genes settings
INPUTS.thr_transition_genes=50; % Set threshold for transition genes determination to 50%
INPUTS.plot_trans_genes=0; % Do not show additional plots
%
% Specify landscape plotting settings
INPUTS.ngrid=300; % Grid size for 2D interpolation
%
% Specify path analysis settings
INPUTS.moving_average_window=10; % Set the size of window (percent of cells in each path) used for the moving averaging
INPUTS.path_auto=1; %Automatically detect paths along the graph based on shortestpath



% Upload and pre-process data 
[DATA,INPUTS]= import_data(INPUTS);

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
% % Save 3D video
% figure(100)
% % daspect([2,2,0.3]);
% OptionZ.FrameRate=20;OptionZ.Duration=20.5;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'CALISTA_lineage_3Dvideo',OptionZ)
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

% CALISTA_landscape_plotting_main(INPUTS,Results);

%% *** 7-PATH ANALYSIS ***
%
% Type 'help CALISTA_path_main' for more information. 

Proceed2=input(' Press 1 if you want to select transition paths and perform additional analysis, 0 otherwise: ');

if ~Proceed2
    return
end

[FileName,PathName,FilterIndex] = uigetfile('*.*');
filename=strcat(PathName, FileName);
selected_genes=importdata(filename);

Results=CALISTA_path_main(DATA,INPUTS,Results,'selected_genes',selected_genes);

% % check final cell types
% [FileName,PathName,FilterIndex] = uigetfile('*.*');
% filename=strcat(PathName, FileName);
% metadata=importdata(filename);
% metadata(DATA.cut_sort.idx2cutCELL,:)=[];
% metadata=metadata(DATA.cut_sort.idx_sorted_cells,:);
% 
% [FileName,PathName,FilterIndex] = uigetfile('*.*');
% filename=strcat(PathName, FileName);
% lineage=importdata(filename);
% 
% metadata2=metadata(find(DATA.timeline==max(DATA.timeline)),:);
% 