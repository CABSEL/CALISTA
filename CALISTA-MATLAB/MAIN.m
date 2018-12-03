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
INPUTS.data_selection=[];%[4:7]; % Include data from all time points
INPUTS.zeros_genes=1; % Remove genes with 100% of zeros 
INPUTS.zeros_cells=1; % Remove cells with 100% of zeros
INPUTS.cut_variable_genes=1000; % Select 1000 most variable genes
INPUTS.top_genes=1; % Retain only top Y the top 100% genes based on CALISTA gene selection step (all genes)
INPUTS.cells_2_cut=0; % No manual removal of cells
INPUTS.cluster_time=0; % Run CALISTA clustering on all datasets
INPUTS.plot=0; %Do not show additional plots
%
% Specify single-cell clustering settings
INPUTS.optimize=0; % Find the optimal number of clusters
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
% [Results_transition_genes]=CALISTA_ordering_main_2_transition_genes(DATA,INPUTS,Results);
%% *** 6-LANDSCAPE PLOTTING ***
%
% Type 'help CALISTA_landscape_plotting_main' for more information.

CALISTA_landscape_plotting_main(INPUTS,Results);
% CALISTA_landscape_plotting_main(INPUTS,Results_transition_genes);
% for k=1:Results_transition_genes.expected_clusters
%     Results_transition_genes.singleCELLclusterDATA{k}=Results_transition_genes.singleCELLclusterDATA{k}(:,Results_transition_genes.GENES.idx_tot_transition_genes);
%     
% end
%% *** 7-PATH ANALYSIS ***
%
% Type 'help CALISTA_path_main' for more information. 

Proceed2=input(' Press 1 if you want to select transition paths and perform additional analysis, 0 otherwise: ');

if ~Proceed2
    return
end

Results=CALISTA_path_main(DATA,INPUTS,Results);


