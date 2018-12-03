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
% Specify landscape plotting settings
INPUTS.ngrid=40; % Grid size for 2D interpolation
%
% Specify path analysis settings
INPUTS.plot_fig=1; % Plot figure of smoothened gene expression along path
INPUTS.hclustering=1; % Perform hierarchical clustering of gene expression for each path 
INPUTS.method=2; % Use pairwise correlation for the gene co-expression network (value_cutoff=0.8, pvalue_cutoff=0.01)
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


