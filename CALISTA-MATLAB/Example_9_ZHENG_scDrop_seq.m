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
% Specify data types and settings for pre-processing
% Specify data types and settings for pre-processing
INPUTS.data_type=5; % snRNA-seq data
INPUTS.format_data=6; % (RECOMMENDED FOR Droplet-based datasets) File .mat containing the following variables:
%                       'expMatrix': NxG expression matrix with rows=cells and columns=genes;
%                       'cellID': cell array/character vector containing cell IDs (not required)
%                       'geneNAMES': cell array containing gene names (not required)
%                       'timeline': 1xN numerical vector containing the time information (not required)
INPUTS.zeros_genes=3; % Remove genes with > 3 zero values
INPUTS.zeros_cells=200; % Remove cells with zero values>200
INPUTS.top_genes=200; % Retain only top 200 most variable genes
%
% Specify single-cell clustering settings
INPUTS.use_drop_prob_in_clustering=1;
INPUTS.optimize=0; % The number of cluster is known a priori
% INPUTS.plot_tsne=0;
% Upload and pre-process data 
[DATA,INPUTS] = import_data(INPUTS);

%% *** 2-SINGLE-CELL CLUSTERING ***
%
% Type 'help CALISTA_clustering_main' for more information. 
[Results, DATA, INPUTS]=CALISTA_clustering_main(DATA,INPUTS);

% cgo = clustergram(W_1','Colormap',redbluecmap,'Standardize','row','Linkage','average')
% Testing
figure
hist(Results.final_groups,1:Results.expected_clusters)
title('Cells in each cluster')
xlabel('Cluster')
ylabel('Num of cells')

% % % [p]=numSubplots(Results.expected_clusters);
% % % figure
% % % for i=1:Results.expected_clusters
% % %     subplot(p(1),p(2),i)
% % %     hist(Results.cell_labels(find(Results.final_groups==i)),1:length(unique(Results.cell_labels)));
% % % end
% % % % title('Original assignments in each cluster from CALISTA')
% % % 
% % % [p]=numSubplots(length(unique(Results.cell_labels)));
% % % figure
% % % for i=1:length(unique(Results.cell_labels))
% % %     subplot(p(1),p(2),i)
% % %     hist(Results.final_groups(find(Results.cell_labels==i)),1:Results.expected_clusters);
% % % end
% % % 
% % % %% Construct graph for Sankey diagram
% % % true_labs=Results.cell_labels;
% % % calista_labs=Results.final_groups';
% % % count=1;
% % % links=[];
% % % for i=1: Results.expected_clusters
% % %     for j=1: Results.expected_clusters
% % %         aa=find(true_labs==i);
% % %         bb=find(calista_labs==j);
% % %         cc=intersect(aa,bb);
% % %         if ~isempty(cc)    
% % %         links{count,1} = sprintf( '%s %2i %s %2i %s %2i', 'Zheng K', i,' [',length(cc),'] Calista K',j );
% % %         count=count+1;
% % %         end       
% % %         
% % %     end
% % % end

% % % %% ARI population and final_groups
% % % 
% % % population=Results.clustering_struct.all.all.population;
% % % count=1;
% % % for i=1: size(population,1)
% % %     [AR(i)]=RandIndex(population(i,:),Results.final_groups);
% % %     [AR2(i)]=RandIndex(population(i,:),Results.cell_labels);
% % %     labels{count}='Pairwise';
% % %     count=count+1;
% % % end
% % % for i=1: size(population,1)
% % %     
% % %     [AR2(i)]=RandIndex(population(i,:),Results.cell_labels);
% % %      labels{count}='Reference';
% % %      count=count+1;
% % % end
% % % 
% % % to_remove=find(AR==1,1,'first');
% % % AR(to_remove)=[];
% % % labels(to_remove)=[];
% % % figure
% % % boxplot([AR'; AR2'],labels)
% % % ylabel('ARI')
