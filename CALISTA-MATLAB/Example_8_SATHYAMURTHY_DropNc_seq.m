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
INPUTS.data_type=5; % snRNA-seq data
INPUTS.format_data=6; % (RECOMMENDED FOR Droplet-based datasets) File .mat containing the following variables:
%                       'expMatrix': NxG expression matrix with rows=cells and columns=genes;
%                       'cellID': cell array/character vector containing cell IDs (not required)
%                       'geneNAMES': cell array containing gene names (not required)
%                       'timeline': 1xN numerical vector containing the time information (not required)
INPUTS.zeros_genes=3; % Remove genes with > 3 zero values
INPUTS.zeros_cells=200; % Remove cells with zero values>200
INPUTS.top_genes=300; % Retain only top 300 most variable genes
%
% Specify single-cell clustering settings
INPUTS.optimize=1; % The number of cluster is known a priori
INPUTS.max_clusters=20; % Max num of clusters
% 

% Upload and pre-process data 
[DATA,INPUTS] = import_data(INPUTS);

%% *** 2-SINGLE-CELL CLUSTERING ***
%
% Type 'help CALISTA_clustering_main' for more information. 
[Results, DATA, INPUTS]=CALISTA_clustering_main(DATA,INPUTS);

figure
hist(Results.final_groups,1:Results.expected_clusters)
title('Cells in each cluster')
xlabel('Cluster')
ylabel('Num of cells')
%% Upload marker genes

[FileName,PathName,FilterIndex] = uigetfile('*.*');
filename=strcat(PathName, FileName);
imported_data=importdata(filename);
[~,idx_marker_genes]=ismember(imported_data,DATA.original_genes);

testing_data=DATA.original_totDATA(:,idx_marker_genes);
mean_exp_marker_genes=[];
for i=1:Results.expected_clusters
    temp_data=testing_data(find(Results.final_groups==i),:);
    mean_exp_marker_genes(i,:)=mean(temp_data);%sum(temp_data, 1) ./ sum(temp_data~=0, 1);%mean(temp_data);%sum(temp_data, 1) ./ sum(temp_data~=0, 1);%
    
end
 mean_exp_marker_genes(isnan( mean_exp_marker_genes))=0;
 % Normalize each gene marker by the max value
 mean_exp_marker_genes=100*mean_exp_marker_genes./repmat(max(mean_exp_marker_genes),Results.expected_clusters,1);
%  cgo = clustergram(mean_exp_marker_genes,'Cluster',1,'Colormap',redbluecmap,'Standardize','column')
figure
imagesc(mean_exp_marker_genes([8 7 3 4 1 2 9 6 5],:))
colormap(flipud(gray))
set(gca, 'XTick', 1:1:length(imported_data), 'XTickLabel',imported_data,'FontSize',15);
xtickangle(45)
clusterID={'Neuron 1','Neuron 2','Oligo','Schwann','Meningeal','Astrocyte','Vascular','OPC','Microglia'};
set(gca, 'YTick', 1:1:Results.expected_clusters, 'YTickLabel',clusterID);
% ytickangle(45)