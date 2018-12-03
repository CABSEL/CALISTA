function Results=CALISTA_path_main(DATA,INPUTS,Results,varargin)
%CALISTA_PATH_MAIN perform post-analysis along developmental path(s) 
% This function gives the top 100 transition genes (or selected genes), performs hiearchical clustering,
% smoothen gene expression and constructs gene co-expression network for
% each path defined by the users.
%
% Usage:
%
% [Results]=CALISTA_path_main(INPUTS,Results)
%
% CALISTA will ask the user to specify the paths for the post-analysis.
% Each path consists of a sequence of clusters in the lineage progression
% graph. 
%
% Inputs:
% INPUTS - a structure containing the settings for the post-analysis with
% following fields:
%
% INPUTS.plot_fig (only for INPUTS.data_type<5)
% 1- plot figure of smoothened gene expression along path 
% 0- otherwise (by default)
%
% INPUTS.hclustering (only for INPUTS.data_type<5)
% 1- perform hierarchical clustering of gene expression for each path 
% 0- otherwise (by default)
%
% INPUTS.method (only for INPUTS.data_type<5)
% 1- use partial correlation for the gene co-expression network (value_cutoff=0.4, pvalue_cutoff=0.05)
% 2- use pairwise correlation for the gene co-expression network
% (value_cutoff=0.8, pvalue_cutoff=0.01) (by default)
%
% INPUTS.moving_average_window - set the size of window (percent of cells in each path)
% used for the moving averaging. 
% INPUTS.moving_average_window=10 (by default)
%
% INPUTS.path_auto
% 1- Automatically detect paths along the graph based on shortestpath (by default)
% 0- Otherwise
% Results - a structure of CALISTA clustering,lineage inference, transition genes detection and cell ordering results
% Run 'CALISTA_clustering_main',
%     'CALISTA_transition_main',
%     'CALISTA_transition_genes_main', and
%     'CALISTA_cell_ordering_main'
% 
% Outputs:
% Results.PATH.path_transition_genes - a 1xP cell array of P paths. The 
% i-th cell contains the transition genes of the i-th path.
% 
% Results.PATH.smoothExpr - a 1xP cell array of P paths. The i-th cell
% gives the moving-averaged single-cell expression matrix where each row
% corresponds to a gene.
%
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. June 1, 2018.


if nargin <3
    error('Not enough input arguments')
end

if INPUTS.data_type>=5
    INPUTS.CV_plot=0;
    INPUTS.plot_fig=0;
    INPUTS.hclustering=0;
end


if ~isfield(INPUTS,'moving_average_window')
    INPUTS.moving_average_window=10;
end

if ~isfield(INPUTS,'CV_plot')
    INPUTS.CV_plot=0;
end

if ~isfield(INPUTS,'plot_fig')
    INPUTS.plot_fig=0;
end

if ~isfield(INPUTS,'hclustering')
    INPUTS.hclustering=0;
end

if ~isfield(INPUTS,'path_auto')
    INPUTS.path_auto=1;
end

if ~isfield(INPUTS,'method')
    INPUTS.method=2;
end

nVarargs=length(varargin);
for k=1:2:nVarargs
    if strcmp(varargin{k},'selected_genes')
        % 1 assignement to calculate optimal rate constants
        selected_genes=varargin{k+1};
        % do not optimize
    end
end


if ~exist('selected_genes','var')
    [Results]=CALISTA_path(Results,INPUTS);
    Results=CALISTA_net_path(Results,INPUTS);
    
    % cytoscape
    for i=1:length(Results.PATH.Theta)
        connectivityMATRIX=tril(Results.PATH.Theta{i});
        genes=Results.GENES.tot_transition_genes;
        fileNAME=[' Table for Cytoscape ' num2str(i) '.xlsx'];
        table_for_cytoscape(connectivityMATRIX,genes,fileNAME);
    end
else
    if length(selected_genes)>100
        error('ERROR: too many genes. Please reduce the number of selected genes')
    end
    %% Scale and normalize raw data
    if INPUTS.cluster_time==1
        totDATA=DATA.TIME_DATA{1}.ORIGINAL_DATA.totDATA_raw_reduced;
        tot_genes=DATA.TIME_DATA{1}.ORIGINAL_DATA.genes_raw_reduced;
    else
        totDATA=DATA.totDATA_raw_reduced;
        tot_genes=DATA.genes_raw_reduced;
    end
    
    [~,idx_selected_genes]=ismember(upper(selected_genes),upper(tot_genes));
    
    if INPUTS.data_type==5
        tot_mRNA_each_cell=sum(totDATA,2);
        new_tot_mRNA_each_cell=median(tot_mRNA_each_cell);
        scaling=repmat(new_tot_mRNA_each_cell./tot_mRNA_each_cell,1,size(totDATA,2));
        scaled_totDATA=totDATA.*scaling;
        totDATA=scaled_totDATA;
    end
    
    selected_data=totDATA(:,idx_selected_genes);
    [DATA, selected_data] = internal_norm_calista(DATA,INPUTS.data_type,selected_data);
    
    
    [Results]=CALISTA_path_selected_genes(Results,INPUTS,selected_data,selected_genes);
    
    Results.PATH.selected_data=selected_data;
    Results.PATH.selected_genes=selected_genes;
end








% % ** CV_plot **
% % 1- calculate and plot the coefficients of variation
% % 0- otherwise (DEFAULT)
% %
% INPUTS.CV_plot=0;
% 
% [Results]=CALISTA_path(Results,INPUTS);
% Results=CALISTA_net_path(Results,INPUTS);
% 
% % cytoscape
% for i=1:length(Results.PATH.Theta)
%     connectivityMATRIX=tril(Results.PATH.Theta{i});
%     genes=Results.GENES.tot_transition_genes;
%     fileNAME=[' Table for Cytoscape ' num2str(i) '.xlsx'];
%     table_for_cytoscape(connectivityMATRIX,genes,fileNAME);
% end
% end

fprintf('\nDONE.\n')