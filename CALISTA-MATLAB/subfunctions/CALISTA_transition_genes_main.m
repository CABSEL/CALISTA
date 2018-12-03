function [Results]=CALISTA_transition_genes_main(DATA,INPUTS,Results,cell_assignments)
%CALISTA_TRANSITION_GENE_MAIN identify the key genes in lineage progression
% Given a lineage progression graph, CALISTA determines the key transition 
% genes for any two connected clusters in the graph, based on the gene-wise 
% likelihood difference between having the cells separately as two clusters 
% and together as a single cluster. Larger differences in the gene-wise
% likelihood point to more informative genes.
% 
% Usage:
% 
% 1- Obtain transition genes after CALISTA lineage inference:
%    [Results]=CALISTA_transition_genes_main(DATA,INPUTS,Results)
% 
% 2- Obtain transition genes using user-defined clustering
%    [Results]=CALISTA_transition_genes_main(DATA,INPUTS,[],cell_assignments)
% CALISTA will ask users to specify sequences of connected clusters, i.e. 
% paths, in the lineage graph. The list of edges in the graph is the union 
% of all edges in the user-specified paths.
%    
% Inputs:   
% DATA - a structure containing preprocessed single-cell expression data
% Use 'import_data' to upload and preprocess single-cell expression values.
% 
% In addition to the specification in 'import_data' and 'CALISTA_clustering_main',
% users need to specify:
%
% ** INPUTS.thr_transition_genes **
% the percentile for the cumulative gene-wise likelihood difference up to
% which genes are included in the set of transition genes.
%
% ** INPUTS.plot_trans_genes
% 1- Plot top transition genes for each inferred edge
% 0- Otherwise (by default)
%
% Results - a structure of CALISTA clustering and lineage inference results
% Run 'CALISTA_clustering_main' and/or 'CALISTA_transition_main'
%
% cell_assignments - 1xN vector of INTEGERS with N = number of cells. 
% The n-th element of cell_assignments contains the cluster assignment of
% the n-th cell of the expression data uploaded. Cluster names must
% be assigned in sequence (e.g. 1,2,3,4 and not 1,2,4).
%
% Outputs:
% Results - a structure containing the results of CALISTA analysis. 
% The most relevant fields containing the transition genes are:
%
% Results.GENES.final_transition_genes - a 1xE cell array (E: the number of 
% edges). The i-th cell contains the names of the transition genes for the
% edge. The edges are defined in Results.TRANSITION.nodes_connection.
% 
% Results.GENES.tot_transition_genes - a cell array containing the names of 
% all transition genes.
% 
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright.  October 1, 2018.

% Error check
if nargin <3
    error('Not enough input arguments')
end

if ~isfield(INPUTS,'transition_new')
     INPUTS.transition_new=1;
end

if ~isfield(INPUTS,'plot_trans_genes')
     INPUTS.plot_trans_genes=0;
end

if isempty(Results)
    if nargin <4
        error('Not enough input arguments')
    end
    if ~isvector(cell_assignments)
        error('Please upload the cell assignments as a vector')
    end
    Results=jump_clustering(DATA,INPUTS,cell_assignments);
    
    if INPUTS.transition_new==1
        [Results]=CALISTA_transition_new(DATA,Results);
    else
        Results=jump_transition(DATA,Results);
        % *** 3.b-Plot mean expressions for each cluster ***
        Results=plot_cluster_mean_exp(DATA,Results);
        
        % *** 3.c-Cell-Cell variability analysis ***
        [Results]=cell_variability(DATA,Results);
    end
    
    
end

if ~isfield(INPUTS,'thr_transition_genes')
    if INPUTS.data_type>=5
        INPUTS.thr_transition_genes=25;
    else
        INPUTS.thr_transition_genes=50;
    end
end

if ~isfield(INPUTS,'plot_trans_genes')
    if INPUTS.data_type>=5
        INPUTS.plot_trans_genes=0;
    else
        INPUTS.plot_trans_genes=1;
    end
end

if ~isfield(INPUTS,'use_drop_prob_in_clustering')
    if INPUTS.data_type<5
        INPUTS.use_drop_prob_in_clustering=0;
    else
        INPUTS.use_drop_prob_in_clustering=1;
    end
end

fprintf('\nCALISTA_transition_genes is running...\n')
Parameters=DATA.Parameters;
thr=INPUTS.thr_transition_genes;
my_results_final=Results.clustering_struct;
hh=Results.TRANSITION.final_graph;
nodes_connection2=[hh.Edges.EndNodes(:,2) hh.Edges.EndNodes(:,1)];

if INPUTS.plot_trans_genes==1
    max_plots_in_each_fig=9;
    count2=1;
    num_edges=size(nodes_connection2,1);
    [p,~]=numSubplots(max_plots_in_each_fig);
    if num_edges>max_plots_in_each_fig
        figs=ceil(num_edges/12);
    else
        figs=1;
    end
end



if INPUTS.use_drop_prob_in_clustering
    p_mat_ma=DATA.Parameters.P;
    opt_lambda=DATA.opt_lambda;
    % invert the probability matrix  %NOT LOG AS IN TRANSITION GENES!
    X=p_mat_ma';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define dropout probability
    epsilon=min(min(X(X>0)));
    mRNA_counts=[0 1:200];
    P_temp=exp(-opt_lambda*mRNA_counts)';
    P=diag([1; 1-P_temp(2:end)]);
    P(:,1)=P_temp;
    Z_temp=X*P;
    Z_temp(Z_temp==0)=epsilon;  %%%%%%% SET ZEROS TO A SMALL VALUE
    Z_temp=log(Z_temp)';
else
    Z_temp=DATA.Parameters.logP;
end

for i=1:size(nodes_connection2,1)
    count=1;
    temp_mRNA_all=[];
    selected_clusters=[];
    for j=1:size(nodes_connection2,2)
        
        [~,idx_final_groups_temp]=find(Results.final_groups==nodes_connection2(i,j));
        temp_mRNA_all=[temp_mRNA_all; DATA.totDATA(idx_final_groups_temp,:)];
        selected_clusters=[selected_clusters count*ones(1,length(idx_final_groups_temp))];
        count=count+1;
    end
    nvars_temp=size(temp_mRNA_all,1);
    [prob_transition_genes(1:size(nodes_connection2,2),:,i)]=get_prob_transition_genes( selected_clusters,Z_temp,temp_mRNA_all,DATA.numGENES,nvars_temp);
    [prob_transition_genes(size(nodes_connection2,2)+1,:,i)]=get_prob_transition_genes( ones(1,nvars_temp),Z_temp,temp_mRNA_all,DATA.numGENES,nvars_temp);
    prob_separated_clusters=sum(prob_transition_genes(1:size(nodes_connection2,2),:,i));
    prob_all_cells=prob_transition_genes(size(nodes_connection2,2)+1,:,i);
    [sorted_gene_prob(i,:),idx_transition_genes(i,:)]=sort(prob_separated_clusters-prob_all_cells,'descend');
    transition_gene_ranking(i,:)=DATA.genes(idx_transition_genes(i,:));
    null_LL(i)=sum(prob_transition_genes(size(nodes_connection2,2)+1,:,i));
    nodes_connection2(i,:)=sort(nodes_connection2(i,:));
    
    cum_LL=cumsum(sorted_gene_prob(i,:));
    idx_thr=find(cum_LL<=(thr*max(cum_LL)/100),1,'last');
    num_transition_genes(i)=idx_thr;
    final_transition_genes{i}=transition_gene_ranking(i,1:num_transition_genes(i));
    
    if INPUTS.plot_trans_genes==1
        
        figure(4000+ceil(i/max_plots_in_each_fig))
%         set(gcf,'units','points','position',[100,100,2000,1000])
        subplot(p(1),p(2),count2)
        bar(sorted_gene_prob(i,1:num_transition_genes(i)))
        %     grid on
        xlim([0 num_transition_genes(i)+1])
        set(gca,'xtick',1:num_transition_genes(i),'xticklabel',final_transition_genes{i})
        ylabel('v^g_{j k}')
        ax = gca;
        ax.XTickLabelRotation = -40;
        title([num2str(nodes_connection2(i,1)) ' - ' num2str(nodes_connection2(i,2))])
        legend([num2str(num_transition_genes(i)) ' transition genes '],'Location', 'northeast')
        if count2==max_plots_in_each_fig
            count2=1;
        else
            count2=count2+1;
        end
    end
end

for transition_num=1:size(nodes_connection2,1)
    [~,idx_final_transition_genes{transition_num}]=ismember(final_transition_genes{transition_num},DATA.genes);
    for i=1:num_transition_genes(transition_num)
        for j=1:size(nodes_connection2,2)
            transition_gene_parameters{transition_num}(:,j,i)=my_results_final.all.all.parameter{nodes_connection2(transition_num,j)}(idx_final_transition_genes{transition_num}(i),:)';
        end
    end
end
idx_tot_transition_genes=unique(cat(2, idx_final_transition_genes{:}));
tot_transition_genes=DATA.genes(idx_tot_transition_genes);

Results.GENES.mRNA_tot_transition_genes=DATA.totDATA(:,idx_tot_transition_genes);

Results.GENES.tot_transition_genes=tot_transition_genes;
Results.GENES.idx_tot_transition_genes=idx_tot_transition_genes;
Results.GENES.thr=thr;
Results.GENES.final_transition_genes=final_transition_genes;
Results.GENES.transition_gene_parameters=transition_gene_parameters;
fprintf('\nDONE.\n')
pause(1)

