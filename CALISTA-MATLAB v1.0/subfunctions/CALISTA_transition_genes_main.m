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
% In addition to the specification in 'import_data', users need to specify:
% INPUTS.thr_transition_genes - the percentile for the cumulative gene-wise
% likelihood difference up to which genes are included in the set of
% transition genes.
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
% Copyright. June 1, 2017.

% Error check
if nargin <3
    error('Not enough input arguments')
end 

if isempty(Results)
   if nargin <4
    error('Not enough input arguments')
   end 
   if ~isvector(cell_assignments)
       error('Please upload the cell assignments as a vector')
   end
      Results=jump_clustering(DATA,cell_assignments);
      Results=jump_transition(DATA,Results);
      % [Results]=CALISTA_transition(DATA,Results);
      
      % *** 3.b-Plot mean expressions for each cluster ***
      Results=plot_cluster_mean_exp(DATA,Results);
      
      % *** 3.c-Cell-Cell variability analysis ***
      [Results]=cell_variability(DATA,Results);

end

fprintf('\nCALISTA_transition_genes is running...\n')
Parameters=DATA.Parameters;
thr=INPUTS.thr_transition_genes;
my_results_final=Results.clustering_struct;
hh=Results.TRANSITION.final_graph;
nodes_connection2=[hh.Edges.EndNodes(:,2) hh.Edges.EndNodes(:,1)];
[p,~]=numSubplots(size(nodes_connection2,1));
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
    [prob_transition_genes(1:size(nodes_connection2,2),:,i)]=get_prob_transition_genes( selected_clusters,Parameters.logP,temp_mRNA_all,DATA.numGENES,nvars_temp);
    [prob_transition_genes(size(nodes_connection2,2)+1,:,i)]=get_prob_transition_genes( ones(1,nvars_temp),Parameters.logP,temp_mRNA_all,DATA.numGENES,nvars_temp);
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
    
    figure(4000)
    set(gcf,'units','points','position',[100,100,2000,1000])
    subplot(p(1),p(2),i)
    bar(sorted_gene_prob(i,1:num_transition_genes(i)))
    %     grid on
    xlim([0 num_transition_genes(i)+1])
    set(gca,'xtick',1:num_transition_genes(i),'xticklabel',final_transition_genes{i})
    ylabel('v^g_{j k}')
    ax = gca;
    ax.XTickLabelRotation = -40;
    title([num2str(nodes_connection2(i,1)) ' - ' num2str(nodes_connection2(i,2))])
    legend([num2str(num_transition_genes(i)) ' transition genes '],'Location', 'northeast')
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
pause(3)

