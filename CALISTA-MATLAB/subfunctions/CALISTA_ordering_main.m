function [Results]=CALISTA_ordering_main(DATA,Results,cell_assignments)
%CALISTA_ORDERING_MAIN perform pseudotemporal ordering of cells
% For pseudotemporal ordering of cells, CALISTA performs maximum
% likelihood optimization for each cell using a linear interpolation of the
% cell likelihoods between any two connected clusters. 
% 
% Usage:
% 
% 1- Run CALISTA cell ordering using results of CALISTA clustering and lineage inference
%    [Results]=CALISTA_ordering_main(DATA,Results);
%  
% 2- Run CALISTA cell ordering using user-defined cell assignments
%    [Results]=CALISTA_ordering_main(DATA,[],cell_assignments)
%
% CALISTA will ask users to specify sequences of connected clusters, i.e. 
% paths, in the lineage graph. The list of edges in the graph is the union 
% of all edges in the user-specified paths.
%    
% Inputs:   
% DATA - a structure containing preprocessed single-cell expression data
% Use 'import_data' to upload and preprocess single-cell expression values.
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
% The most relevant fields containing the cell ordering are:
%
% Results.ORDERING.normed_cell_ordering - a vector contaning the pseudotime
% of the cells.
% 
% Results.ORDERING.idx_sorted_cells - a vector containing the ordered cell
% indices in increasing order of the cell pseudotime.
% 
% Results.ORDERING.idx_actual_edge: a 1xE cell array with the i-th cell 
% containing the vector of cell indices assigned in the i-th edge. The
% edges are defined in Results.TRANSITION.nodes_connection.
% 
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. June 1, 2017.

% Error check
if nargin <2
    error('Not enough input arguments')
end 

if isempty(Results)
   if nargin <3
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


fprintf('\nCALISTA_ordering is running...\n')

Parameters=DATA.Parameters;
n_points_interp=100;
hh=Results.TRANSITION.final_graph;
[cell_prob_3D]=get_3D_cell_prob(Results,Parameters,DATA);
nodes_connection3=[hh.Edges.EndNodes(:,1) hh.Edges.EndNodes(:,2)];
n_edges=size(nodes_connection3,1);
cell_ordering=[];
for CELL=1:DATA.nvars
    actual_cluster=Results.final_groups(CELL);
    edges_2_check=find(nodes_connection3(:,1)==actual_cluster | nodes_connection3(:,2)==actual_cluster);
    max_log_P_each_edge=[];
    time_max_log_P_each_edge=[];
    for i=1:length(edges_2_check)
        clust1=find(Results.cluster_predicted==nodes_connection3(edges_2_check(i),1));
        clust2=find(Results.cluster_predicted==nodes_connection3(edges_2_check(i),2));
        t1=Results.cluster_progression(clust1);
        t2=Results.cluster_progression(clust2);
        same_t=0;
        if t1==t2
            %             error('Clusters with the same time info')
            same_t=1;
            tequi=t1;
            t1=0;
            t2=1;
            t_interp = linspace(t1,t2,n_points_interp);
        else
            t_interp = linspace(t1,t2,n_points_interp);
        end
        sum_log_P_genes=zeros(1,length(t_interp));
        for j=1:DATA.numGENES
            int1=squeeze(cell_prob_3D(CELL,clust1,j));
            int2=squeeze(cell_prob_3D(CELL,clust2,j));
            xnew=interp1([t1 t2],[int1 int2],t_interp);
            sum_log_P_genes=sum_log_P_genes+log2(xnew);
        end
        [max_log_P_each_edge(i),idx_max_log_P_each_edge]=max(sum_log_P_genes);
        if same_t
            time_max_log_P_each_edge(i)=tequi;
        else
            time_max_log_P_each_edge(i)=t_interp(idx_max_log_P_each_edge);
        end
        
    end
    [max_log_P_cell(CELL),idx_max_log_P_cell]=max(max_log_P_each_edge);
    if same_t
        cell_ordering(CELL,1)= tequi;
    else
        cell_ordering(CELL,1)= time_max_log_P_each_edge(idx_max_log_P_cell);
    end
    cell_ordering(CELL,2)=edges_2_check(idx_max_log_P_cell);
end

normed_cell_ordering=(cell_ordering(:,1)-min(cell_ordering(:,1)))/max(cell_ordering(:,1)-min(cell_ordering(:,1)));

x_center=zeros(Results.expected_clusters,1);
y_center=Results.TRANSITION.y_center;
z_center=Results.TRANSITION.z_center;
for i=1:Results.expected_clusters
    x_center(i)=mean(normed_cell_ordering(Results.final_groups==i));    
end
figure
az=-37.5000;
el=30;
for k=1:Results.expected_clusters
    scatter3(normed_cell_ordering(Results.final_groups==Results.cluster_predicted(k))',Results.score3(Results.final_groups==Results.cluster_predicted(k),1),Results.score3(Results.final_groups==Results.cluster_predicted(k),2),30,Results.colorMARK_calista(k,:),'fill')
    hold on
end
grid on
p=plot(hh);
view([az,el])
p.MarkerSize=15;
p.EdgeAlpha=1;
p.NodeColor=Results.colorMARK_calista;
p.LineWidth=2.5;
p.XData=x_center;
p.YData=y_center;
p.ZData=z_center;
p.NodeLabel=[];
title('Cell Ordering')
xlabel('Cell Ordering')
ylabel('PC1')
zlabel('PC2')
% 
% legend(Results.legendInfo_calista,'Location', 'northeast')
% title('Transition states')
% xlabel('Cell ordering')
% ylabel('COMP1')
% zlabel('COMP2')

Results.ORDERING.max_log_P_cell=max_log_P_cell;
Results.ORDERING.cell_ordering=cell_ordering;
Results.ORDERING.normed_cell_ordering=normed_cell_ordering;
for i=1:n_edges
    % Define cells in edge
    idx_actual_edge{i}=find(Results.ORDERING.cell_ordering(:,2)==i);
    [~,bb]=sort(Results.ORDERING.cell_ordering(idx_actual_edge{i},1)); % sort based on cell ordering
    idx_actual_edge{i}=idx_actual_edge{i}(bb);
    cells_assigned_to_edge(i)=length(idx_actual_edge{i});
end
Results.ORDERING.cells_assigned_to_edge=cells_assigned_to_edge;
Results.ORDERING.idx_actual_edge=idx_actual_edge;
[~,Results.ORDERING.idx_sorted_cells]=sort(Results.ORDERING.cell_ordering(:,1));
% Results.ORDERING.totDATA_sorted=DATA.totDATA(Results.ORDERING.idx_sorted_cells,:);

pause(1)