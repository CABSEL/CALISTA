function [Results]=CALISTA_ordering_main_2(DATA,INPUTS,Results,cell_assignments)
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
% Copyright.  October 1, 2017.

% Error check
if nargin <3
    error('Not enough input arguments')
end

if ~isfield(INPUTS,'transition_new')
     INPUTS.transition_new=1;
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

fprintf('\nCALISTA_ordering is running...\n')
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
    Z_temp=Z_temp';
    Results.ORDERING.drop_prob=P;
%     Z_temp=log(Z_temp)';
else
    Z_temp=DATA.Parameters.P;
end


cell_ordering=zeros(DATA.nvars,2);
max_log_P_cell_all=zeros(DATA.nvars,1);
placement=zeros(DATA.nvars,1);
edges_of_cell=zeros(DATA.nvars,2);
n_points_interp=100;
% p=Results.TRANSITION.final_graph_force_layout;
h=Results.TRANSITION.final_graph;
nodes_connection3=[h.Edges.EndNodes(:,1) h.Edges.EndNodes(:,2)];
n_edges=size(nodes_connection3,1);

for clust=1:Results.expected_clusters
    max_log_P_each_edge=[];
    time_max_log_P_each_edge=[];
    idx_max_log_P_each_edge=[];
    actual_cluster=clust;
    edges_2_check=find(nodes_connection3(:,1)==actual_cluster | nodes_connection3(:,2)==actual_cluster);
    idx_cells_actual_clust=find(Results.final_groups==clust)';    
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
        
        sum_log_P_cell=zeros(n_points_interp,length(idx_cells_actual_clust),1);
        for j=1:DATA.numGENES
            int1=Z_temp(:, Results.clustering_struct.all.all.param_idx(clust1,j))';
            int2=Z_temp(:, Results.clustering_struct.all.all.param_idx(clust2,j))';
            xnew=interp1([t1 t2],[int1; int2],t_interp);
            exp_cells_actual_clust_actual_gene=DATA.totDATA(idx_cells_actual_clust,j)+1;
            
            sum_log_P_cell=sum_log_P_cell+log2(xnew(:,exp_cells_actual_clust_actual_gene));
        end        
        
        [max_log_P_each_edge(i,:),idx_max_log_P_each_edge(i,:)]=max(sum_log_P_cell);
        if same_t
            time_max_log_P_each_edge(i,:)=repmat(tequi,1,length(idx_cells_actual_clust));
        else
            time_max_log_P_each_edge(i,:)=t_interp(idx_max_log_P_each_edge(i,:));
        end
        
     end
     if length(edges_2_check)==1
         max_log_P_cell=max_log_P_each_edge;
         idx_max_log_P_cell=ones(1,length(max_log_P_cell));
         cell_ordering(idx_cells_actual_clust,1)=time_max_log_P_each_edge';
          placement(idx_cells_actual_clust)=idx_max_log_P_each_edge';
     else
         [max_log_P_cell,idx_max_log_P_cell]=max(max_log_P_each_edge);
         A=time_max_log_P_each_edge';
         A2=idx_max_log_P_each_edge';
         B=idx_max_log_P_cell;
         I = (1 : size(A, 1)) .';
         J = reshape(B, [], 1);
         k = sub2ind(size(A), I, J);
         C=A(k);
         cell_ordering(idx_cells_actual_clust,1) = C; 
         placement(idx_cells_actual_clust) = A2(k); 
     end
         max_log_P_cell_all(idx_cells_actual_clust)=max_log_P_cell;
         cell_ordering(idx_cells_actual_clust,2)=edges_2_check(idx_max_log_P_cell');
         edges_of_cell(idx_cells_actual_clust,:)=nodes_connection3(cell_ordering(idx_cells_actual_clust,2),:);
%          edges_of_cell(idx_cells_actual_clust,2)=nodes_connection3(edges_2_check(idx_max_log_P_cell'),2);
     
    
end
% cell_ordering=in_silico_pseudo;
% cell_ordering=monocle_pseudo;
normed_cell_ordering=(cell_ordering(:,1)-min(cell_ordering(:,1))/max(cell_ordering(:,1)-min(cell_ordering(:,1))));

x_center=zeros(Results.expected_clusters,1);
y_center=Results.TRANSITION.y_center;
z_center=Results.TRANSITION.z_center;
for i=1:Results.expected_clusters
    x_center(i)=mean(normed_cell_ordering(Results.final_groups==i));    
end



if ishandle(104)
clf(104)
end
hfig=figure(104);
% set(hfig,'Position',[500, 600, 650, 650])
p = plot(h);
layout(p,'force','XStart',p.XData,'YStart',1-x_center);
hold on
p.EdgeAlpha=1;
p.NodeColor=Results.colorMARK_calista;
p.LineWidth=0.001;
% p.EdgeColor=[0 0 0];
p.NodeLabel=[];
p.MarkerSize=10;
p.EdgeLabel=[];%abs(h.Edges.Weight);
p.ZData=1-x_center;
for i=1:Results.expected_clusters
    text(p.XData(i)+0.1,p.YData(i)+0.1,p.ZData(i)+0.1,num2str(i),'FontSize',15);
end
title('Cell Ordering Graph')
set(gca,'xticklabel',{[]}) 
set(gca,'yticklabel',{[]}) 

% Place cells on the progression graph
xnew=interp1([1 100],[p.YData(edges_of_cell(:,1)); p.YData(edges_of_cell(:,2))],1:100);
A=xnew';
B=placement';
I = (1 : size(A, 1)) .';
J = reshape(B, [], 1);
k = sub2ind(size(A), I, J);
Ycoordinate=A(k);

xnew=interp1([1 100],[p.XData(edges_of_cell(:,1)); p.XData(edges_of_cell(:,2))],1:100);
A=xnew';
Xcoordinate=A(k);


xnew=interp1([1 100],[p.ZData(edges_of_cell(:,1)); p.ZData(edges_of_cell(:,2))],1:100);
A=xnew';
Zcoordinate=A(k);


% Genenrate colormap based on pseudotime
[unique_pseudotime,~,idx_colormap]=unique(normed_cell_ordering);
colormap_pseudotime=parula(length(unique_pseudotime));
Results.ORDERING.orderingCOLORMARK=colormap_pseudotime(idx_colormap,:);
Results.ORDERING.colormap_pseudotime=colormap_pseudotime;
% figure
% plot(Xcoordinate,Ycoordinate,'o')
scatter3(Xcoordinate,Ycoordinate,Zcoordinate, 30, Results.ORDERING.orderingCOLORMARK, 'fill')

fprintf('\nWaiting for plots...\n')

if ishandle(103)
clf(103)
end
hfig=figure(103);
% az=-37.5000;
% el=30;
for k=1:Results.expected_clusters
    scatter3(normed_cell_ordering(Results.final_groups==Results.cluster_predicted(k))',Results.score3(Results.final_groups==Results.cluster_predicted(k),1),Results.score3(Results.final_groups==Results.cluster_predicted(k),2),30,Results.colorMARK_calista(k,:),'fill')
    hold on
end
grid on
p=plot(h);
% view([az,el])
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
Results.ORDERING.Xcoordinate=Xcoordinate;
Results.ORDERING.Ycoordinate=Ycoordinate;
Results.ORDERING.Zcoordinate=Zcoordinate;
Results.ORDERING.max_log_P_cell=max_log_P_cell_all;
Results.ORDERING.cell_ordering=cell_ordering;
Results.ORDERING.normed_cell_ordering=normed_cell_ordering;
Results.ORDERING.edges_of_cell=edges_of_cell;
Results.ORDERING.placement=placement;
Results.ORDERING.x_center=x_center;
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
fprintf('\nDONE.\n')
pause(1)