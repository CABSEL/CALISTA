function [Results]=CALISTA_transition_new(DATA,Results)

fprintf('\nCALISTA_transition is running...\n\n')
if nargin <2
    error('Not enough input variables.')
end

if length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0 % no time info, don't plot cluster pseudotime here
    Results.legendInfo_calista_transition=Results.legendInfo_calista;
else
    stt1='Cluster:';
    stt2='Cluster pseudotime:';
    for clust=1:Results.expected_clusters
        Results.legendInfo_calista_transition{clust} = sprintf( '%s %4i %s %1.2f', stt1, clust,stt2,Results.cluster_progression(clust));
    end
end

Results.TRANSITION.clusterTEXT=0; %for CALISTA_path
my_results_final=Results.clustering_struct;
my_results_final.all.all.distance(my_results_final.all.all.distance==0)=inf;
[~,neighbour]=min(my_results_final.all.all.distance,[],2);
all_in_one=[neighbour Results.final_groups'];
all_in_one_sorted=sortrows(all_in_one,2);



for i=1:Results.expected_clusters
    for j=1:Results.expected_clusters
        aa(i,j)=sum(all_in_one_sorted(all_in_one_sorted(:,2)==i,1)==j);
    end
end

my_mean=zeros(Results.expected_clusters,Results.expected_clusters);
xxx=my_results_final.all.all.cell_prob./repmat(sum(my_results_final.all.all.cell_prob,2),1,Results.expected_clusters);
my_mean=zeros(Results.expected_clusters);
mean_prob_in_out_cluster=zeros(Results.expected_clusters,2);
for i=1:Results.expected_clusters
    aaa=find(Results.final_groups==i);
    my_mean(i,:)=mean(my_results_final.all.all.cell_prob(aaa,:));
    mean_prob_in_out_cluster(i,1)=my_mean(i,i);
    temp_my_mean=my_mean(i,:);
    temp_my_mean(i)=[];
    mean_prob_in_out_cluster(i,2)=mean(temp_my_mean);
    hold on
end

distance=(my_mean-diag(my_mean)*ones(1,Results.expected_clusters))./DATA.numGENES;

Results.TRANSITION.distance=abs(distance);


%%%%%%%%%%%%%%%%%%%%%%%%
%% Conditional transition

if  length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0 % no time info
    error('*** New CALISTA_transition is running. Please add time info and re-launch CALISTA ***')
else % time info present
    fprintf('\nConditional transition activated...\n\n')
    mask=-inf*ones(Results.expected_clusters);
    unique_cluster_progression=unique(Results.cluster_progression);
    for i=1:length(unique_cluster_progression)
        % remove edges between cell clusters at the same pseudotime
        cluster_selectected=find(Results.cluster_progression==unique_cluster_progression(i));
        distance(cluster_selectected,cluster_selectected)=-inf;
        % remove edges between cell clusters at abs(pseudotime)>2x
        allowed_edges_based_on_clust_progression=unique_cluster_progression(max(i-2,1):min(i+2,length(unique_cluster_progression)));
        [cluster_selectected2]=find(ismember(Results.cluster_progression,allowed_edges_based_on_clust_progression)==1);
        mask(cluster_selectected2,cluster_selectected2)=0;
        
    end
    distance=distance+mask;
    
    distance(logical(eye(size(distance)))) = 0;
    
end
%%%%%%%%%%%%%%%%%

%% Sort distances


[sorted_dist, idx_dist]=sort(distance(:),'descend');
idx_2_remove=find(sorted_dist==-inf | sorted_dist==0);
sorted_dist(idx_2_remove)=[];
idx_dist(idx_2_remove)=[];

I=zeros(length(idx_dist),3);
[I(:,2),I(:,1)] = ind2sub(Results.expected_clusters,idx_dist);
I(:,3)=sorted_dist;
nodes_all=I;
nodes_all(:,1)=min(I(:,1:2),[],2);
nodes_all(:,2)=max(I(:,1:2),[],2);
[~, ind_dir ]=unique(nodes_all(:,1:2),'rows');
nodes=nodes_all(ind_dir,:);
nodes=sortrows(nodes,-3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selected_temp_idx=[];
final_nodes=[];
count=1;
for i=1:Results.expected_clusters %for each cluster
    % select the IN-edge
    temp_idx=find(nodes(:,1)==i,1,'first');    
    if ~isempty(temp_idx) && isempty(find(selected_temp_idx==temp_idx))
        selected_temp_idx(count)=temp_idx;
        final_nodes(count,:)=nodes(temp_idx,:);
        count=count+1;
    end
    % select the OUT-edge
    temp_idx=find(nodes(:,2)==i,1,'first');    
    if ~isempty(temp_idx) && isempty(find(selected_temp_idx==temp_idx))
        selected_temp_idx(count)=temp_idx;
        final_nodes(count,:)=nodes(temp_idx,:);
        count=count+1;
    end
end
nodes=sortrows(final_nodes,-3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cluster_distance=nodes;
cluster_distance(:,3)=abs(cluster_distance(:,3));


x_center=zeros(Results.expected_clusters,1);
z_center=x_center;
y_center=x_center;
for i=1:Results.expected_clusters
    x_center(i)=mean(Results.cell_cluster_progression(Results.final_groups==i));
    y_center(i)=mean(Results.score3(Results.final_groups==i,1));
    z_center(i)=mean(Results.score3(Results.final_groups==i,2));
    
end
%
Results.TRANSITION.nodes=cluster_distance;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ishandle(100)
    clf(100)
end
hfig=figure(100);
% set(hfig,'Position',[600, 600, 600, 600])

[h,p]=AddGraph_new(nodes,Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista_transition);
figure(100)
title('Lineage Progression Graph')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove unwanted edges before additional analysis

hh=h;
proceed=input('Press 1 if you want to remove edges, 0 otherwise: ');
if ~isempty(proceed) & proceed==1
    
    figure(100)
   
    nodes_connection=[hh.Edges.EndNodes(:,1) hh.Edges.EndNodes(:,2)];
    add_paths=true;
    count=0;
    while add_paths
        
        
        fprintf('\n ******************************************************** \n')
        remove_edges=input(' Specify the node pairs (e.g. [4 5]: ');
        hh_temp = rmedge(hh,remove_edges(:,1),remove_edges(:,2));
        mininal_graph=length(dfsearch(hh_temp,1));
        if  mininal_graph==size(hh.Nodes,1)
            count=count+1;
            hh=hh_temp;
        else
            fprintf('WARNING: the chosen edge can not be removed. All nodes should be connected by at least one edge \n')
        end
        continue_add=input(' Press 1 to remove another edge, 0 otherwise: ');
        if continue_add~=1
            add_paths=false;
        end
    end
    fprintf('%4i %s',count,'edge to removed ')
    
    if ishandle(100)
        clf(100)
    end
    hfig=figure(100);
    % set(hfig,'Position',[500, 600, 650, 650])
    p = plot(hh);
    layout(p,'force','XStart',p.XData,'YStart',1-x_center);
    % sources=find(Results.cluster_progression==min(Results.cluster_progression));
    % sinks=find(Results.cluster_progression==max(Results.cluster_progression));
    % layout(p,'layered','Sources',sources,'Sinks',sinks)
    p.EdgeAlpha=1;
    p.NodeColor=Results.colorMARK_calista;
    p.LineWidth=2.5;
    % p.EdgeColor=[0 0 0];
    p.NodeLabel=[];
    p.MarkerSize=20;
    p.EdgeLabel=[];%abs(h.Edges.Weight);
    p.ZData=1-x_center;
    for i=1:Results.expected_clusters
        text(p.XData(i)+0.1,p.YData(i)+0.1,p.ZData(i)+0.07,num2str(i),'FontSize',15);
    end
    title('Lineage Progression Graph')
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    h=hh;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results.TRANSITION.x_center=x_center;
Results.TRANSITION.y_center=y_center;
Results.TRANSITION.z_center=z_center;

Results.TRANSITION.final_graph=h;
Results.TRANSITION.final_graph_force_layout=p;
hh=Results.TRANSITION.final_graph;
nodes_connection3=[hh.Edges.EndNodes(:,1) hh.Edges.EndNodes(:,2)];
Results.TRANSITION.nodes_connection=nodes_connection3;
% Results.TRANSITION.cluster_distance_tot=cluster_distance_tot;
Results.TRANSITION.cluster_distance=cluster_distance;
pause(2)
end