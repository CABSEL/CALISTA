function [h,NumberOfEdges,nodes]=AddGraph(nodes,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,az,el,cluster_progression,DATA,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if size(varargin,1)>0
    h=varargin{1};
    NumberOfEdges=size(h.Edges,1);
else
    h=[];
    h=graph;
    non=length(unique(nodes(:,1:2)));
    h=addnode(h,non);
%     NumberOfEdges=non-1;
    
    
   
    % If there's no time info, start the graph with one edge
    if length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0
        NumberOfEdges=1;
        for j=1:NumberOfEdges
            h=addedge(h, nodes(j,1),nodes(j,2),nodes(j,3));
        end
    else % add one incoming edge for each descendant node
        TF = isoutlier(nodes(:,3));
        nodes_temp=nodes;
        nodes_added=[];
        NumberOfEdges=0;
        unique_cluster_progression=unique(cluster_progression);
        jj_tot=[];
        normalized_timeline=(DATA.timeline-min(DATA.timeline))/(max(DATA.timeline)-min(DATA.timeline));
        bin_tot_time_edge=unique(normalized_timeline);
        for i=1:expected_clusters
            idx_cells_in_k=find( ClusterGroup2==i);
            labels_in_k=normalized_timeline(idx_cells_in_k);
            occurrence_tot_time(i,:) = hist( labels_in_k,bin_tot_time_edge);
        end
        for j=2:length(unique_cluster_progression)
            idx_actual_clusters=find(cluster_progression==unique_cluster_progression(j));
            for k=1:length(idx_actual_clusters)
                cells_in_cluster=find(ClusterGroup2==idx_actual_clusters(k));
                actual_times_in_cluster=normalized_timeline(cells_in_cluster);
                bin_edge=unique(actual_times_in_cluster);
%                 occurrence=hist(actual_times_in_cluster,bin_edge);
%                 % find time info for cells that occurr at least >5% of the total cells in the cluster
%                 selected_actual_times_in_cluster=bin_edge(find((occurrence/length(cells_in_cluster))>.10));
                % remove clusters at actual/future cluster pseudotime
                idx_excluded_clusters=find(cluster_progression>=unique_cluster_progression(j));
%                 % find the cluster that has more cells at the previous time
%                 % point
                [~, index_closest_time_point] = min(abs(prctile(actual_times_in_cluster,5)-bin_edge));    % time point at 5 percentile    
                exit_while=0;
                temp_occurrence_tot_time=occurrence_tot_time(:,max(1,find(bin_tot_time_edge==bin_edge(index_closest_time_point))-1));
                temp_occurrence_tot_time(idx_excluded_clusters)=0;
                while(~exit_while)                    
                    [~,idx_thr]=max(temp_occurrence_tot_time);
                    [~, idx_edge_to_check]=ismember([idx_thr idx_actual_clusters(k)],nodes(:,1:2),'rows');
                    if TF(idx_edge_to_check)==0
                        exit_while=1;
                    else
                         temp_occurrence_tot_time(idx_thr)=0; %find the next max
                    end
                end
                thr=find(unique_cluster_progression==cluster_progression(idx_thr));
                idx_previous_clusters=find(cluster_progression<unique_cluster_progression(j) & cluster_progression>=unique_cluster_progression(max(1,thr)));
                target=find(nodes(:,2)==idx_actual_clusters(k));
                sources=find(ismember(nodes(:,1),idx_previous_clusters)==1);
                idx_intersected_nodes=intersect(sources,target);
                [~,idx_edge_2_add]=max(nodes(idx_intersected_nodes,3));
                jj=idx_intersected_nodes(idx_edge_2_add);
                h=addedge(h, nodes(jj,1),nodes(jj,2),nodes(jj,3));
                jj_tot=[jj_tot jj];
                NumberOfEdges=NumberOfEdges+1;
                nodes_added=[ nodes_added; nodes(jj,:)];
            end
            
        end
        nodes_temp(jj_tot,:)=[];
        nodes=[nodes_added; nodes_temp];
    end
end


for k=1:expected_clusters
    scatter3(normed_pseudotime(ClusterGroup2==k)',score3(ClusterGroup2==k,1),score3(ClusterGroup2==k,2),30,colorMARK2(k,:),'fill')
    hold on
end
legend(legendInfo,'Location', 'northeast')
xlabel('Cluster pseudotime')
ylabel('COMP1')
zlabel('COMP2')
grid on
az=-26;
el=8;
p=plot(h,'EdgeLabel',h.Edges.Weight);
view([az,el])
p.MarkerSize=1;
% p.ArrowSize=15;
p.EdgeAlpha=1;
p.NodeColor=colorMARK2;
p.LineWidth=2.5;
p.EdgeColor=[0 0 0];
p.XData=x_center;
p.YData=y_center;
p.ZData=z_center;
p.NodeLabel=[];
p.EdgeLabel=abs(h.Edges.Weight);
% colorMARK3=bone(size(nodes,1));
% [~,IDXCOLORMARK3]=ismember(h.Edges.Weight,nodes(:,3));
% colorMARK3=colorMARK3(1:size(h.Edges,1),:);
% colorMARK3=flipud(colorMARK3);
% p.EdgeColor=[0 0 0];%colorMARK3(IDXCOLORMARK3,:);

end

