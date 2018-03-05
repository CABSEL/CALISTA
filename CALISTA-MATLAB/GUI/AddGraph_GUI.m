function [h,NumberOfEdges,IDXCOLORMARK3]=AddGraph_GUI(cluster_distance,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,checked_boxes)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
h=[];
if ~isempty(checked_boxes)
    
    h=graph;
    unique_nodes=unique(cluster_distance(:,1:2));
    non=length(unique_nodes);
    h=addnode(h,non);
    cl_dist=cluster_distance(checked_boxes,:);
    NumberOfEdges=size(cl_dist,1);
    %     NumberOfEdges=1;
    
    for j=1:NumberOfEdges
        h=addedge(h, cl_dist(j,1),cl_dist(j,2),cl_dist(j,3));
    end
    
    
    
%     for k=1:expected_clusters
%         scatter3(normed_pseudotime(ClusterGroup2==k)',score3(ClusterGroup2==k,1),score3(ClusterGroup2==k,2),30,colorMARK2(k,:),'fill')
%         hold on
%     end
%     legend(legendInfo,'Location', 'northeast')
%     xlabel('Cluster pseudotime')
%     ylabel('COMP1')
%     zlabel('COMP2')
%     grid on
    
    p=plot(h,'EdgeLabel',h.Edges.Weight);
    % view([az,el])
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
%     colorMARK3=pink(size(cluster_distance,1));
%     [~,IDXCOLORMARK3]=ismember(h.Edges.Weight,cluster_distance(:,3));
%     colorMARK3=colorMARK3(1:size(h.Edges,1),:);
    % colorMARK3=flipud(colorMARK3);
%     p.EdgeColor=colorMARK3(IDXCOLORMARK3,:);
    
end
end