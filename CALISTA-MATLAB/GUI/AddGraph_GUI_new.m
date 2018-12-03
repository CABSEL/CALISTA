function [h,NumberOfEdges]=AddGraph_GUI_new(Results,checked_boxes)

cluster_distance=Results.TRANSITION.cluster_distance;
colorMARK2=Results.colorMARK_calista;
x_center=Results.TRANSITION.x_center;

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
        
    p = plot(h);
    layout(p,'force','XStart',p.XData,'YStart',1-x_center);
%     az=-26;
%     el=8;
%     view([az,el])
    % sources=find(Results.cluster_progression==min(Results.cluster_progression));
    % sinks=find(Results.cluster_progression==max(Results.cluster_progression));
    % layout(p,'layered','Sources',sources,'Sinks',sinks)
    p.EdgeAlpha=1;
    p.NodeColor=colorMARK2;
    p.LineWidth=2.5;
    % p.EdgeColor=[0 0 0];
    p.NodeLabel=[];
    p.MarkerSize=20;
    p.EdgeLabel=[];%abs(h.Edges.Weight);
    p.ZData=1-x_center;
    for i=1:Results.expected_clusters
        text(p.XData(i)+0.1,p.YData(i)+0.1,p.ZData(i)+0.07,num2str(i),'FontSize',15);
    end
%     title('Lineage Progression Graph')
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    grid off
    axis off
    
end
end