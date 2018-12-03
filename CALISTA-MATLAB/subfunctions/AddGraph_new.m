function [h,p]=AddGraph_new(nodes,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,varargin)
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
%     NumberOfEdges=1;
    h=addedge(h, nodes(:,1),nodes(:,2),abs(nodes(:,3)));
   
end

% colorMARK3=bone(size(nodes,1));
% [~,IDXCOLORMARK3]=ismember(h.Edges.Weight,nodes(:,3));
% colorMARK3=colorMARK3(1:size(h.Edges,1),:);
% colorMARK3=flipud(colorMARK3);
% p.EdgeColor=[0 0 0];%colorMARK3(IDXCOLORMARK3,:);


% set(hfig,'Position',[500, 600, 650, 650])
p = plot(h);
layout(p,'force','XStart',p.XData,'YStart',1-x_center);
az=-26;
    el=8;
view([az,el])
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
for i=1:expected_clusters
    text(p.XData(i)+0.1,p.YData(i)+0.1,p.ZData(i)+0.07,num2str(i),'FontSize',15);
end
set(gca,'xticklabel',{[]}) 
set(gca,'yticklabel',{[]}) 
grid off
axis off
if ishandle(101)
    clf(101)
end
hfig=figure(101);
    
for k=1:expected_clusters
    scatter3(normed_pseudotime(ClusterGroup2==k)',score3(ClusterGroup2==k,1),score3(ClusterGroup2==k,2),30,colorMARK2(k,:),'fill')
    hold on
end
legend(legendInfo,'Location', 'northeast')
xlabel('Cluster pseudotime')
ylabel('COMP1')
zlabel('COMP2')
grid on

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
p.EdgeLabel=h.Edges.Weight;
title('Lineage Progression')



end

