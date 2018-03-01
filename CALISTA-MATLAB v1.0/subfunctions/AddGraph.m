function [h,NumberOfEdges,IDXCOLORMARK3]=AddGraph(nodes,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,az,el,varargin)
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
    NumberOfEdges=non-1;
    NumberOfEdges=1;
    
    for j=1:NumberOfEdges
        h=addedge(h, nodes(j,1),nodes(j,2),nodes(j,3));
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

p=plot(h,'EdgeLabel',h.Edges.Weight);
view([az,el])
p.MarkerSize=1;
% p.ArrowSize=15;
p.EdgeAlpha=1;
p.NodeColor=colorMARK2;
p.LineWidth=2.5;
% p.EdgeColor=[0 0 0];
p.XData=x_center;
p.YData=y_center;
p.ZData=z_center;
p.NodeLabel=[];
p.EdgeLabel=abs(h.Edges.Weight);
colorMARK3=bone(size(nodes,1));
[~,IDXCOLORMARK3]=ismember(h.Edges.Weight,nodes(:,3));
colorMARK3=colorMARK3(1:size(h.Edges,1),:);
% colorMARK3=flipud(colorMARK3);
p.EdgeColor=colorMARK3(IDXCOLORMARK3,:);

end

