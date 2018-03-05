function  h=temp_RmEdge(h,nodes,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,az,el)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for k=1:expected_clusters
    scatter3(normed_pseudotime(ClusterGroup2==k)',score3(ClusterGroup2==k,1),score3(ClusterGroup2==k,2),30,colorMARK2(k,:),'fill')
    hold on
end
legend(legendInfo,'Location', 'northeast')
xlabel('Cluster progression')
ylabel('COMP1')
zlabel('COMP2')
grid on

h=rmedge(h, nodes(1),nodes(2));
EdgeColor=zeros(size(h.Edges,1),3);
% [~, Ind]=min(h.Edges{:,2});
% if ~isempty(Ind)
%     EdgeColor(Ind,:)=[1 0 0];
% end

p=[];
p=plot(h,'EdgeLabel',h.Edges.Weight);
view([az,el])
p.MarkerSize=1;
% p.ArrowSize=15;
p.EdgeAlpha=1;
p.NodeColor=colorMARK2;
p.LineWidth=2.5;
p.EdgeColor=EdgeColor;
p.XData=x_center;
p.YData=y_center;
p.ZData=z_center;
p.NodeLabel=[];

end