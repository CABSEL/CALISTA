function [ h] = AddNewEdge(h,nodes,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,az,el)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


        for k=1:expected_clusters
            scatter3(normed_pseudotime(ClusterGroup2==k)',score3(ClusterGroup2==k,1),score3(ClusterGroup2==k,2),30,colorMARK2(k,:),'fill')
            hold on
        end
        legend(legendInfo,'Location', 'northeast')
        xlabel('Cluster pseudotime')
        ylabel('COMP1')
        zlabel('COMP2')
        grid on
        
        EdgeColor=zeros(size(h.Edges,1)+1,3);
        h=addedge(h, nodes(1),nodes(2),nodes(3));
        check=sum(h.Edges{:,1}==[nodes(1) nodes(2)],2)==2;
        EdgeColor(check,:)=[1 0 0];
        p=[];
        p=plot(h,'EdgeLabel',h.Edges.Weight);
%         view(3)
        p.MarkerSize=1;
        % p.ArrowSize=15;
        p.EdgeAlpha=1;
        view([az,el]);
        p.NodeColor=colorMARK2;
        p.LineWidth=2.5;
        p.EdgeColor=EdgeColor;
        p.XData=x_center;
        p.YData=y_center;
        p.ZData=z_center;
        p.NodeLabel=[];
        p.EdgeLabel=abs(h.Edges.Weight);



end

