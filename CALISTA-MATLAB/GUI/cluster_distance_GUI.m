function [Results]=cluster_distance_GUI(DATA,Results,handles_graphtextwarning,handles_axes1,handles_uipanel1)

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

% arrayfun(@cla,handles_axes1) %refresh
axes(handles_axes1)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance=(my_mean-diag(my_mean)*ones(1,Results.expected_clusters))./DATA.numGENES;
[sorted_dist idx_dist]=sort(distance(:),'descend');
I=zeros(Results.expected_clusters*Results.expected_clusters,3);
[I(:,2),I(:,1)] = ind2sub(Results.expected_clusters,idx_dist);
I(:,3)=sorted_dist;
nodes_n=I;
nodes_n(1:Results.expected_clusters,:)=[];
nodes_all=nodes_n;
nodes_all(:,1)=min(nodes_n(:,1:2),[],2);
nodes_all(:,2)=max(nodes_n(:,1:2),[],2);
[~, ind_dir ]=unique(nodes_all(:,1:2),'rows');
nodes=nodes_all(ind_dir,:);
nodes=sortrows(nodes,-3);
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
% Results.TRANSITION.nodes=nodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

az=-26;
el=8;
[h,i]=AddGraph(cluster_distance,Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista_transition,az,el);


MaxNumberOfEdges=size(cluster_distance,1);

connected=false;
while connected==false
    NumberOfConnectedNodes=length(dfsearch(h,1));
    
    if NumberOfConnectedNodes==size(h.Nodes,1)
        connected=true;
    else
        [MaxEdgesAdded]=CheckNumberOfEdges(i+1,MaxNumberOfEdges);
        if MaxEdgesAdded==false
            i=i+1;
            h=addedge(h, cluster_distance(i,1),cluster_distance(i,2),cluster_distance(i,3));
        end
    end
    for k=1:Results.expected_clusters
        scatter3(Results.cell_cluster_progression(Results.final_groups==k)',Results.score3(Results.final_groups==k,1),Results.score3(Results.final_groups==k,2),30,Results.colorMARK_calista(k,:),'fill')
        hold on
    end
    legend(Results.legendInfo_calista_transition,'Location', 'northeast')
    xlabel('Cluster pseudotime')
    ylabel('COMP1')
    zlabel('COMP2')
    grid on
    
    
    
    p=[];
    p=plot(h,'EdgeLabel',h.Edges.Weight);
%     view(3)
    p.MarkerSize=1;
    % p.ArrowSize=15;
    p.EdgeAlpha=1;
    p.NodeColor=Results.colorMARK_calista;
    p.LineWidth=2.5;
    p.EdgeColor=[0 0 0];
    p.XData=x_center;
    p.YData=y_center;
    p.ZData=z_center;
    p.NodeLabel=[];
    p.EdgeLabel=abs(h.Edges.Weight);
    set(handles_graphtextwarning,'Visible','on')  % GRAPH CONNECTED
end

for i=1:Results.expected_clusters
    text(x_center(i),y_center(i),z_center(i),num2str(i),'FontSize',50);
end
Results.TRANSITION.x_center=x_center;
Results.TRANSITION.y_center=y_center;
Results.TRANSITION.z_center=z_center;

%% Add check boxes

% set up model
jCBList = com.jidesoft.swing.CheckBoxList;
jCBList.setModel(javax.swing.DefaultListModel); % dynamic items
% jCBList.ClickInCheckBoxOnly = true; % it's default
% jCBList.setSelectionMode(0); % I need single selection

% Now display onscreen:
jScrollPane = com.mathworks.mwswing.MJScrollPane(jCBList);
% Set relative position
hParent = uicontainer('Parent',handles_uipanel1);
parentPixelPos = getpixelposition(hParent);
pos = [1,1,parentPixelPos(3),parentPixelPos(4)]; % fill the parent uicontainer completely
[~, hContainer] = javacomponent(jScrollPane, pos, hParent);
set(hContainer, 'Units', 'normalized');  % better behavior on resizing


% javacomponent(jScrollPane, [630 200 230 280], gcf);
% set(hButton,'units','norm', 'position',[0.2,0.3,0.1,0.05]);
% add some items and set check/selecttion
for k=1:size(nodes,1)
    list_box_text{k}=sprintf('%1i %s %1i %s %4.2f',cluster_distance(k,1),'-',cluster_distance(k,2), ' Cluster distance: ',cluster_distance(k,3));
    jCBList.getModel.addElement(list_box_text{k});
end

% Update some items' state programmatically

for k=1:size(h.Edges,1)
    jCBList.addCheckBoxListSelectedIndex(k-1);
end

Results.mean_prob_in_out_cluster=mean_prob_in_out_cluster;
Results.TRANSITION.final_graph=h;
Results.cluster_distance=cluster_distance;
Results.jCBList= jCBList;
Results.list_box_text=list_box_text;
end