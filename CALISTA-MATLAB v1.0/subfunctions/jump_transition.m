function Results=jump_transition(DATA,Results)


ClusterGroup2=Results.final_groups;
expected_clusters=Results.expected_clusters;
normed_pseudotime=Results.cell_cluster_progression;
score3=Results.score3;
legendInfo=Results.legendInfo_calista;
colorMARK2=Results.colorMARK_calista;
x_center=zeros(Results.expected_clusters,1);
z_center=x_center;
y_center=x_center;
for i=1:Results.expected_clusters
    x_center(i)=mean(Results.cell_cluster_progression(Results.final_groups==i));
    y_center(i)=mean(Results.score3(Results.final_groups==i,1));
    z_center(i)=mean(Results.score3(Results.final_groups==i,2));
    
end

if ishandle(101)
    clf(101)
end
hfig=figure(101);
set(hfig,'Position',[500, 500, 500, 500])
az=-37.5000;
el=30;


for k=1:expected_clusters
    scatter3(normed_pseudotime(ClusterGroup2==k)',score3(ClusterGroup2==k,1),score3(ClusterGroup2==k,2),30,colorMARK2(k,:),'fill')
    hold on
end
legend(legendInfo,'Location', 'northeast')
xlabel('Cluster pseudotime')
ylabel('COMP1')
zlabel('COMP2')
grid on

for i=1:Results.expected_clusters
    text(x_center(i),y_center(i),z_center(i),num2str(i),'FontSize',50);
end


% Select the lineage manually
add_paths=true;
count=1;
tot_nodes=[];
while add_paths
    
    fprintf('%s %4i','Path num: ',count')
    fprintf('\n ******************************************************** \n')
    CELL_path{count}=input(' Key the clusters in the path based on the progression (e.g. [1 2 3 4]: ');
    nodes=[CELL_path{count}(1:end-1)' CELL_path{count}(2:end)'];
    tot_nodes=[tot_nodes; nodes];
    continue_add=input(' Press 1 to add another path, 0 otherwise: ');
    if continue_add
        
        count=count+1;
    else
        add_paths=false;
    end
end


tot_nodes=unique(tot_nodes,'rows');
h=graph;
non=length(unique(tot_nodes(:,1:2)));
h=addnode(h,non); 

% Construct the graph
NumberOfEdges=non-1;
for j=1:NumberOfEdges
    h=addedge(h, tot_nodes(j,1),tot_nodes(j,2),1); %same weight for all edges
end

% Check to have a connected graph
NumberOfConnectedNodes=length(dfsearch(h,1));
if NumberOfConnectedNodes ~=size(h.Nodes,1)
    error('The graph must be CONNECTED. Please run CALISTA again')
end

% Plot the graph
p=plot(h);
p.MarkerSize=15;
% p.ArrowSize=15;
p.EdgeAlpha=1;
p.NodeColor=colorMARK2;
p.LineWidth=2.5;
p.EdgeColor=[0 0 0];
p.XData=x_center;
p.YData=y_center;
p.ZData=z_center;
p.NodeLabel=[];
p.EdgeLabel=[];
title('Cluster ordering')

Results.PATH.CELL_path=CELL_path;
Results.TRANSITION.x_center=x_center;
Results.TRANSITION.y_center=y_center;
Results.TRANSITION.z_center=z_center;
Results.TRANSITION.final_graph=h;
hh=Results.TRANSITION.final_graph;
nodes_connection3=[hh.Edges.EndNodes(:,1) hh.Edges.EndNodes(:,2)];
Results.TRANSITION.nodes_connection=nodes_connection3;


%% Plot each predicted cluster separately
p=numSubplots(length(unique(Results.cluster_progression)));
stages_name=unique(Results.cluster_progression);
figure('position', [500, 500, 800, 800])
for K=1:length(stages_name)
    subplot(p(1),p(2),K)
    
    actual_stage=find(Results.cluster_progression==stages_name(K));
    for ii=1:length(actual_stage)
        temp_plot=find(Results.final_groups==actual_stage(ii));
        scatter3(Results.score3(temp_plot,1), Results.score3(temp_plot,2), Results.score3(temp_plot,3),30, Results.colorMARK_calista(actual_stage(ii),:), 'fill');
        hold on
    end
    grid on
    
    [~,idx_remove]=ismember(actual_stage,Results.cluster_predicted);
    
    cluster_plot=Results.cluster_predicted;
    cluster_plot(idx_remove)=[];
    for iii=1:length(cluster_plot)
        temp_plot=find(Results.final_groups==cluster_plot(iii));
        scatter3(Results.score3(temp_plot,1), Results.score3(temp_plot,2), Results.score3(temp_plot,3),30, [0.8 0.8 0.8], 'fill');
        title(sprintf('CALISTA cluster pseudotime %g',stages_name(K)/max(Results.cluster_progression)))
    
        xlabel('PC1')
        ylabel('PC2')
        zlabel('PC3')
    end
    
    
    hold off
end

Results.stages_name=stages_name;
%% Pie charts

if ~(length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0) || isfield(Results,'cell_labels')
    p=numSubplots(Results.expected_clusters);
    figure
    for i=1:Results.expected_clusters
        idx_cells_in_k=find( Results.final_groups==Results.cluster_predicted(i));
        if isfield(Results,'cell_labels')
            labels_in_k=Results.cell_labels(idx_cells_in_k);
            [aa,~,cc]=unique(labels_in_k);
            occurrence = (hist(cc,length(aa)))*100/length(labels_in_k);
            cut=find(occurrence<=5);
            occurrence(cut)=[];
            aa(cut)=[];
            
            ax=subplot(p(1),p(2),i);
            [~,ee]=ismember(aa,Results.unique_cell_labels);
            pp=pie(ax,occurrence,strcat(aa));
            hp = findobj(pp, 'Type', 'patch');
            for j=1:length(occurrence)
                
                set(hp(j), 'FaceColor', Results.colorMARK_labels(ee(j),:));
            end
            
            title(ax,sprintf( '%s %4i', 'Cluster', Results.cluster_predicted(i)))
            
        else
            
            labels_in_k=DATA.timeline(idx_cells_in_k);
            [aa,~,cc]=unique(labels_in_k);
            occurrence = (hist(cc,length(aa)))*100/length(labels_in_k);
            
            cut=find(occurrence<=5);
            occurrence(cut)=[];
            aa(cut)=[];
            
            ax=subplot(p(1),p(2),i);
            [~,ee]=ismember(aa,DATA.time);
            pp=pie(ax,occurrence,cellstr(strcat('t ',num2str(aa))));
            hp = findobj(pp, 'Type', 'patch');
            for j=1:length(occurrence)
                
                set(hp(j), 'FaceColor', Results.colorMARK_time(ee(j),:));
            end
            
            title(ax,sprintf( '%s %4i', 'Cluster', Results.cluster_predicted(i)))
            
        end
    end
end
