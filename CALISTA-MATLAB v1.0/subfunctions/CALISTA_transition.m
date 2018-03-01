function [Results]=CALISTA_transition(DATA,Results)

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
if ishandle(101)
clf(101)
end
hfig=figure(101);
set(hfig,'Position',[500, 500, 500, 500])
az=-37.5000;
el=30;

[h,i]=AddGraph(nodes,Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista_transition,az,el);

MaxNumberOfEdges=size(nodes,1);

connected=false;
while connected==false
    NumberOfConnectedNodes=length(dfsearch(h,1));
    
    if NumberOfConnectedNodes==size(h.Nodes,1)
        connected=true;
    else
        [MaxEdgesAdded]=CheckNumberOfEdges(i+1,MaxNumberOfEdges);
        if MaxEdgesAdded==false
            i=i+1;
            h=addedge(h, nodes(i,1),nodes(i,2),nodes(i,3));
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
    view(3)
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
    title('Lineage Progression')
    colorMARK3=bone(size(nodes,1));
    [~,IDXCOLORMARK3]=ismember(h.Edges.Weight,nodes(:,3));
%     colorMARK3=colorMARK3(1:size(hh.Edges,1),:);
    % colorMARK3=flipud(colorMARK3);
    p.EdgeColor=colorMARK3(IDXCOLORMARK3,:);
    
end





stop_adding=false;
add_edge=true;
MaxNumberOfEdges=size(nodes,1);

while stop_adding==false
    
    NumberOfConnectedNodes=length(dfsearch(h,1));
    fprintf('%d edge(s) have been added and the graph is connected. \nIf you want to add another edge press "p" \nIf you want to remove an edge press "m" \nIf you want to continue with the next step press "enter" \n(Please make sure that figure 101 is in the foreground and no additional tools are selected (e.g. zooming, rotation, data cursor etc.)) \n\n----------------------------------------------\n\n',i);
        
    t=getkey;
    true_key=false;
    % 13->enter
    while true_key==false
        if t==13 %&& NumberOfConnectedNodes==size(h.Nodes,1)
            if NumberOfConnectedNodes ~=size(h.Nodes,1)
                warn_diag1=warndlg('The graph must be CONNECTED ','!! Warning !!')
                set(warn_diag1,'paperposition',[0.1,0,50,50])
                %add_edge=[]
                break;
            end
            close(101);
            stop_adding=true;
            figure(102)
            [h,i]=AddGraph(nodes,Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista_transition,az,el,h);
            title('Lineage Progression')
        end
        %43 -> +
        if t==80 || t==112
            add_edge=true;
            %45 -> -
        elseif t==77 || t==109
            add_edge=false;
        end
        if ismember(t, [13 80 112 77 109])
            true_key=true;
        else
            warn_diag2=warndlg('Please press ''p'',''m'',or ''Enter'' ','!! Warning !!')
            set(warn_diag2,'paperposition',[0,0,50,50])
            t=getkey
        end
        
    end
    
    if add_edge==true && stop_adding==false && true_key==true
        figure(101)
        [az, el]=view;
        hold off
        [MaxEdgesAdded]=CheckNumberOfEdges(i+1,MaxNumberOfEdges);
        
        if MaxEdgesAdded==false
            i=i+1;
            h=AddNewEdge(h,nodes(i,:),Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista_transition,az,el);
            title('Lineage Progression')
        else
            fprintf('You have already added the maximum number of edges \n')
            
        end
    end
    
    if add_edge==false & stop_adding==false & true_key==true
        h_temp=temp_RmEdge(h,nodes(i,:),Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista_transition,az,el);
        NumberOfConnectedNodes2=length(dfsearch(h_temp,1));
        if NumberOfConnectedNodes2==size(h_temp.Nodes,1)
            figure(101)
            [az, el]=view;
            hold off
            h=RmEdge(h,nodes(i,:),Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista_transition,az,el);
            i=i-1;
            title('Lineage Progression')
        else
            warn_diag1=warndlg('The edge can not be removed. The graph must be CONNECTED ','!! Warning !!')
            set(warn_diag1,'paperposition',[0.1,0,50,50])
            
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there's no time info, reorder clusters based on the connections
% h=hh;
hh=h;
if length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0
    Results.TRANSITION.clusterTEXT=1;
    UG=zeros(Results.expected_clusters);
    for k=1:size(h.Edges,1)
        UG(h.Edges.EndNodes(k,1),h.Edges.EndNodes(k,2))=h.Edges.Weight(k);
    end
    UG=abs(tril(UG + UG'));
    UG=sparse(UG);
    for i=2:Results.expected_clusters
        [dist(i),path{i}] = graphshortestpath(UG,1,i,'directed',false);
    end
    
    [cluster_progression,bb]=sort(dist);
    final_groups2=Results.final_groups;
    c2=Results.c;
    order=1;
    for i=2:Results.expected_clusters
        order(i)=path{bb(i)}(end);
        find_temp_idx=find(Results.final_groups==order(i));
        final_groups2(find_temp_idx)=Results.cluster_predicted(i);
        c2(find_temp_idx,:)=repmat(Results.colorMARK_calista(i,:),length(find_temp_idx),1);
        
    end
    h = reordernodes(h,order);
    
    Results.final_groups=final_groups2;
    Results.c=c2;
    Results.cluster_progression=(cluster_progression-min(cluster_progression))/(max(cluster_progression)-min(cluster_progression));
    
    Results.cell_cluster_progression=zeros(length(Results.final_groups),1);
    for i=1:Results.expected_clusters
        Results.cell_cluster_progression(find(Results.final_groups==Results.cluster_predicted(i)))= Results.cluster_progression(i);
    end
%     Results.cell_cluster_progression=(Results.cell_cluster_progression-min(Results.cell_cluster_progression))/max(Results.cell_cluster_progression-min(Results.cell_cluster_progression));
    
    
    for i=1:Results.expected_clusters
        x_center(i)=mean(Results.cell_cluster_progression(Results.final_groups==i));
        y_center(i)=mean(Results.score3(Results.final_groups==i,1));
        z_center(i)=mean(Results.score3(Results.final_groups==i,2));
        
    end
    
    close(102)
    figure(102)
    for k=1:Results.expected_clusters
        scatter3(Results.cell_cluster_progression( Results.final_groups==k)',Results.score3( Results.final_groups==k,1),Results.score3( Results.final_groups==k,2),30,Results.colorMARK_calista(k,:),'fill')
        hold on
    end
    
    stt1='Cluster:';
    stt2='Cluster pseudotime:';
    for clust=1:Results.expected_clusters
        Results.legendInfo_calista_transition{clust} = sprintf( '%s %4i %s %1.2f', stt1, clust,stt2,Results.cluster_progression(clust));
    end
    
    legend(Results.legendInfo_calista_transition,'Location', 'northeast')
    xlabel('Cluster pseudotime')
    ylabel('COMP1')
    zlabel('COMP2')
    grid on
    
    p=plot(h,'EdgeLabel',h.Edges.Weight);
    view([az,el])
    p.MarkerSize=1;
    % p.ArrowSize=15;
    p.EdgeAlpha=1;
    p.NodeColor=Results.colorMARK_calista;
    p.LineWidth=2.5;
    % p.EdgeColor=[0 0 0];
    p.XData=x_center;
    p.YData=y_center;
    p.ZData=z_center;
    p.NodeLabel=[];
    p.EdgeLabel=abs(hh.Edges.Weight);
    colorMARK3=bone(size(nodes,1));
    [~,IDXCOLORMARK3]=ismember(hh.Edges.Weight,nodes(:,3));
    colorMARK3=colorMARK3(1:size(hh.Edges,1),:);
    % colorMARK3=flipud(colorMARK3);
    p.EdgeColor=colorMARK3(IDXCOLORMARK3,:);
    
    title(' Plot after cluster relabelling')
    for i=1:Results.expected_clusters
        text(x_center(i),y_center(i),z_center(i),num2str(i),'FontSize',50);
    end
    
    % Replot clustering colors
    figure(1000)
    cla(Results.h2)
    subplot(122);
    for k=1:Results.expected_clusters
        scatter3(Results.score3( Results.final_groups==k,1),Results.score3( Results.final_groups==k,2),Results.score3( Results.final_groups==k,3),30,Results.colorMARK_calista(k,:),'fill')
        hold on
        title('Cell Clustering after relabelling')
        xlabel('PC1')
        ylabel('PC2')
        zlabel('PC3')
        grid on
       %     legendInfo_calista{i} = sprintf( '%s %4i', 'Cluster ', i);
    end
    legend(Results.legendInfo_calista,'Location', 'northeast')
    pause(3)
end

% Remove unwanted edges before additional analysis

hh=h;
proceed=input('Press 1 if you want to remove edges, 0 otherwise: ');
if ~isempty(proceed) & proceed==1
    Results.TRANSITION.clusterTEXT=1;
    figure(102)
    for i=1:Results.expected_clusters
        text(x_center(i),y_center(i),z_center(i),num2str(i),'FontSize',50);
    end
    
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
    
    close(102)
    hfig=figure(102);
    set(hfig,'position', [500, 500, 500, 500])
    az=-37.5000;
    el=30;
    for k=1:Results.expected_clusters
        scatter3(Results.cell_cluster_progression(Results.final_groups==k)',Results.score3(Results.final_groups==k,1),Results.score3(Results.final_groups==k,2),30,Results.colorMARK_calista(k,:),'fill')
        hold on
    end
    grid on
    p=plot(hh);
    view([az,el])
    p.MarkerSize=15;
    % p.ArrowSize=15;
    p.EdgeAlpha=1;
    p.NodeColor=Results.colorMARK_calista;
    p.LineWidth=2.5;
    % p.EdgeColor=EdgeColor;
    p.XData=x_center;
    p.YData=y_center;
    p.ZData=z_center;
    p.NodeLabel=[];
    p.EdgeLabel=abs(hh.Edges.Weight);
    colorMARK3=bone(size(nodes,1));
    [~,IDXCOLORMARK3]=ismember(hh.Edges.Weight,nodes(:,3));
%     colorMARK3=colorMARK3(1:size(hh.Edges,1),:);
    % colorMARK3=flipud(colorMARK3);
    p.EdgeColor=colorMARK3(IDXCOLORMARK3,:);
    
    legend(Results.legendInfo_calista_transition,'Location', 'northeast')
    
    title('Lineage Progression')
    xlabel('Cluster pseudotime')
    ylabel('PC1')
    zlabel('PC2')
    for i=1:Results.expected_clusters
        text(x_center(i),y_center(i),z_center(i),num2str(i),'FontSize',50);
    end
    pause(3)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results.TRANSITION.x_center=x_center;
Results.TRANSITION.y_center=y_center;
Results.TRANSITION.z_center=z_center;
Results.TRANSITION.final_graph=hh;

if exist('order','var')
    Results.mean_prob_in_out_cluster=mean_prob_in_out_cluster(order,:);
end
Results.mean_prob_in_out_cluster=mean_prob_in_out_cluster;


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
        title(sprintf('CALISTA cluster pseudotime %1.2f',stages_name(K)/max(Results.cluster_progression)))
    
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
hh=Results.TRANSITION.final_graph;
nodes_connection3=[hh.Edges.EndNodes(:,1) hh.Edges.EndNodes(:,2)];
Results.TRANSITION.nodes_connection=nodes_connection3;

Results.TRANSITION.cluster_distance=cluster_distance;
pause(2)
end