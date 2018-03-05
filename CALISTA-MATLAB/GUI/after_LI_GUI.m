function Results= after_LI_GUI(DATA,Results)


hh=Results.TRANSITION.final_graph;
nodes_connection3=[hh.Edges.EndNodes(:,1) hh.Edges.EndNodes(:,2)];
nodes=[nodes_connection3 hh.Edges.Weight];
x_center=Results.TRANSITION.x_center;
y_center=Results.TRANSITION.y_center;
z_center=Results.TRANSITION.z_center;
figure(102)
az=-37.5000;
el=30;
[hh,i]=AddGraph(nodes,Results.colorMARK_calista,x_center,y_center,z_center,Results.final_groups,Results.expected_clusters,Results.cell_cluster_progression,Results.score3,Results.legendInfo_calista_transition,az,el,hh);
title('Lineage Progression')


mean_prob_in_out_cluster=Results.mean_prob_in_out_cluster;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% If there's no time info, reorder clusters based on the connections
% h=hh;

if length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0
    Results.TRANSITION.clusterTEXT=1;
    UG=zeros(Results.expected_clusters);
    for k=1:size(hh.Edges,1)
        UG(hh.Edges.EndNodes(k,1),hh.Edges.EndNodes(k,2))=hh.Edges.Weight(k);
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
    hh = reordernodes(hh,order);
    
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
    
    p=plot(hh,'EdgeLabel',hh.Edges.Weight);
    view([az,el])
    p.MarkerSize=1;
    % p.ArrowSize=15;
    p.EdgeAlpha=1;
    p.NodeColor=Results.colorMARK_calista;
    p.LineWidth=2.5;
%     % p.EdgeColor=[0 0 0];
    p.XData=x_center;
    p.YData=y_center;
    p.ZData=z_center;
    p.NodeLabel=[];
    p.EdgeLabel=abs(hh.Edges.Weight);
%     colorMARK3=bone(size(nodes,1));
%     [~,IDXCOLORMARK3]=ismember(hh.Edges.Weight,nodes(:,3));
%     colorMARK3=colorMARK3(1:size(hh.Edges,1),:);
    % colorMARK3=flipud(colorMARK3);
%     p.EdgeColor=colorMARK3(IDXCOLORMARK3,:);
    
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
   pause(1)
end

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

% Results.TRANSITION.cluster_distance=cluster_distance;
pause(2)