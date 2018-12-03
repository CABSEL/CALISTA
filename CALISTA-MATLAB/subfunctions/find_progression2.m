function  [Results,DATA]= find_progression2(Results,DATA,method)
final_groups=Results.final_groups;
nvars=length(final_groups);
final_num_clusters=Results.expected_clusters;
cluster_progression=zeros(1,final_num_clusters);
timeline=DATA.timeline;



if length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0 % no time info
    
    if (~isfield(Results,'starting_method')) || (Results.starting_method~=1 && Results.starting_method~=2 && Results.starting_method~=3)
        fprintf('\nNo time info found. Please enter the starting cell or the marker gene whenever available\n')
        Results.starting_method=input('Press 1 to enter the starting cell, 2 to enter the marker gene, 3 to enter the starting cluster, 0 otherwise: ');
    end
    
    switch Results.starting_method
        case 1
            if ~isfield(Results,'starting_cell')
                Results.starting_cell=input('\nKey the index of the starting cell (e.g. 1): ');
                Results.starting_cell=round(Results.starting_cell); % in case it is not integer, just round the input
            end
            idx_abs_cell=find(DATA.cut_sort.abs_indices==Results.starting_cell);
            if isempty(idx_abs_cell)
                error(' Invalid index. Please re-run the section and key a correct startint cell.');
                Results=rmfield(Results, 'starting_cell');
            else
                new_starting_cluster=unique(final_groups(Results.starting_cell));
                final_groups2=final_groups';
                Results.starting_cell=idx_abs_cell;
                if final_groups(Results.starting_cell)~=1 %define the cluster with the starting cell as cluster 1
                    final_groups2(find(final_groups==final_groups(Results.starting_cell)))=1;
                    final_groups2(find(final_groups==1))=final_groups(Results.starting_cell);
                end
            end
        case 2
            if ~isfield(Results,'markerGENE') % up- or down- regulated in all differentiation paths
                Results.markerGENE=input('\nEnter the marker gene chosen (e.g. ''Erg'', case insensitive): ');
            end
            if ~isfield(Results,'UPorDOWN')
                Results.UPorDOWN=input('\nPress 1 if it is upregulated, 2 if it is downregulated (e.g. 2): ');
            end
            [tf, idx_markerGENE]=ismember(lower(Results.markerGENE), lower(DATA.genes));
            if tf==0
                error(' Marker gene not found. Please re-run the section and enter a correct gene''s name');
                Results=rmfield(Results, 'markerGENE');
            else
                for i=1:Results.expected_clusters
                    mean_markerGENE(i)=mean(DATA.totDATA(find(Results.final_groups==i),idx_markerGENE));
                end
                
                switch Results.UPorDOWN
                    case 1
                        [~,initial_cluster]=min(mean_markerGENE);
                    case 2
                        [~,initial_cluster]=max(mean_markerGENE);
                    otherwise
                        error('Please re-run the section and enter 1 for up- or 2 for down- regulated marked gene');
                        Results=rmfield(Results, 'UPorDOWN');
                end
                final_groups2=final_groups';
                if initial_cluster~=1 %define the cluster with the starting cell as cluster 1
                    final_groups2(find(final_groups==initial_cluster))=1;
                    final_groups2(find(final_groups==1))=initial_cluster;
                    new_starting_cluster=initial_cluster;
                end
            end
        case 3
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            totDATA=DATA.totDATA.^(1/2);
            zeros_cells=sum(totDATA==0,2)*100/size(totDATA,2);
            dot_size=30*ones(size(totDATA,1),1);
            cluster_predicted=unique(final_groups);
            ClusterGroup_calista = final_groups;
            colorMARK_calista=jet(length(cluster_predicted));
            % ********** PCA **********
            [~, score3] = pca(zscore(totDATA));
            figure_temp=figure(10);
            subplot(121)
            for i=1:length(cluster_predicted)
                temp=dot_size(ClusterGroup_calista==cluster_predicted(i));
                scatter3(score3(ClusterGroup_calista==cluster_predicted(i),1), score3(ClusterGroup_calista==cluster_predicted(i),2), score3(ClusterGroup_calista==cluster_predicted(i),3),temp, repmat(colorMARK_calista(i,:),size(temp,1),1), 'fill');
                title('Cell Clustering')
                xlabel('PC1')
                ylabel('PC2')
                zlabel('PC3')
                grid on
                hold on
                legendInfo_calista{i} = sprintf( '%s %4i', 'Cluster ', i);
            end
            
            legend(legendInfo_calista,'Location', 'northeast')
            %Plot with true labels
            fprintf('\nPlease select file with true labels/cell ID ')
            [FileName,PathName,FilterIndex] = uigetfile('*.*');
            filename=strcat(PathName, FileName);
            cell_labels=importdata(filename);
            if ~isempty(DATA.cut_sort.cell_2_cut)
                cell_labels(DATA.cut_sort.cell_2_cut)=[]; % remove labels for removed clusters
            end
            
            if ~isempty(DATA.cut_sort.idx2cutCELL)
                cell_labels(DATA.cut_sort.idx2cutCELL)=[]; % remove labels for removed cells based on %zeros
            end
            
            if length(DATA.cut_sort.idx_sorted_cells)>DATA.nvars
            else
                cell_labels=cell_labels(DATA.cut_sort.idx_sorted_cells);
            end
            [unique_cell_labels,IA,IC]=unique(cell_labels,'rows','stable');
            
            colorMARK_labels=parula(length(unique_cell_labels));
            
            subplot(122)
            for i=1:length(unique_cell_labels)
                scatter3(score3(IC==i,1), score3(IC==i,2),score3(IC==i,3),dot_size(IC==i),colorMARK_labels(i,:),'fill')
                xlabel('PC1')
                ylabel('PC2')
                zlabel('PC3')
                title('True labels')
                grid on
                hold on
            end
            
            if iscell(unique_cell_labels)
                legend(cellstr(unique_cell_labels),'Location', 'northeast')
            else
                legend(num2str(unique_cell_labels),'Location', 'northeast')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ~isfield(Results,'starting_cluster')
                Results.starting_cluster=input('\nKey the index of the starting cluster (e.g. 2): ');
                Results.starting_cluster=round(Results.starting_cluster); % in case it is not integer, just round the input
            end
            
            idx_cells_in_starting_cluster=find(Results.final_groups==Results.starting_cluster);
            if isempty(idx_cells_in_starting_cluster)
                error(' Invalid index. Please re-run the section and key a correct startint cluster.');
                Results=rmfield(Results, 'starting_cluster');
            else
                final_groups2=final_groups';
                if Results.starting_cluster~=1 %define the cluster with the starting cell as cluster 1
                    new_starting_cluster=Results.starting_cluster;
                    final_groups2(idx_cells_in_starting_cluster)=1;
                    final_groups2(find(final_groups==1))=Results.starting_cluster;
                end
            end
            close(figure_temp)
        case 0
            fprintf('\nSkip relabelling step\n')
            final_groups2=final_groups';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update all variables
    
    stages=1:final_num_clusters;
    stt1='Cluster  ';
    stt2='pseudo-stage ';
    for clust=1:final_num_clusters
        idx_final_groups=find(final_groups2==clust);
        final_groups_progression(idx_final_groups)=stages(clust);
        legendInfo{clust} = sprintf( '%s %4i', stt1, clust);
    end
    
else
    for clust=1:final_num_clusters
        idx_final_groups=find(final_groups==clust);
        switch method
            case 1
                cluster_progression(clust)=median(timeline(idx_final_groups));
            case 2
                cluster_progression(clust)=mean(timeline(idx_final_groups));
            case 3
                cluster_progression(clust)=mode(timeline(idx_final_groups));
        end
        
        MEAN_progression(clust)=mean(timeline(idx_final_groups));
        SD_progression(clust)=std(timeline(idx_final_groups));
    end
    [stages,idx_progression]=sort(cluster_progression);
    final_groups2=zeros(nvars,1);
    stt1='Cluster  ';
    stt2='pseudo-stage ';
    for clust=1:final_num_clusters
        new_clust_order=find(idx_progression==clust);
        idx_final_groups=find(final_groups==clust);
        final_groups2(idx_final_groups)=new_clust_order;
        final_groups_progression(idx_final_groups)=stages(clust);
        legendInfo{clust} = sprintf( '%s %4i', stt1, clust);
    end
end
Results.final_groups=final_groups2';
Results.cluster_progression=(stages-min(stages))/(max(stages)-min(stages)); %normalized
Results.legendInfo_calista=legendInfo;
Results.final_groups_progression=final_groups_progression;
Results.cluster_predicted=unique(Results.final_groups);
DATA.my_stages=final_groups_progression';

Results.cell_cluster_progression=zeros(length(Results.final_groups),1);
for i=1:Results.expected_clusters
    Results.cell_cluster_progression(find(Results.final_groups==Results.cluster_predicted(i)))= Results.cluster_progression(i);
end
% Results.cell_cluster_progression=(Results.cell_cluster_progression-min(Results.cell_cluster_progression))/max(Results.cell_cluster_progression-min(Results.cell_cluster_progression));

end