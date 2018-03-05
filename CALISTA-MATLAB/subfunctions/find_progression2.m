function  [Results,DATA]= find_progression2(Results,DATA,method)
final_groups=Results.final_groups;
nvars=length(final_groups);
final_num_clusters=Results.expected_clusters;
cluster_progression=zeros(1,final_num_clusters);
timeline=DATA.timeline;



if length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0 % no time info
    
    if (~isfield(Results,'starting_method')) || (Results.starting_method~=1 && Results.starting_method~=2)
        fprintf('\nNo time info found. Please enter the starting cell or the marker gene whenever available\n')
        Results.starting_method=input('Press 1 to enter the starting cell, 2 to enter the marker gene, 0 otherwise: ');
    end
    
    switch Results.starting_method
        case 1
            if ~isfield(Results,'starting_cell')
                Results.starting_cell=input('\nKey the index of the starting cell (e.g. 1): ');
                Results.starting_cell=round(Results.starting_cell); % in case it is not integer, just round the input
            end
            if Results.starting_cell>DATA.nvars || Results.starting_cell<1
                error(' Invalid index. Please re-run the section and key a correct startint cell.');
                Results=rmfield(Results, 'starting_cell');
            else
                final_groups2=final_groups';
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
                end
            end
            case 0
                fprintf('\nSkip relabelling step\n')
                final_groups2=final_groups';
    end
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