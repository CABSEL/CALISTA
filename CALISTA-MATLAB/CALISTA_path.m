function [Results]=CALISTA_path(Results,INPUTS)

% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. June 1, 2017.
if ~isfield(INPUTS,'moving_average_window')
    INPUTS.moving_average_window=10;
end

if ~isfield(INPUTS,'CV_plot')
    INPUTS.CV_plot=0;
end

if ~isfield(INPUTS,'plot_fig')
    INPUTS.plot_fig=0;
end

if ~isfield(INPUTS,'hclustering')
    INPUTS.hclustering=0;
end


%% Define cell paths
%%%%%%%%%%%%%%%%%%%%%%
if INPUTS.path_auto
      h=Results.TRANSITION.final_graph;
      h.Edges.Weight=abs(h.Edges.Weight);
    selected_edges=[];
    fprintf('\nAutomatically detect paths along the graph based on shortestpath..\n')
    ending_clusters=find(Results.cluster_progression==max(Results.cluster_progression));
    for i =1: length(ending_clusters)
        [shortest_path{i},~,edgepath] = shortestpath(h,1,ending_clusters(i));
        selected_edges=[selected_edges edgepath];
    end
    CELL_path=shortest_path;
    selected_edges=unique(selected_edges);
    %     highlight(p,'Edges',selected_edges)
    refined_nodes(:,1)=h.Edges.EndNodes(selected_edges,1);
    refined_nodes(:,2)=h.Edges.EndNodes(selected_edges,2);
    refined_nodes(:,3)=h.Edges.Weight(selected_edges);
    h_refined=[];
    h_refined=graph;
    non=length(unique(refined_nodes(:,1:2)));
    h_refined=addnode(h_refined,non);
    %     NumberOfEdges=non-1;
    %     NumberOfEdges=1;
    h_refined=addedge(h_refined, refined_nodes(:,1),refined_nodes(:,2),abs(refined_nodes(:,3)));
    Results.PATH.h_refined=h_refined;
    Results.PATH.refined_nodes=refined_nodes;
    
else
    figure(102)
    add_paths=true;
    count=1;
    while add_paths
        
        fprintf('%s %4i','Path num: ',count')
        fprintf('\n ******************************************************** \n')
        fprintf('\nKey the clusters in the path based on the progression (e.g. [1 2 3 4]): ')
        CELL_path{count}=input('');
        fprintf('\nPress 1 to add another path, 0 otherwise: ');
        continue_add=input('');
        if continue_add
            figure(102)
            count=count+1;
        else
            add_paths=false;
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_paths=length(CELL_path);
Results.num_cells_each_cluster=hist(Results.final_groups,1:Results.expected_clusters);
for i=1:n_paths
    
    nodes_connections=[CELL_path{i}(1:end-1)' CELL_path{i}(2:end)'];
    idx_swap_nodes=find(nodes_connections(:,1)>nodes_connections(:,2)); % swap nodes if node1>node2
    nodes_connections(idx_swap_nodes,:)=[nodes_connections(idx_swap_nodes,2) nodes_connections(idx_swap_nodes,1)];
    [~,idx_actual_edges_in_path]=ismember(nodes_connections,Results.TRANSITION.final_graph.Edges.EndNodes,'rows');
    idx_cells_in_path=find((ismember(Results.final_groups,CELL_path{i}))==1)';
    
    % Define cells in path
    cells_in_path=Results.ORDERING.cell_ordering(idx_cells_in_path,1);
    % Define transition genes of the actual path
    actual_transition_genes=Results.GENES.final_transition_genes(idx_actual_edges_in_path);
    actual_transition_genes=unique(cat(2,actual_transition_genes{:}));
    [~,idx_actual_transition_genes]=ismember(actual_transition_genes,Results.GENES.tot_transition_genes);
    n_transition_genes=length(idx_actual_transition_genes);
    [idx_sorted_cells,idxSORTED_cells_in_path]=sort(cells_in_path);
    path_mRNA_all{i}=Results.GENES.mRNA_tot_transition_genes(idx_cells_in_path,:);
    path_mRNA_all{i}= path_mRNA_all{i}(idxSORTED_cells_in_path,:);
    radius_width=min(round(max(length(cells_in_path)*INPUTS.moving_average_window/100,10)),100);
    smoothExpr{i}=movmean(path_mRNA_all{i},radius_width,1,'Endpoints','discard');
    if INPUTS.CV_plot
        smoothExpr_std{i}=movstd(path_mRNA_all{i},radius_width,1,'Endpoints','discard');
        CV{i}=smoothExpr_std{i}./smoothExpr{i};
        Results.PATH.CV{i}=CV{i};
    end
    windowCenters{i}=1:size(smoothExpr{i},1);
    
    if INPUTS.plot_fig
        colors=[ 1 0 0; 0 0.6 0 ];
        % Plot one figure with all genes
        hfig=figure(5000);
        %             set(hfig,'Position', [50, 50, 800, 800]);
        p=numSubplots(length(Results.GENES.tot_transition_genes));
        % Define colormap
        colorMARK_smooth=[];
        cell_in_clusters=Results.num_cells_each_cluster(CELL_path{i});
        cell_in_clusters=[cell_in_clusters(1) round(cell_in_clusters(2:end-1)/2) cell_in_clusters(end)];
        counts_cells_in_path_edges=sum([cell_in_clusters(1:end-1); cell_in_clusters(2:end)]',2);
        for j=1:size(nodes_connections,1)
            % Mesh colormap based on edge
            color1=Results.colorMARK_calista(nodes_connections(j,1),:);
            color2=Results.colorMARK_calista(nodes_connections(j,2),:);
            R=[]; G=[]; B=[];
            R(:,1)=linspace(color1(1),color2(1),counts_cells_in_path_edges(j));
            G(:,1)=linspace(color1(2),color2(2),counts_cells_in_path_edges(j));
            B(:,1)=linspace(color1(3),color2(3),counts_cells_in_path_edges(j));
            temp=[R G B];
            colorMARK_smooth=[colorMARK_smooth; temp];
        end
        % down sampling
        idx_subsampling=round(linspace(1, length(cells_in_path), size(smoothExpr{i},1)));
        colorMARK_smooth=colorMARK_smooth(idx_subsampling,:);
        
        for index = 1:length(Results.GENES.tot_transition_genes)
            subplot(p(1),p(2),index)
            gg=gscatter(windowCenters{i}, smoothExpr{i}(:,index),1:length(windowCenters{i}), colorMARK_smooth,'o',5,'off');
            for n = 1:length(windowCenters{i})
                set(gg(n), 'MarkerFaceColor', colorMARK_smooth(n,:));
            end
            st=Results.GENES.tot_transition_genes{index};
            title(st)
            grid on
            hold on
            xlabel('Cell Ordering')
            ylabel('Mean Expr')
            
            %                 suptitle([' Hierarchical clustering for path num ' num2str(i)])
            %         ylim([-15 15])
        end
        
        %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             ffig=[7 13 14];
        %             colors=[ 1 0 0;0 0.6 0];
        %             figure(8000+index)
        %             for index = 1:3
        %                     subplot(1,3,index)
        %
        %                 gg=gscatter(windowCenters{i}, smoothExpr{i}(:,ffig(index)),1:length(windowCenters{i}),colorMARK_smooth,'o',5,'off');
        %                 for n = 1:length(windowCenters{i})
        %                     set(gg(n), 'MarkerFaceColor', colorMARK_smooth(n,:));
        %                 end
        %                 st=Results.GENES.tot_transition_genes{ffig(index)};
        %                 %     title(st)
        %                 grid on
        %                 hold on
        %                 xlabel('Cell Ordering')
        %                 ylabel('Mean Expr')
        %                 %     set(gca,'color','none')
        %                 %                     suptitle([' Hierarchical clustering for path num ' num2str(i)])
        % %                             xlim([0 2000])
        %             end
        %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if INPUTS.CV_plot
            figure(6000)
            for index = 1:length(Results.GENES.tot_transition_genes)
                subplot(p(1),p(2),index)
                gg=gscatter(windowCenters{i}, CV{i}(:,index),1:length(windowCenters{i}),colorMARK_smooth,'o',5,'off');
                hold on
                st=Results.GENES.tot_transition_genes{index};
                title(st)
                grid on
                hold on
                xlabel('Cell Ordering')
                ylabel('Mean Expr')
                
            end
        end
        
    end
    
    
    if INPUTS.hclustering
        
        fprintf('\nHierarchical clustering based on cell ordering...\n')
        color_hierarchical=Results.c(idx_cells_in_path,:);
        
        ColumnLabelsColorValue.Labels=cellstr(num2str(idx_sorted_cells));
        ColumnLabelsColorValue.Colors= num2cell(color_hierarchical(idxSORTED_cells_in_path,:),2);
        
        transition_expression= Results.GENES.mRNA_tot_transition_genes(idx_cells_in_path,idx_actual_transition_genes);
        transition_expression= transition_expression(idxSORTED_cells_in_path,:);
        transition_expression=zscore(transition_expression); % stardardization
        transition_expression=zscore(transition_expression');
        cgo = clustergram(transition_expression,'Cluster',1,'Standardize','none','RowPDist','spearman','Linkage','average','Dendrogram',1);
        %     set(cgo,)
        set(cgo,'RowLabels',actual_transition_genes)
        set(cgo,'ColumnLabels',cellstr(num2str(idx_sorted_cells)),'ColumnLabelsColor',ColumnLabelsColorValue,'LabelsWithMarkers',true);
        %     addXLabel(cgo, ' Cell Ordering ')
        %     addYLabel(cgo, ' Hierarchical clustering ')
        %             set(0,'ShowHiddenHandles','on')
        %             % Get all handles from root
        %             allhnds = get(0,'Children');
        %             % Find hearmap axis and change the font size
        %             h = findall(allhnds, 'Tag', 'HeatMapAxes');
        %             set(h, 'FontSize', 20)
        addTitle(cgo,[' Hierarchical clustering for path num ' num2str(i)])
        clustergrams{i}=cgo;
        
    end
    path_transition_genes{i}=actual_transition_genes;
    
end

if exist('clustergrams','var')
    Results.PATH.hclustering=clustergrams;
end
Results.PATH.CELL_path=CELL_path;
Results.PATH.smoothExpr=smoothExpr;
Results.PATH.windowCenters=windowCenters;
Results.PATH.path_transition_genes=path_transition_genes;
