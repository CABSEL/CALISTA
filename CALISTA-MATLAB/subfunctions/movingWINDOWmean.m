function [Results]=movingWINDOWmean(Results,DATA,radius_width,plot_fig,trans_genes)

fprintf('\nCalculating mean expression based on cell ordering...\n')
nodes_connection=[Results.TRANSITION.final_graph.Edges.EndNodes(:,1) Results.TRANSITION.final_graph.Edges.EndNodes(:,2)];
n_edges=size(nodes_connection,1);
t0=0;
t1=1;
idx_actual_edge=Results.ORDERING.idx_actual_edge;
for i=1:n_edges
    % Define cells in edge
    cells_in_edge=Results.ORDERING.cell_ordering(idx_actual_edge{i},1);
    cells_in_edge=(cells_in_edge-min(cells_in_edge))/(max(cells_in_edge)-min(cells_in_edge));
    
    % Define transition genes of the actual edge
    switch trans_genes
        case 1
            actual_transition_genes=Results.GENES.final_transition_genes{i};
        case 2
            actual_transition_genes=Results.GENES.tot_transition_genes;
        case 3
            actual_transition_genes=DATA.genes;
    end
    Results.GENES.actual_transition_genes{i}=actual_transition_genes;
    [~,idx_actual_transition_genes]=ismember(actual_transition_genes,DATA.genes);
    n_transition_genes=length(idx_actual_transition_genes);
    %     % Sliding windows
    %     nWindows = round(length(cells_in_edge)/3);
    %     centersEQUI=linspace(t0,t1,nWindows);
    %     for j=1:nWindows
    %         [~, index] = min(abs(cells_in_edge-centersEQUI(j)));
    %         windowCenters{i}(j) = cells_in_edge(index); % Finds first one only!
    %     end
    %     windowRadius = radius_width*(t1-t0);
    %     edge_mRNA_all=DATA.totDATA(idx_actual_edge{i},idx_actual_transition_genes);
    %     smoothExpr{i} = calculateWindowMoving_pseudotime(windowCenters{i}(j), windowRadius, cells_in_edge,cells_in_edge, edge_mRNA_all);
    %
    
    edge_mRNA_all{i}=DATA.totDATA(idx_actual_edge{i},idx_actual_transition_genes);
    [cells_in_edgeSORTED,idxSORTED_cells_in_edge]=sort(cells_in_edge);
    edge_mRNA_all{i}=edge_mRNA_all{i}(idxSORTED_cells_in_edge,:);
    %     int1=1;
    %     int2=radius_width+1;
    %     count=1;
    %     while(int2<=length(idx_actual_edge{i}))
    %         windowCenters{i}(count)=count;
    % %         windowCenters{i}(count)=mean(cells_in_edgeSORTED(int1:int2));
    %         smoothExpr{i}(count,:)=mean(edge_mRNA_all{i}(int1:int2,:));
    %         int1=int1+1;
    %         int2=int2+1;
    %         count=count+1;
    %     end
    
    smoothExpr{i}=movmean(edge_mRNA_all{i},radius_width,1,'Endpoints','discard');
    windowCenters{i}=1:size(smoothExpr{i},1);
    if plot_fig
        % Plot one figure with all genes
        hfig=figure(1000+i);
        set(hfig,'Position', [50, 50, 800, 800]);
        [p,n]=numSubplots(n_transition_genes);
        % Mesh colormap based on edge
        color1=Results.colorMARK_calista(nodes_connection(i,1),:);
        color2=Results.colorMARK_calista(nodes_connection(i,2),:);
        R=[]; G=[]; B=[];
        R(:,1)=linspace(color1(1),color2(1),length(windowCenters{i}));
        G(:,1)=linspace(color1(2),color2(2),length(windowCenters{i}));
        B(:,1)=linspace(color1(3),color2(3),length(windowCenters{i}));
        colorMARK_smooth{i}=[R G B];
        for index = 1:n_transition_genes
            subplot(p(1),p(2),index)
            gg=gscatter(windowCenters{i}, smoothExpr{i}(:,index),1:length(windowCenters{i}),colorMARK_smooth{i},'o',10,'off');
            for n = 1:length(windowCenters{i})
                set(gg(n), 'MarkerFaceColor', colorMARK_smooth{i}(n,:));
            end
            st=actual_transition_genes{index};
            title(st)
            grid on
            hold on
            xlabel('Cell Ordering')
            ylabel('Mean Expr')
            xlim([0 windowCenters{i}(end)])
            
        end
        suptitle([num2str(Results.TRANSITION.nodes_connection(i,1)) ' - ' num2str(Results.TRANSITION.nodes_connection(i,2))])
        pause(2)
    end
end

Results.ORDERING.transition_genes_Expr= smoothExpr;
Results.ORDERING.edge_mRNA_all=edge_mRNA_all;
pause(3)