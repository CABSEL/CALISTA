function [Results]=visualization(reduction,DATA,Results)

fprintf('\nPlotting...\n')
expected_clusters=Results.expected_clusters;
final_groups=Results.final_groups;
% totDATA=log2(mRNA_all+1);
totDATA=DATA.totDATA.^(1/2);
switch reduction
    case 1
        % ********** t-SNE **********
        %set parameters
        initial_dims=DATA.numGENES;
        perplexity=30; %usually between 5-50
        no_dims=3;
        % all the data
        dataANALYSING=totDATA;
        score3= tsne2(dataANALYSING, [], no_dims, initial_dims, perplexity);
        
    case 2
        % ********** PCA **********
        [coeff, score3] = pca(zscore(totDATA));
        Results.coeff=coeff;
    case 3
        % ********** DIFFUSION MAP **********
        %without zscore
        no_dims = 4;
        t       = 10;
        S   = 50;
        [score3] = diffusion_maps(zscore(totDATA), no_dims, t, S); %or UPandDOWNdata
end

% *** Colormap for the original cell info ***
colorMARK_time=cool(DATA.nc);

% *** Colormap after clustering ***
% final_groups=U_star_max;
cluster_predicted=Results.cluster_predicted;
ClusterGroup_calista = final_groups;
% colorMARK1=hsv(length(cluster_predicted));
colorMARK_calista=jet(length(cluster_predicted));
c=zeros(DATA.nvars,3);
for i=1:expected_clusters
    idx_temp=find(ClusterGroup_calista==Results.cluster_predicted(i));
    c(idx_temp,:)=repmat(colorMARK_calista(Results.cluster_predicted(i),:),length(idx_temp),1);
end

% Plot
hfig=figure(1000);
set(hfig,'position', [500, 500, 1200, 400]) 
subplot(121)
for i=1:DATA.num_time_points
    scatter3(score3(DATA.timeline==DATA.time(i),1), score3(DATA.timeline==DATA.time(i),2),score3(DATA.timeline==DATA.time(i),3),30,colorMARK_time(i,:),'fill')
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    title('Original time/cell stage info')
    grid on
    hold on
    
    legendInfo_time{i} = sprintf( '%s %4i', 'Time/Stage', DATA.time(i));
end
legend(legendInfo_time,'Location', 'northeast')

Results.h2=subplot(122);
for i=1:length(cluster_predicted)
    scatter3(score3(ClusterGroup_calista==cluster_predicted(i),1), score3(ClusterGroup_calista==cluster_predicted(i),2), score3(ClusterGroup_calista==cluster_predicted(i),3),30, colorMARK_calista(i,:), 'fill');
    title('Cell Clustering')
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    grid on
    hold on
%     legendInfo_calista{i} = sprintf( '%s %4i', 'Cluster ', i);
end
legend(Results.legendInfo_calista,'Location', 'northeast')
pause(3)

Results.colorMARK_calista=colorMARK_calista;
Results.colorMARK_time=colorMARK_time;
Results.score3=score3;
Results.c=c;


% %% Plot with true labels
% true_labels=input('\nKey 1 for true labels visualization, 0 otherwise :');
% 
% if true_labels
%     [FileName,PathName,FilterIndex] = uigetfile('*.*');
%     filename=strcat(PathName, FileName);
%     cell_labels=importdata(filename);
%     if ~isempty(DATA.cut_sort.idx2cutCELL)
%         cell_labels(DATA.cut_sort.idx2cutCELL)=[]; % remove labels for removed cells
%     end
%     
%     if length(DATA.cut_sort.idx_sorted_cells)>DATA.nvars
%     else
%     cell_labels=cell_labels(DATA.cut_sort.idx_sorted_cells);
%     end
%     [unique_cell_labels,IA,IC]=unique(cell_labels,'rows','stable');
%     
% %     colorMARK_labels=[42 183 127;
% %         51 153 255;
% %         204 0 204;
% %                 252 144 201;
% %         153 153 0]/255;
%      colorMARK_labels=parula(length(unique_cell_labels));
%     c_true_labels=colorMARK_labels(IC,:);%zeros(length(cell_labels),3);
%     
%     % Plot
%     figure('position', [500, 500, 400, 400]) 
%     for i=1:length(unique_cell_labels)
%         scatter3(score3(IC==i,1), score3(IC==i,2),score3(IC==i,3),30,colorMARK_labels(i,:),'fill')
%         xlabel('PC1')
%         ylabel('PC2')
%         zlabel('PC3')
%         title('True labels')
%         grid on
%         hold on
%         
%         
%     end
%     legend(unique_cell_labels,'Location', 'northeast')
%     Results.cell_labels=cell_labels;
%     Results.c_true_labels=c_true_labels;
%     Results.colorMARK_labels=colorMARK_labels;
%     Results.unique_cell_labels=unique_cell_labels;
% end
% 
% 
