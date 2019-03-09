function [Results]=visualization(INPUTS,DATA,Results)

fprintf('\nPlotting...\n')
expected_clusters=Results.expected_clusters;
final_groups=Results.final_groups;
% totDATA=log2(mRNA_all+1);
totDATA=DATA.totDATA.^(1/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeros_cells=sum(totDATA==0,2)*100/size(totDATA,2);
% [aa,bb]=unique(round(zeros_cells));
dot_size=30*ones(size(totDATA,1),1);%10+round(zeros_cells/2);%30*ones(size(totDATA,2),1);%round(zeros_cells/2);%30*ones(size(totDATA,2),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

    
% ********** PCA **********
[coeff, score3,latent] = pca(zscore(totDATA));
Results.coeff=coeff;
explainedVar=latent*100/sum(latent);

figure;
bar(explainedVar)
title('Explained Variance in %')
grid on
xlabel('PC')
ylabel('Variance in %')

if INPUTS.plot_tsne==1
    % ********** t-SNE **********
    %set parameters
    %         initial_dims=DATA.numGENES;
    disp('******************************************************************************')
    disp('*** ATTENTION: for tsne plot please visit https://github.com/KlugerLab/FIt-SNE')
    disp('and install FFTW from http://www.fftw.org ***')
    disp('******************************************************************************')
    opts.perplexity=INPUTS.tsne_opts.perplexity;
    if INPUTS.tsne_opts.input_data == 0
        score2=fast_tsne(totDATA,opts);
    else
        score2=fast_tsne(score3(:,1:min(size(score3,2),INPUTS.tsne_opts.input_data)),opts);
    end
    Results.score2=score2;
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % *** Colormap for K-means ***
% seed=round(DATA.nvars/3*expected_clusters/3); % for reproducibility
% rng(seed); % For reproducibility 
% dataTOcluster_stded=zscore(totDATA);
% ClusterGroup_kmeans_struct.final_groups= kmeans(dataTOcluster_stded,expected_clusters);
% ClusterGroup_kmeans_struct.expected_clusters=expected_clusters;
% % colorMARK1=winter(eva.OptimalK);
% colorMARK_kmeans=jet(expected_clusters);
% figure
% for i=1:expected_clusters
%     scatter3(score3(ClusterGroup_kmeans_struct.final_groups==i,1), score3(ClusterGroup_kmeans_struct.final_groups==i,2), score3(ClusterGroup_kmeans_struct.final_groups==i,3),30, colorMARK_kmeans(i,:), 'fill');
%     xlabel('PC1')
%     ylabel('PC2')
%     zlabel('PC3')
%     title('Kmeans')
%     grid on
%     hold on
%     legendInfo_kmeans{i} = sprintf( '%s %4i', 'Cluster ', i);
% end
% legend(legendInfo_kmeans,'Location', 'northeast')
% Results.kmeans=ClusterGroup_kmeans_struct.final_groups;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot

hfig=figure(1000);
if ishandle(hfig)
clf(1000)
end
% set(hfig,'position', [500, 500, 1200, 400]) 
subplot(121)
for i=1:DATA.num_time_points
    scatter3(score3(DATA.timeline==DATA.time(i),1), score3(DATA.timeline==DATA.time(i),2),score3(DATA.timeline==DATA.time(i),3),dot_size(DATA.timeline==DATA.time(i)),colorMARK_time(i,:),'fill')
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
    temp=dot_size(ClusterGroup_calista==cluster_predicted(i));
    scatter3(score3(ClusterGroup_calista==cluster_predicted(i),1), score3(ClusterGroup_calista==cluster_predicted(i),2), score3(ClusterGroup_calista==cluster_predicted(i),3),temp, repmat(colorMARK_calista(i,:),size(temp,1),1), 'fill');
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

if INPUTS.plot_tsne==1
    
    hfig=figure(1001);
    if ishandle(hfig)
        clf(1001)
    end
    % set(hfig,'position', [500, 500, 1200, 400])
    subplot(121)
    for i=1:DATA.num_time_points
        scatter(score2(DATA.timeline==DATA.time(i),1), score2(DATA.timeline==DATA.time(i),2),dot_size(DATA.timeline==DATA.time(i)),colorMARK_time(i,:),'fill')
        xlabel('TSNE1')
        ylabel('TSNE2')
        %zlabel('PC3')
        title('Original time/cell stage info')
        grid on
        hold on
        
        legendInfo_time{i} = sprintf( '%s %4i', 'Time/Stage', DATA.time(i));
    end
    legend(legendInfo_time,'Location', 'northeast')
    
    Results.h2=subplot(122);
    for i=1:length(cluster_predicted)
        temp=dot_size(ClusterGroup_calista==cluster_predicted(i));
        scatter(score2(ClusterGroup_calista==cluster_predicted(i),1), score2(ClusterGroup_calista==cluster_predicted(i),2), temp, repmat(colorMARK_calista(i,:),size(temp,1),1), 'fill');
        title('Cell Clustering')
        xlabel('TSNE1')
        ylabel('TSNE2')
        grid on
        hold on
        %     legendInfo_calista{i} = sprintf( '%s %4i', 'Cluster ', i);
    end
    
    legend(Results.legendInfo_calista,'Location', 'northeast')
    pause(3)
    
end





Results.colorMARK_calista=colorMARK_calista;
Results.colorMARK_time=colorMARK_time;
Results.score3=score3;
Results.c=c;

%% Plot time snapshots
% if ~(length(unique(DATA.timeline))==1 && unique(DATA.timeline)==0) 
% figure
% for i=1:length(cluster_predicted)
%     temp=dot_size(ClusterGroup_calista==cluster_predicted(i));
%     scatter3(DATA.timeline(ClusterGroup_calista==cluster_predicted(i),1), score3(ClusterGroup_calista==cluster_predicted(i),1), score3(ClusterGroup_calista==cluster_predicted(i),3),temp, repmat(colorMARK_calista(i,:),size(temp,1),1), 'fill');
%     title('Cell Clustering')
%     xlabel('timeline')
%     ylabel('PC1')
%     zlabel('PC3')
%     grid on
%     hold on
% %     legendInfo_calista{i} = sprintf( '%s %4i', 'Cluster ', i);
% end
% end
% %% Plot with true labels
% true_labels=input('\nKey 1 for true labels visualization, 0 otherwise :');
% 
% if true_labels
%     [FileName,PathName,FilterIndex] = uigetfile('*.*');
%     filename=strcat(PathName, FileName);
%     cell_labels=importdata(filename);
%     if ~isempty(DATA.cut_sort.cell_2_cut)
%         cell_labels(DATA.cut_sort.cell_2_cut)=[]; % remove labels for removed clusters
%     end
%     
%     if ~isempty(DATA.cut_sort.idx2cutCELL)
%         cell_labels(DATA.cut_sort.idx2cutCELL)=[]; % remove labels for removed cells based on %zeros
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
%     figure%('position', [500, 500, 400, 400]) 
%     for i=1:length(unique_cell_labels)
%         scatter3(score3(IC==i,1), score3(IC==i,2),score3(IC==i,3),dot_size(IC==i),colorMARK_labels(i,:),'fill')
%         xlabel('PC1')
%         ylabel('PC2')
%         zlabel('PC3')
%         title('True labels')
%         grid on
%         hold on
%         
%         
%     end
%     figure
%     for i=1:length(unique_cell_labels)
%         scatter(score2(IC==i,1), score2(IC==i,2),dot_size(IC==i),colorMARK_labels(i,:),'fill')
%         xlabel('TSNE1')
%         ylabel('TSNE2')
%         
%         title('True labels')
%         grid on
%         hold on
%         
%         
%     end
%     
%     if iscell(unique_cell_labels)
%         legend(cellstr(unique_cell_labels),'Location', 'northeast')
%     else
%         legend(num2str(unique_cell_labels),'Location', 'northeast')
%     end
%     Results.cell_labels=cell_labels;
%     Results.c_true_labels=c_true_labels;
%     Results.colorMARK_labels=colorMARK_labels;
%     Results.unique_cell_labels=unique_cell_labels;
% end
% 
% 
