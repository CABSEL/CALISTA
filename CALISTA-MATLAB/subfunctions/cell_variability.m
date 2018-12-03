function [Results]=cell_variability(DATA,Results)

fprintf('\nCalculating cell to cell variability...\n')
% Entropy
[~,idx_sorted_groups]=sort(Results.final_groups);
mRNA_all=DATA.totDATA(idx_sorted_groups,:);
cutDIMENSION= nonzeros(histcounts(Results.final_groups));

DDD= mat2cell(mRNA_all',DATA.numGENES,cutDIMENSION);

min_num_cells=10;
[ H] = entropy_calculation( DDD,min_num_cells );

figure(111)
subplot(121)
boxplot(H')
xlabel('Cluster')
ylabel('Entropy')
title('Boxplot for the entropy')
set(gca,'xtick',Results.cluster_predicted)
figure(111)
subplot(122)
plot(mean(H'),'c')
hold on
plot(median(H'),'m')
title('Entropy Mean and Median')
legend('MeanEntropy','MedianEntropy')
xlabel('Cluster')
ylabel('Entropy value')
set(gca,'xtick',Results.cluster_predicted)
% subplot(223)
% plot(sum(H,2),'k')
% title('Total Entropy')
% xlabel('Cluster')
% ylabel('Entropy value')
% [dd(:,1), dd(:,2)]=max(H);
% set(gca,'xtick',Results.cluster_predicted)

% Cell-Cell correlation

for i=1:Results.expected_clusters
    aaa=find(Results.final_groups==Results.cluster_predicted(i));
    cell_cell_correlation{i}=corr(DATA.totDATA(aaa,:)');
    temp_corr=sort(cell_cell_correlation{i}(:),'descend');
    temp_corr=temp_corr(length(aaa)+1:end);
    mean_cell_cell_correlation(i)=mean(abs(temp_corr));   
end

% subplot(223)
% plot(mean_cell_cell_correlation)
% title('Mean cell-cell correlation in each cluster')
% xlabel('Cluster')
% ylabel('Correlation')
% set(gca,'xtick',Results.cluster_predicted)
% 
% Results.cluster_entropy=H;
% Results.mean_cell_cell_correlation(i)=mean_cell_cell_correlation(i);
% 
% subplot(224)
% plot(1:Results.expected_clusters,Results.mean_prob_in_out_cluster(:,1))
% set(gca,'XTick',(1:1:Results.expected_clusters))
% title('Mean Probability')
% xlabel('Cluster')
% ylabel('Log P')

pause(1)