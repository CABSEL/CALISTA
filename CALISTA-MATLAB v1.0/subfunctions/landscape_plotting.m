function Results=landscape_plotting(DATA,Results)

cell_likelihood=abs(Results.ORDERING.max_log_P_cell);

%% Plot likelihood colormap
[sorted_values,index_sorted]=sort(cell_likelihood);
figure
colorMARK_likelihood=winter(length(sorted_values));
% colorMARK_likelihood=colorMARK_likelihood(index_sorted,:);
st1=1;
st2=2;
COMP1=Results.score3(:,st1);
COMP2=Results.score3(:,st2);
% Based on pseudotime colors
for k=1:size(colorMARK_likelihood,1)
    scatter3(COMP1(cell_likelihood==sorted_values(k)),COMP2(cell_likelihood==sorted_values(k)),1,30,colorMARK_likelihood(k,:),'fill')
    hold on
end
% title('Landscape - Pseudotime plot')
xlabel(sprintf('%s%1i', 'PC',st1))
ylabel(sprintf('%s%1i', 'PC',st2))
%% LANDSCAPE PLOTTING

% Remove outliers
mean_likelihood=abs(mean(cell_likelihood));
std_likelihood=abs(std(cell_likelihood));
thr1=mean_likelihood+2*std_likelihood;
thr2=mean_likelihood-2*std_likelihood;
not_outlier_index=find(cell_likelihood<thr1 & cell_likelihood>thr2);
st1=1;
st2=2;
COMP1=Results.score3(not_outlier_index,st1);
COMP2=Results.score3(not_outlier_index,st2);
xq_grid=linspace(min(COMP1),max(COMP1),20);
yq_grid=linspace(min(COMP2),max(COMP2),20);
[xq,yq] = meshgrid(xq_grid, yq_grid);
x=COMP1;
y=COMP2;
v=cell_likelihood(not_outlier_index)';
vq = griddata(x,y,v,xq,yq,'v4');
figure
% subplot(122)
s=surf(xq,yq,vq,'FaceAlpha',1);
s.EdgeColor = 'none';
spz = interp2(xq,yq,vq, x,y);
% Build mask of points which are above surface
mask = v < spz;

v(mask)=spz(mask)+0.02*(spz(mask));  % translate points up to the surface plot
colormap bone
hold on

% Define clormap for cell ordering
normed_cell_ordering_no_outliers=Results.ORDERING.normed_cell_ordering(not_outlier_index);
unique_pseudotimes=unique(normed_cell_ordering_no_outliers);
colorMARK_landscape_cells=parula(length(unique_pseudotimes));

% Based on pseudotime colors
for k=1:length(unique_pseudotimes)
    scatter3(COMP1(normed_cell_ordering_no_outliers==unique_pseudotimes(k)),COMP2(normed_cell_ordering_no_outliers==unique_pseudotimes(k)),v(normed_cell_ordering_no_outliers==unique_pseudotimes(k))',30,colorMARK_landscape_cells(k,:),'fill')
    hold on
end
title('Landscape - Pseudotime plot')
xlabel(sprintf('%s%1i', 'PC',st1))
ylabel(sprintf('%s%1i', 'PC',st2))
zlabel('Negative log likelihood')
figure
% subplot(121)
s=surf(xq,yq,vq,'FaceAlpha',1);
s.EdgeColor = 'none';
colormap bone
hold on
% Based on cluster colors
for k=1:Results.expected_clusters
      scatter3(COMP1(Results.final_groups(not_outlier_index)==Results.cluster_predicted(k)),COMP2(Results.final_groups(not_outlier_index)==Results.cluster_predicted(k)),v(Results.final_groups(not_outlier_index)==Results.cluster_predicted(k))',30,Results.colorMARK_calista(k,:),'fill')
      hold on
end
grid on
title('Landscape - Cluster plot')
xlabel(sprintf('%s%1i', 'PC',st1))
ylabel(sprintf('%s%1i', 'PC',st2))
zlabel('Negative log likelihood')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
st1=1;
st2=3;
COMP1=Results.score3(not_outlier_index,st1);
COMP2=Results.score3(not_outlier_index,st2);
xq_grid=linspace(min(COMP1),max(COMP1),20);
yq_grid=linspace(min(COMP2),max(COMP2),20);
[xq,yq] = meshgrid(xq_grid, yq_grid);
x=COMP1;
y=COMP2;
v=cell_likelihood(not_outlier_index)';
vq = griddata(x,y,v,xq,yq,'v4');
% figure
% subplot(122)
% s=surf(xq,yq,vq,'FaceAlpha',1);
% s.EdgeColor = 'none';
% spz = interp2(xq,yq,vq, x,y);
% % Build mask of points which are above surface
% mask = v < spz;
% 
% v(mask)=spz(mask)+0.02*(spz(mask));  % translate points up to the surface plot
% colormap winter
% hold on
% 
% % Define clormap for cell ordering
% normed_cell_ordering_no_outliers=Results.ORDERING.normed_cell_ordering(not_outlier_index);
% unique_pseudotimes=unique(normed_cell_ordering_no_outliers);
% colorMARK_landscape_cells=parula(length(unique_pseudotimes));
% 
% % Based on pseudotime colors
% for k=1:length(unique_pseudotimes)
%     scatter3(COMP1(normed_cell_ordering_no_outliers==unique_pseudotimes(k)),COMP2(normed_cell_ordering_no_outliers==unique_pseudotimes(k)),v(normed_cell_ordering_no_outliers==unique_pseudotimes(k))',30,colorMARK_landscape_cells(k,:),'fill')
%     hold on
% end
% title('Landscape - Pseudotime plot')
% xlabel(sprintf('%s%1i', 'PC',st1))
% ylabel(sprintf('%s%1i', 'PC',st2))
% zlabel('abs Log likelihood')
% subplot(121)
% s=surf(xq,yq,vq,'FaceAlpha',1);
% s.EdgeColor = 'none';
% colormap winter
% hold on
% % Based on cluster colors
% for k=1:Results.expected_clusters
%       scatter3(COMP1(Results.final_groups(not_outlier_index)==Results.cluster_predicted(k)),COMP2(Results.final_groups(not_outlier_index)==Results.cluster_predicted(k)),v(Results.final_groups(not_outlier_index)==Results.cluster_predicted(k))',30,Results.colorMARK_calista(k,:),'fill')
%       hold on
% end
% grid on
% title('Landscape - Cluster plot')
% xlabel(sprintf('%s%1i', 'PC',st1))
% ylabel(sprintf('%s%1i', 'PC',st2))
% zlabel('abs Log likelihood')
end