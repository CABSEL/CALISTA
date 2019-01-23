function CALISTA_landscape_plotting_main(INPUTS,Results,varargin)

% CALISTA_landscape_plotting_main returns the landscape plots of single
% cells in the dataset based on cell-likelihood values calculated after cell
% ordering.
% 
% Usage:
% 
% CALISTA_landscape_plotting_main(INPUTS,Results)
% 
% Inputs: 
%
% In addition to the specification in the previous steps users need to specify:
%
% ** INPUTS.ngrid **
% Grid size for 2D interpolation
% INPUTS.ngrid=40 (by default)
%
% ** INPUTS.gridshift **
% Shift for cell representation on the landscape
% INPUTS.gridshift=3 (by default)
%
% Results - a structure of CALISTA clustering and lineage inference results
% Run 'CALISTA_clustering_main', 'CALISTA_transition_main' and CALISTA_ordering_main_2
%
%
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright.  October 1, 2018.

% Error check
if nargin <2
    error('Not enough input arguments')
end

if ~isfield(INPUTS,'ngrid')
     INPUTS.ngrid=40; 
end

if ~isfield(INPUTS,'gridshift')
INPUTS.gridshift=3;
end

st1=1;
st2=2;
nVarargs=length(varargin);
for k=1:2:nVarargs
     if strcmp(varargin{k},'comp')
         st=varargin{k+1};
         st1=st(1);
         st2=st(2);
     end
end
cell_likelihood=abs(Results.ORDERING.max_log_P_cell);
% cell_likelihood=exp(cell_likelihood);
%% Plot likelihood colormap
[sorted_values,index_sorted]=sort(cell_likelihood);
figure
colorMARK_likelihood=winter(length(sorted_values));
% colorMARK_likelihood=colorMARK_likelihood(index_sorted,:);

st3=3;
% COMP1=Results.ORDERING.Xcoordinate;
% COMP2=Results.ORDERING.Zcoordinate;
% COMP1=COMP1(index_sorted);
% COMP2=COMP2(index_sorted);
% scatter3(COMP1,COMP2,ones(length(COMP1),1),30,colorMARK_likelihood,'fill')
COMP1=Results.score3(:,st1);
COMP2=Results.score3(:,st2);
COMP3=Results.score3(:,st3);

COMP1=COMP1(index_sorted);
COMP2=COMP2(index_sorted);
COMP3=COMP3(index_sorted);
scatter3(COMP1,COMP2,COMP3,30,colorMARK_likelihood,'fill')


% title('Landscape - Pseudotime plot')
xlabel(sprintf('%s%1i', 'PC',st1))
ylabel(sprintf('%s%1i', 'PC',st2))
zlabel(sprintf('%s%1i', 'PC',st3))

%% LANDSCAPE PLOTTING

% Remove outliers
mean_likelihood=abs(mean(cell_likelihood));
std_likelihood=abs(std(cell_likelihood));
thr1=mean_likelihood+3*std_likelihood;
thr2=mean_likelihood-3*std_likelihood;
not_outlier_index=find(cell_likelihood<thr1 & cell_likelihood>thr2);
COMP1=Results.score3(not_outlier_index,st1);
COMP2=Results.score3(not_outlier_index,st2);
% COMP1=Results.ORDERING.Xcoordinate(not_outlier_index);
% COMP2=Results.ORDERING.Zcoordinate(not_outlier_index);
xq_grid=linspace(min(COMP1),max(COMP1),INPUTS.ngrid);
yq_grid=linspace(min(COMP2),max(COMP2),INPUTS.ngrid);

x=COMP1;
y=COMP2;
if size(cell_likelihood,1)>=size(cell_likelihood,2)
    v=cell_likelihood(not_outlier_index);
else
    v=cell_likelihood(not_outlier_index)';
end
vq=gridfit(x,y,v,xq_grid,yq_grid);
figure
set(0,'DefaultAxesColor','none')
% subplot(122)
s=surf(xq_grid,yq_grid,vq,'FaceAlpha',1);
s.EdgeColor = 'none';
colormap bone
s.EdgeColor = 'none';
% shading interp
% lightangle(-45,30)
lightangle(45,30)
s.FaceLighting = 'gouraud';
s.AmbientStrength = 0.3;
s.DiffuseStrength = 0.8;
s.SpecularStrength = 0.9;
s.SpecularExponent = 25;
s.BackFaceLighting = 'unlit';
hold on
% Now for each point find the closest point on the surface
[~,closetest_z(:,1)]=min((repmat(x,1,length(xq_grid))-repmat(xq_grid,length(x),1)).^2,[],2);
[~,closetest_z(:,2)]=min((repmat(y,1,length(yq_grid))-repmat(yq_grid,length(y),1)).^2,[],2);
% x_new=xq_grid(closetest_z(:,1));
% y_new=yq_grid(closetest_z(:,2));
linearInd = sub2ind(size(vq),closetest_z(:,2),closetest_z(:,1));
v_new=vq(linearInd)+INPUTS.gridshift;


% Define colormap for cell ordering
colorMARK_landscape_cells=Results.ORDERING.orderingCOLORMARK(not_outlier_index,:);
scatter3(x,y,v_new,50,colorMARK_landscape_cells,'fill')
title('Landscape - Pseudotime plot')
xlabel(sprintf('%s%1i', 'PC',st1))
ylabel(sprintf('%s%1i', 'PC',st2))
zlabel('Negative log likelihood')

figure
colorMARK_CALISTA_C=Results.c(not_outlier_index,:);
% subplot(121)
s=surf(xq_grid,yq_grid,vq,'FaceAlpha',1);
colormap bone
s.EdgeColor = 'none';
% shading interp
% lightangle(-45,30)
lightangle(45,30)
s.FaceLighting = 'gouraud';
s.AmbientStrength = 0.3;
s.DiffuseStrength = 0.8;
s.SpecularStrength = 0.9;
s.SpecularExponent = 25;
s.BackFaceLighting = 'unlit';

hold on
scatter3(x,y,v_new,50,colorMARK_CALISTA_C,'fill')
% lightangle(-45,30)
% f.FaceLighting = 'gouraud';
% f.AmbientStrength = 0.3;
% f.DiffuseStrength = 0.8;
% f.SpecularStrength = 0.9;
% f.SpecularExponent = 25;
% f.BackFaceLighting = 'unlit';

grid on
title('Landscape - Cluster plot')
xlabel(sprintf('%s%1i', 'PC',st1))
ylabel(sprintf('%s%1i', 'PC',st2))
zlabel('Negative log likelihood')

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % h=Results.TRANSITION.final_graph;
% % h.Edges.Weight=abs(h.Edges.Weight);
% % selected_edges=[];
% % fprintf('\nAutomatically detect paths along the graph based on shortestpath..\n')
% % ending_clusters=find(Results.cluster_progression==max(Results.cluster_progression));
% % other_ending_clusters=histcounts((h.Edges.EndNodes(:)),Results.expected_clusters);
% % ending_clusters=unique([ending_clusters other_ending_clusters]);
% % ending_clusters(ending_clusters==1)=[];
% % if length(ending_clusters)==1
% %     
% %     CELL_path{1}=Results.cluster_predicted;
% % else
% %     for i =1: length(ending_clusters)
% %         [shortest_path{i},~,edgepath] = shortestpath(h,1,ending_clusters(i));
% % % TR=shortestpathtree(h,1);
% %         selected_edges=[selected_edges edgepath];
% %     end
% %     CELL_path=shortest_path;
% % end
% 
% 
% add_paths=true;
%     count=1;
%     while add_paths
%         
%         fprintf('%s %4i','Path num: ',count')
%         fprintf('\n ******************************************************** \n')
%         fprintf('\nKey the clusters in the path based on the progression (e.g. [1 2 3 4]): ')
%         CELL_path{count}=input('');
%         fprintf('\nPress 1 to add another path, 0 otherwise: ');
%         continue_add=input('');
%         if continue_add
%             figure(102)
%             count=count+1;
%         else
%             add_paths=false;
%         end
%     end
% n_paths=length(CELL_path);
% % figure
% for path=1:n_paths
%     
% %     if length(ending_clusters)==inf
% %         idx_cells_in_path=1:length(Results.final_groups);
% %     else
%         idx_cells_in_path=[];
%         for idx_cluster_in_path=1:length(CELL_path{path})
%             idx_cells_in_path=[idx_cells_in_path find(Results.final_groups==CELL_path{path}(idx_cluster_in_path))];
%         end
% %     end
%     normed_cell_ordering_cells_in_path= Results.ORDERING.normed_cell_ordering(idx_cells_in_path);
%     cell_likelihood_cells_in_path=cell_likelihood(idx_cells_in_path);
%     figure
%     [~,idx_ordering]=sort(normed_cell_ordering_cells_in_path);
%     radius_width=min(200,round(length(cell_likelihood_cells_in_path)/10));
%     smoothExpr=movmean(cell_likelihood_cells_in_path(idx_ordering),radius_width,1,'Endpoints','discard');
%     smoothExpr_std=movstd(cell_likelihood_cells_in_path(idx_ordering),radius_width,1,'Endpoints','discard');
%     y=smoothExpr;
%     dy=smoothExpr_std;
%     x=(1:length(smoothExpr))';
%     line(x,y,'LineWidth',2)
%     hold on
%     fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 .9 .9],'linestyle','none');
%     line(x,y,'LineWidth',2)
%     xlim([0 length(smoothExpr)])
%     % plot(smoothExpr);
%     title('Moving averaged cell-logP after cell ordering')
%     ylabel('Negative log likelihood')
%     xlabel('Cell Ordering')
%     legendInfo = sprintf( '%s %4i %4i %4i %4i %4i %4i %4i %4i %4i', 'Clusters in path: ', CELL_path{path});
%     legend(legendInfo)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% st1=1;
% st2=3;
% COMP1=Results.score3(not_outlier_index,st1);
% COMP2=Results.score3(not_outlier_index,st2);
% xq_grid=linspace(min(COMP1),max(COMP1),20);
% yq_grid=linspace(min(COMP2),max(COMP2),20);
% [xq,yq] = meshgrid(xq_grid, yq_grid);
% x=COMP1;
% y=COMP2;
% v=cell_likelihood(not_outlier_index)';
% vq = griddata(x,y,v,xq,yq,'v4');
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