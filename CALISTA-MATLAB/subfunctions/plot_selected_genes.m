function [Results]= plot_selected_genes(Results)
Xcoordinate=Results.ORDERING.Xcoordinate;
Ycoordinate=Results.ORDERING.Ycoordinate;
Zcoordinate=Results.ORDERING.Zcoordinate;

normed_cell_ordering=Results.ORDERING.normed_cell_ordering;
h=Results.TRANSITION.final_graph;
x_center=zeros(Results.expected_clusters,1);
% y_center=Results.TRANSITION.y_center;
% z_center=Results.TRANSITION.z_center;
for i=1:Results.expected_clusters
    x_center(i)=mean(normed_cell_ordering(Results.final_groups==i));
end




selected_data=Results.PATH.selected_data;
selected_genes=Results.PATH.selected_genes;
cluster_text=[1 find(Results.cluster_progression==max(Results.cluster_progression))];
max_plots_in_each_fig=26;
count2=1;
num_genes_selected=length(selected_genes);
if num_genes_selected>max_plots_in_each_fig
    figs=ceil(num_genes_selected/12);
    [sub,~]=numSubplots(max_plots_in_each_fig);
else
    figs=1;
end

for g=1:num_genes_selected
    figure(6000+ceil(g/max_plots_in_each_fig))
%     set(gcf,'units','points','position',[100,100,1000,1000])
    subplot(7,4,count2)

    colormark=zeros(size(selected_data,1),3);
%     figure
    p = plot(h);
    layout(p,'force','XStart',p.XData,'YStart',1-x_center);
    hold on
    p.EdgeAlpha=1;
    % p.NodeColor=Results.colorMARK_calista;
    p.LineWidth=0.001;
    % p.EdgeColor=[0 0 0];
    p.NodeLabel=[];
    p.MarkerSize=1;
    p.EdgeLabel=[];%abs(h.Edges.Weight);
    p.ZData=1-x_center;
%     for i=1:length(cluster_text)
%         text(p.XData(cluster_text(i))+0.1,p.YData(cluster_text(i))+0.1,p.ZData(cluster_text(i))+0.1,num2str(cluster_text(i)),'FontSize',15);
%     end
    title(selected_genes(g))
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    hold on
    temp_data=selected_data(:,g);
    idx_gray=find(temp_data==0);
    idx_not_gray=find(temp_data~=0);
    exp_temp_data=temp_data(idx_not_gray);
    unique_exp_temp_data=unique(exp_temp_data);
    unique_c_temp=linspace(1,0,length(unique_exp_temp_data))';%flipud(pink(length(unique_exp_temp_data)));
    unique_c=[ones(length(unique_exp_temp_data),1) unique_c_temp unique_c_temp];
    [~,pick_color]=ismember(exp_temp_data,unique_exp_temp_data);
    temp_colormark= unique_c(pick_color,:);
    colormark(idx_gray,:)=repmat([0.8 0.8 0.8],length(idx_gray),1);
    colormark(idx_not_gray,:)=temp_colormark;%repmat([1 0 0],length(idx_not_gray),1);
    markersize(idx_gray,:)=30;
    markersize(idx_not_gray,:)=40;
    % Set colormap
    scatter3(Xcoordinate,Ycoordinate,Zcoordinate, markersize, colormark, 'fill')
    
    if count2==max_plots_in_each_fig
        count2=1;
    else
        count2=count2+1;
    end
end
end

