function [Results]=hierarchical(Results,DATA)

fprintf('\nHierarchical clustering based on cell ordering...\n')
color_hierarchical=Results.c;

ColumnLabelsColorValue.Labels=cellstr(num2str(Results.ORDERING.idx_sorted_cells));
ColumnLabelsColorValue.Colors= num2cell(color_hierarchical(Results.ORDERING.idx_sorted_cells,:),2);

transition_expression=DATA.totDATA(Results.ORDERING.idx_sorted_cells,Results.GENES.idx_tot_transition_genes);
transition_expression=zscore(transition_expression); % stardardization
transition_expression=zscore(transition_expression');
cgo = clustergram(transition_expression,'Cluster',1,'Standardize','none','RowPDist','spearman','Linkage','average','Dendrogram',1);
%     set(cgo,)
set(cgo,'RowLabels',Results.GENES.tot_transition_genes)
set(cgo,'ColumnLabels',Results.ORDERING.idx_sorted_cells,'ColumnLabelsColor',ColumnLabelsColorValue,'LabelsWithMarkers',true);
%     addXLabel(cgo, ' Cell Ordering ')
%     addYLabel(cgo, ' Hierarchical clustering ')
addTitle(cgo,' Hierarchical clustering ')
Results.ORDERING.cgo=cgo;
NumCluster=input('Key the number of gene clusters: ');
Z = linkage(transition_expression,'average','spearman');
figure
dendrogram(Z,0,'Orientation','left');
c = cluster(Z,'maxclust',NumCluster);
color = Z(end-NumCluster+2,3)-eps;
[~,~,perm] = dendrogram(Z, 0, 'Orientation','left','colorthreshold', color);
permGENES=Results.GENES.tot_transition_genes(perm);
c=c(perm);
for i= 1:NumCluster
    Results.GENES.hierarchicalK{i}= permGENES(find(c==i))';
end