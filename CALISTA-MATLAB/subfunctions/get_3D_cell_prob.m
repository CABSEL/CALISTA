function [cell_prob_3D]=get_3D_cell_prob(Results,Parameters,DATA)


mRNA_all=DATA.totDATA;
final_groups=Results.final_groups;
log_p_mat_ma=Parameters.logP;
p_mat_ma=Parameters.P;
nvars=DATA.nvars;
n_genes=DATA.numGENES;

% mRNA counts of lookup table
[max_mRNA_counts, ~]=size(log_p_mat_ma);
% zero is not included, therefore -1
max_mRNA_counts=max_mRNA_counts-1;
% define array
num_cells_mRNA=zeros(max_mRNA_counts+1,n_genes);

% sort clusters
[sorted_as,sortedIDX]=sort(final_groups);
% sort mRNA based on clusters
data_sorted=mRNA_all(sortedIDX,:);
% calculate unique clusters
[clusters,lastIDX]=unique(sorted_as,'last');
% calculate number of clusters
num_clusters=length(clusters);
% bounds of each clusters
bounds=[1 lastIDX(1:end-1)'+1; lastIDX']';
% invert the probability matrix
X=log_p_mat_ma';
% define array
opt_idx_clusters=zeros(num_clusters,n_genes);
%loop over clusters
for clust=1:num_clusters
    %  clust
    cells_in_each_cluster=data_sorted(bounds(clust,1):bounds(clust,2),:);
    % if cluster has only one entry
    if bounds(clust,1)==bounds(clust,2)
        idx = sub2ind(size(num_cells_mRNA),cells_in_each_cluster+1, 1:n_genes);
        num_cells_mRNA(idx)=1;
    else
        % number of cells with mRNA count
        num_cells_mRNA(:,:)=hist(cells_in_each_cluster,0:1:max_mRNA_counts);
    end
    % Matrixmultiplication
    Y=sparse(num_cells_mRNA);
    Y=Y/sum(Y(:,1));
    Z=X*Y;
    % get the maximum for each gene
    [~,idx_max_L]=max(Z);
    opt_idx_clusters(clust,:)=idx_max_L;
end

% loop over all cells
for i=1:nvars
    
    % loop over cluster
    for clust=1:num_clusters
        % loop over genes
        for j=1:n_genes
            % calculate the lop prob of gene in cell belonging
            % to cluster
            opt_param_each_gene=p_mat_ma(:,opt_idx_clusters(clust,j));
            % sum all genes together
            cell_prob_3D(i,clust,j)=opt_param_each_gene(mRNA_all(i,j)+1);
        end
    end
end

            
            