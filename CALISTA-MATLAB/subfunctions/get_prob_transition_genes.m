function [prob_each_gene,opt_idx_clusters]=get_prob_transition_genes( selected_clusters,Z_temp,temp_mRNA_all,n_genes,nvars_temp)  

% mRNA counts of lookup table   
    [max_mRNA_counts, ~]=size(Z_temp); 
    % zero is not included, therefore -1
    max_mRNA_counts=max_mRNA_counts-1;
    % define array
    num_cells_mRNA=zeros(max_mRNA_counts+1,n_genes);
   
        % sort clusters
        [sorted_as,sortedIDX]=sort(selected_clusters);
        % sort mRNA based on clusters
        data_sorted=temp_mRNA_all(sortedIDX,:);
        % calculate unique clusters
        [clusters,lastIDX]=unique(sorted_as,'last');
        % calculate number of clusters
        num_clusters=length(clusters);
        % bounds of each clusters
        bounds=[1 lastIDX(1:end-1)'+1; lastIDX']';
        % invert the probability matrix
        X=Z_temp';
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
        
         % vector for new index
        idx_max_cell_prob=zeros(1,nvars_temp);
       
        
        % vector for the cell probabilities in each cluster
        cell_prob=zeros(nvars_temp,num_clusters);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        opt_param_each_gene=reshape(Z_temp(:,opt_idx_clusters'),max_mRNA_counts+1,n_genes,num_clusters);
            
            
                       
            for i=1:n_genes
                for j=1:num_clusters
%                 prob_each_gene(j,i)=sum(squeeze(cell_prob(find(selected_clusters==j),j,i)));
               prob_each_gene(j,i)= sum(squeeze(opt_param_each_gene(temp_mRNA_all(find(selected_clusters==j),i)+1,i,j)));
                end
            end