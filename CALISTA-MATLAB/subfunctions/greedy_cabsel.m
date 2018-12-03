function [my_results ] = greedy_cabsel( as_all,log_p_mat_ma,k_new,mRNA_all,n_genes ,max_iter,nvars,optimize,opt_idx_a,p,sum_prob_tot,population,loops,expected_clusters,algorithm,display)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
if display && loops>1
    fprintf('\nCALISTA_clustering is running...\n');
    fprintf('Progress:');
    fprintf(['\n' repmat('.',1,loops) '\n\n']);
end

parfor (jj=1:loops,p)
    startIteration(jj) = tic;
    % for jj=1:loops
    as=as_all(jj,:);
    % mRNA counts of lookup table
    [max_mRNA_counts, ~]=size(log_p_mat_ma);
    % zero is not included, therefore -1
    max_mRNA_counts=max_mRNA_counts-1;
    % define array
    num_cells_mRNA=zeros(max_mRNA_counts+1,n_genes);
    % define logical variable
    my_break=false;
    for iterations=1:max_iter
        
        if iterations==max_iter
            fprintf(' Max number of iterations reached \n')
        end
        % sort clusters
        [sorted_as,sortedIDX]=sort(as);
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
        logP_each_gene_each_clust=zeros(num_clusters,n_genes);
        %loop over clusters
        for clust=1:num_clusters
            %  clust
            cells_in_each_cluster=data_sorted(bounds(clust,1):bounds(clust,2),:);
            % if cluster has only one entry
            if bounds(clust,1)==bounds(clust,2)
                idx = sub2ind(size(num_cells_mRNA),cells_in_each_cluster+1, 1:n_genes);
                num_cells_mRNA=zeros(max_mRNA_counts+1,n_genes);
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
            [max_L,idx_max_L]=max(Z);
            opt_idx_clusters(clust,:)=idx_max_L;
            logP_each_gene_each_clust(clust,:)=max_L;
        end
%         weights=-logP_each_gene_each_clust';
%         weights=(weights./repmat(sum(weights,2),1,n_genes))';
        % Check if converged or number of maxiterations
        if my_break==true
            break
        end
        % vector for new index
        idx_max_cell_prob=zeros(1,nvars);
        % check if user want to optimize
        
        % vector for the cell probabilities in each cluster
        cell_prob=zeros(nvars,num_clusters);
        my_distance=cell_prob;
%                 % loop over all cells
%                 for i=1:nvars
%         
%                     % loop over cluster
%                     for clust=1:num_clusters
%                         % loop over genes
%                         for j=1:n_genes
%                             % calculate the lop prob of gene in cell belonging
%                             % to cluster
%                             opt_param_each_gene=log_p_mat_ma(:,opt_idx_clusters(clust,j));
%                             % sum all genes together
%                             cell_prob(i,clust)=cell_prob(i,clust)+opt_param_each_gene(mRNA_all(i,j)+1);
%                         end
%                     end
%         
%                     [sum_prob,idx_max_cell_prob(i)]=max(cell_prob(i,:));
%                     sum_prob_tot(jj)=sum_prob_tot(jj)+sum_prob;
%                     maxrel=1;
%         
%                 end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp_mRNA_all=mRNA_all+1;
        opt_param_each_gene=reshape(log_p_mat_ma(:,opt_idx_clusters'),max_mRNA_counts+1,n_genes,num_clusters);
%         opt_param_each_gene_P=reshape(p_mat_ma(:,opt_idx_clusters'),max_mRNA_counts+1,n_genes,num_clusters);
%         mask=(logP_each_gene_each_clust>-1.5)';
        for j=1:n_genes
%             cell_prob=cell_prob-(squeeze(opt_param_each_gene(temp_mRNA_all(:,j),j,:)).*repmat(weights(j,:),nvars,1));
%             masked_cell_prob3D=zeros(nvars,num_clusters);
%             masked_cell_prob3D=squeeze(opt_param_each_gene(temp_mRNA_all(:,j),j,:));
%             masked_cell_prob3D(:,find(mask(j,:)==0))=0;
%              cell_prob=cell_prob+masked_cell_prob3D;
%             temp=cell2mat(splitapply(@(x1){mean(x1)},cell_prob3D,as'));            
%             aaa(j)=max(nanstd(temp,[],2)./nanmean(temp,2));
            cell_prob=cell_prob+squeeze(opt_param_each_gene(temp_mRNA_all(:,j),j,:));
%             
        end
%         cell_prob=cell_prob./repmat(sum(mask),nvars,1);
%         figure
%         bar(aaa)
%         xlabel('gene')
%         ylabel('max SD')
        
        switch algorithm
            % just go for the maximum
            case 'greedy_cabsel'
                [sum_prob,idx_max_cell_prob]=max(cell_prob,[],2);
                idx_max_cell_prob=idx_max_cell_prob';
                sum_prob_tot(jj)=sum(sum_prob);
                maxrel=1;
                % go for the maxiumum difference between clusters
            case 'greedy_maxdiff'
                [sum_prob,idx_max_cell_prob]=max(cell_prob-repmat(sum(cell_prob,2),1,expected_clusters),[],2);
                idx_max_cell_prob=idx_max_cell_prob';
                sum_prob_tot(jj)=sum(sum_prob);
                maxrel=1;
                % maximum difference and flip coin, aneal, etc. code sabec
                          
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % check how many swaps
        comp=as==idx_max_cell_prob;
        rel=sum(comp);
        % calculate the relative number of changes
        
        if rel/nvars>=maxrel || optimize==false
            % final assignmnent
            population(jj,:)=as;
            my_break=true;
            if optimize==false
                for m=1:expected_clusters
                    
                    aaa=find(as==m);
                    my_distance(aaa,:)=abs(repmat(cell_prob(aaa,m),1,expected_clusters)-cell_prob(aaa,:));%./cell_prob(aaa,m);
                end
                
                clusterprobabilities=2.^cell_prob;
                clusterprobabilities=clusterprobabilities./repmat(sum(clusterprobabilities,2),1,expected_clusters);
                my_results(jj).clusterprobabilities=clusterprobabilities;
                my_results(jj).distance=my_distance;
                my_results(jj).cell_prob=cell_prob;
                my_results(jj).best=as;
                
                
                
            end
            
            
        end
        as=idx_max_cell_prob;
        
    end
    
    
    opt_idx_a(jj).run=opt_idx_clusters;
    if display
        if loops>1
            fprintf('\b|\n')
        end
    end
    endIteration(jj) = toc(startIteration(jj));
end

if optimize==true
    if max(population(1,:))~=max(population(2,:))
        my_results.sum_prob_tot=sum_prob_tot;
    end
    
    
    [~,indx_best]=max(sum_prob_tot);
    best_as=population(indx_best,:);
    my_results.population=population;
    my_results.best=best_as;
    my_results.overallsum=sum_prob_tot;
    %     my_results.clusterprobabilities=clusterprobabilities;
    
    opt_id=opt_idx_a(indx_best).run;
    for j=1:size(opt_id,1)
        optpar(j,:)={k_new(opt_id(j,:),:)};
    end
    my_results.parameter=optpar;
    my_results.runtime_eachloop=endIteration;
    my_results.param_idx=opt_id;
else
    my_results.population=[];
    %     my_results.best=as;
    %     my_results.clusterprobabilities=clusterprobabilities;
    opt_id=opt_idx_a.run;
    for j=1:size(opt_id,1)
        optpar(j,:)={k_new(opt_id(j,:),:)};
    end
    my_results.parameter=optpar;
    my_results.param_idx=opt_id;
    %     my_results.optimal_parameter=optpar;
end

end

