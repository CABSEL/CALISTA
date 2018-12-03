function consensus=get_consensus(population,p)

[runs,nvars]=size(population);
expected_clusters=length(unique(population(1,:)));
consensus=sparse(nvars,nvars);

% for cycle=1:runs
%     for j = 1:nvars            % only loop over all columns
%         for i = 1:j        % only the top half
%             if population(cycle,j)==population(cycle,i) && i~=j
%                 consensus(i,j) = consensus(i,j)+1;
%                 consensus(j,i) = consensus(j,i)+1;% copy the value from the vector to the matrix
%             end
%         end
%     end
% end
% consensus=consensus+runs*eye(nvars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor (cycle=1:runs,p)
%     cycle
    for clust=1:expected_clusters
        cell_idx_in_cluster=find(population(cycle,:)==clust);
%         consensus(cell_idx_in_cluster,cell_idx_in_cluster)=1;%consensus(cell_idx_in_cluster,cell_idx_in_cluster)+1;
%         %or
        temp_consensus=sparse(nvars,nvars);
        temp_consensus(cell_idx_in_cluster,cell_idx_in_cluster)=1;
        consensus=consensus+temp_consensus;
%                         
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
