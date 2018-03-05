function consensus=get_consensus(population)

[runs,nvars]=size(population);

consensus=zeros(nvars);

for cycle=1:runs
    for j = 1:nvars            % only loop over all columns
        for i = 1:j        % only the top half
            if population(cycle,j)==population(cycle,i) && i~=j
                consensus(i,j) = consensus(i,j)+1;
                consensus(j,i) = consensus(j,i)+1;% copy the value from the vector to the matrix
            end
        end
    end
end

consensus=consensus+runs*eye(nvars);