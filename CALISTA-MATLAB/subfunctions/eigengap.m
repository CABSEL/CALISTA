function [expected_clusters]=eigengap(consensus,max_clust,varargin)

min_clust=1;
auto=0;
nVarargs=length(varargin);
for k=1:2:nVarargs
      
    % check if user only wants to get the optimal rate constants for a
    % given assignment
    if strcmp(varargin{k},'min_clust')
        % 1 assignement to calculate optimal rate constants
        min_clust=varargin{k+1};
       
    end
    
    if strcmp(varargin{k},'auto')
        % 1 assignement to calculate optimal rate constants
        auto=varargin{k+1};
       
    end
end




% Reorder consensus matrix
consensus=sparse(consensus);
p = symrcm(consensus);
consensus=consensus(p,p);
% calculate degree matrix
tic
degs = sum(consensus, 2);
D    = sparse(1:size(consensus, 1), 1:size(consensus, 2), degs);

% compute unnormalized Laplacian
L = D - consensus;

% compute normalized Laplacian

% avoid dividing by zero
degs(degs == 0) = eps;
% calculate D^(-1/2)
D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));

% calculate normalized Laplacian
L = D * L * D;

% compute the eigenvectors corresponding to the max_clust smallest
% eigenvalues
% max_clust=20;
diff   = eps;
[~,V] = eigs(L, max_clust, diff);
toc
v=diag(V);
v=sort(v,'ascend');
selected_clust=min_clust:max_clust;
all_diff=get_diff(v(selected_clust));
[~,ind_maxdiff]=sort(all_diff,'descend');

if auto==0
    figure()
    title('Eigengap')
    xlabel('Eigenvalues')
    plot(selected_clust,v(selected_clust),'*');
    hold on
    h1=plot(selected_clust(ind_maxdiff(1)),v(selected_clust(ind_maxdiff(1))),'r*');
end
if length(selected_clust)<=5
    top_clust=selected_clust(ind_maxdiff(1));
    if auto==0
        legend(h1,sprintf( '%s %4i', 'First Eigengap: ', selected_clust(ind_maxdiff(1))))
        ylim([0 1])
        hold off
        title('Eigengap values')
        xlabel('Number of clusters')
        promt=sprintf('\nOptimal number of cluster according to max. eigenvalue: %d  \nIf you want to use this value press enter, \nelse provide desired number of cluster: ',top_clust);
    end
    
    
    
else
    top_clust=max(selected_clust(ind_maxdiff(1:3)));
    if auto==0
        h2=plot(selected_clust(ind_maxdiff(2)),v(selected_clust(ind_maxdiff(2))),'m*');
        h3=plot(selected_clust(ind_maxdiff(3)),v(selected_clust(ind_maxdiff(3))),'g*');
        legend([h1 h2 h3],{sprintf( '%s %4i', 'First Eigengap: ', selected_clust(ind_maxdiff(1))),sprintf( '%s %4i', 'Second Eigengap: ', selected_clust(ind_maxdiff(2))),sprintf('%s %4i', 'Third Eigengap: ', selected_clust(ind_maxdiff(3)))})
        
        ylim([0 1])
        hold off
        title('Eigengap values')
        xlabel('Number of clusters')
        promt=sprintf('\nOptimal number of cluster according to max. of 3 eigenvalues: %d  \nIf you want to use this value press enter, \nelse provide desired number of cluster: ',top_clust);
    end
    
    
end

if auto==0
    input1=input(promt);
    
    if isempty(input1)
        expected_clusters=top_clust;
    else
        expected_clusters=input1;
    end
else
    fprintf('\nOptimal number of clusters: ',top_clust);
    expected_clusters=top_clust;
end


% % in case of the Jordan-Weiss algorithm, we need to normalize
% % the eigenvectors row-wise
%
% U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));


end