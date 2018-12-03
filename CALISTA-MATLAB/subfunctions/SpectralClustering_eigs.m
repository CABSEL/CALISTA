function [U] = SpectralClustering_eigs(consensus,k)
%SPECTRALCLUSTERING Executes spectral clustering algorithm
%   Executes the spectral clustering algorithm defined by
%   Type on the consensus matrix and returns the cluster
%   assignments.
%   If L and U are also called, the (normalized) Laplacian and
%   eigenvectors will also be returned.
%
%   Spectral clustering algorithm used:
%   Normalized according to Jordan and Weiss (2002)
%
%   References:
%   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering",
%     Statistics and Computing 17 (4), 2007
%

% calculate degree matrix
degs = sum(consensus, 2);
D    = sparse(1:size(consensus, 1), 1:size(consensus, 2), degs);

% compute unnormalized Laplacian
L = D - consensus;


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
[U,~] = eigs(L, k, diff);

% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise

U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));


end