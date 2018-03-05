function [expected_clusters]=eigengap(consensus,max_clust)

    D=eye(length(consensus)).*sum(consensus);
    L=D-consensus;
    L=D^(-1/2)*L*D^(-1/2);
    [U, V]=eig(L);
    v=diag(V);
    v=sort(v,'ascend');
    figure()
    title('Eigengap')
    xlabel('Eigenvalues')
    all_diff=get_diff(v(1:max_clust));
    [v_diff,ind_maxdiff]=max(all_diff);
    plot(v(1:max_clust),'*');
    hold on
    plot(ind_maxdiff,v(ind_maxdiff),'r*');
    ylim([0 1])
    hold off
    title('Eigengap values')
    xlabel('Number of clusters')
    promt=sprintf('\nOptimal number of cluster according to max. eigenvalue %d:  \nIf you want to use this value press enter, \nelse provide desired number of cluster: ',ind_maxdiff);
    input1=input(promt);
    
    if isempty(input1)
        expected_clusters=ind_maxdiff;
    else
        expected_clusters=input1;
    end
end