get_prob_transition_genes<-function(selected_clusters,log_p_mat_ma,temp_mRNA_all,n_genes,nvar_temp){
  max_mRNA_counts=nrow(log_p_mat_ma)
  max_mRNA_counts=max_mRNA_counts-1
  num_cells_mRNA=matrix(0,max_mRNA_counts+1,n_genes)
  sorted_as=sort(selected_clusters)
  sortedIDX=order(selected_clusters)
  data_sorted=temp_mRNA_all[sortedIDX,]
  Clusters=unique(sorted_as)
  lastIDX=length(sorted_as)-match(Clusters,rev(sorted_as))+1
  #calculate number of clusters
  num_clusers=length(Clusters)
  #bouds of each Clusters
  bounds=as.matrix(data.frame(c(1,lastIDX[1:length(lastIDX)-1]+1),lastIDX))
  #invert the prabalitity matrix
  X=t(log_p_mat_ma)  
  #define array
  opt_idx_clusters=matrix(0,num_clusers,n_genes)
  #loop over each cluster
  for (clust in 1:num_clusers){
    cells_in_each_cluster=data_sorted[bounds[clust,1]:bounds[clust,2],]
    if(bounds[clust,1] == bounds[clust,2]){
      idx=((1:n_genes)-1)*nrow(num_cells_mRNA)+cells_in_each_cluster+1
      num_cells_mRNA[idx]=1
    }else{
      num_cells_mRNA=apply(cells_in_each_cluster,2,function(x){
        hist(x,seq(-0.5,max_mRNA_counts+0.5,1),plot = FALSE)[[2]]
      })
    }
    Y=Matrix(num_cells_mRNA,sparse = TRUE)
    Y=Y /sum(Y[,1])
    Z=X %*% Y
    idx_max_L=apply(Z,2,function(x){
      which(x==max(x))[1]
    })
    opt_idx_clusters[clust,]=idx_max_L
  }
  #idx_max_cell_prob=matrix(0,1,nvar_temp)
  cell_prob_3D=array(0,c(nvar_temp,num_clusers,n_genes))
  for(i in 1:nvar_temp){
    for(clust in 1:num_clusers){
      for(j in 1:n_genes){
        opt_param_each_gene=log_p_mat_ma[,opt_idx_clusters[clust,j]]
        cell_prob_3D[i,clust,j]=opt_param_each_gene[temp_mRNA_all[i,j]+1]
      }
    }
  }
  prob_each_gene=matrix(0,num_clusers,n_genes)
  for(i in 1:n_genes){
    for(j in 1:num_clusers){
      prob_each_gene[j,i]=sum(cell_prob_3D[selected_clusters==j,j,i])
    }
  }
  
  return(prob_each_gene)
}