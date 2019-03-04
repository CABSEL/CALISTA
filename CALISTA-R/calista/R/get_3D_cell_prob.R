#' A calista Function
#'
#' @param Results,Parameters,DATA
#' @keywords calista
#' @export
#' @examples
#' get_3D_cell_prob()






get_3D_cell_prob<-function(Results,Parameters,DATA){

  mRNA_all=DATA$totDATA
  final_groups=Results$final_groups
  log_p_mat_ma=Parameters$Parameters[[3]]
  p_mat_ma=Parameters$Parameters[[2]]
  nvars=DATA$nvars
  n_genes=DATA$numGENES
  #mRNA counts of lookup table
  max_mRNA_counts=nrow(log_p_mat_ma)
  #zero is not included, therefore -1
  max_mRNA_counts=max_mRNA_counts-1
  #define array
  num_cells_mRNA=matrix(0,max_mRNA_counts+1,n_genes)
  #sort cluster
  sorted_as=sort(final_groups)
  sortedIDX=order(final_groups)
  #sort mRNA based on clusters
  data_sorted=mRNA_all[sortedIDX,]
  #calculate unique clusters
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
  cell_prob_3D=array(0,c(nvars,num_clusers,n_genes))
  #loop over all cells
  for(i in 1:nvars){
    for(clust in 1:num_clusers){
      for(j in 1:n_genes){
        opt_param_each_gene=p_mat_ma[,opt_idx_clusters[clust,j]]
        cell_prob_3D[i,clust,j]=opt_param_each_gene[mRNA_all[i,j]+1]
      }
    }
  }
  return(cell_prob_3D)

}
