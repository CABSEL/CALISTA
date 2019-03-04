#' A calista Function
#'
#' @param as_all,log_p_mat_ma,k_new,mRNA_all,n_genes,max_iter,nvars,optimize,p,sum_prob_tot,population,loops,expected_clusters,algorithm,display,opt_idx_a
#' @keywords calista
#' @export
#' @examples
#'
#' greedy_cabsel_f()




greedy_cabsel_f = function(as_all,log_p_mat_ma,k_new,mRNA_all,n_genes,max_iter,nvars,optimize,p,
                           sum_prob_tot,population,loops,expected_clusters,algorithm,display,opt_idx_a) {
  #system("R CMD SHLIB ./src/greedy.f95")
  #print('aaaaaaaa')
  #R CMD SHLIB *.f95
    if (display && loops>1){
    writeLines('\nCALISTA_clustering is running...\n')}

 # dyn.load("greedy.so")

n_genes = as.integer(n_genes)
loops = as.integer(loops)
max_iter = as.integer(max_iter)
nvars = as.integer(nvars)
p = as.integer(p)
expected_clusters = as.integer(expected_clusters)
display = as.integer(display)
optimize = as.integer(optimize)

as_all = as.integer(as.matrix(as_all))
log_p_mat_ma = as.matrix(log_p_mat_ma)
k_new = as.matrix(k_new)
mRNA_all = as.integer(as.matrix(mRNA_all))

sum_prob_tot = as.matrix(sum_prob_tot)
population = as.matrix(population)

  my_results = .Fortran("greedy", as_all, log_p_mat_ma, k_new, mRNA_all, n_genes, max_iter, nvars, optimize, p,
                        sum_prob_tot, population, loops, expected_clusters, display,
                        population = as.integer(matrix(0,loops,nvars)), best = as.integer(matrix(0,1,nvars)),
                        overallsum = matrix(0,1,loops), runtime_eachloop = matrix(0,1,loops),
                        param_idx = as.integer(matrix(0, expected_clusters,n_genes)),
                        parameter = array(0, dim = c(expected_clusters,n_genes,3)),
                        clusterprobabilities = as.double(matrix(1,nvars,expected_clusters)),
                        distance = as.double(matrix(1,nvars,expected_clusters)),
                        cell_prob = matrix(0,nvars,expected_clusters),
                        opt_idx_a = as.integer(matrix(0, expected_clusters,n_genes,loops)))

  my_results$opt_idx_a = matrix(my_results$param_idx, expected_clusters, n_genes,loops)

  my_results$population = matrix(as.double(my_results$population), loops, nvars)
  my_results$best = array(as.double(my_results$best), dim = c(nvars))
  my_results$overallsum = matrix(my_results$overallsum, 1, loops)
  my_results$runtime_eachloop = matrix(my_results$runtime_eachloop, 1, loops)
  my_results$param_idx = matrix(my_results$param_idx, expected_clusters, n_genes)
  my_results$parameter = array(my_results$parameter, dim = c(expected_clusters,n_genes,3))
  my_results$parameter = lapply(seq(dim(my_results$parameter)[1]), function(z) my_results$parameter[z,,])
  my_results$cluster = expected_clusters

  my_results$clusterprobabilities = matrix(as.double(my_results$clusterprobabilities),nvars,expected_clusters)
  my_results$distance = matrix(as.double(my_results$distance),nvars,expected_clusters)
  my_results$cell_prob = matrix(as.double(my_results$cell_prob),nvars,expected_clusters)

  return(my_results)

}
