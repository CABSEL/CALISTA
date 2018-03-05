#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector up_cell(NumericVector cell,NumericMatrix opt, NumericMatrix log,NumericMatrix mRNA_all,int num_clusters, int n_genes,int ii){
  int nrow=log.nrow();
  NumericVector opt_pram(nrow);
  //Rcout <<cell;
  // for(clust in 1:num_clusters){
  //   for (j in 1:n_genes) {
  //     opt_param_each_gene=log_p_mat_ma[,opt_idx_clusters[clust,j]]
  //     cell_prob[i,clust]=cell_prob[i,clust]+opt_param_each_gene[mRNA_all[i,j]+1]
  //   }
  // }
  for(int clust =0; clust<num_clusters; clust++){
    for(int j=0; j<n_genes; j++){
      //Rcout<<opt(clust,j)-1;
      for(int log_n=0;log_n<nrow;log_n++){
        opt_pram[log_n]=log(log_n,opt(clust,j)-1);
      }
      //Rcout<<opt_pram
      //Rcout<<mRNA_all(ii,j);
      //Rcout<< opt_pram[mRNA_all(ii,j)]; // here  mRNA_all(ii,j), not small in R to 1
      cell[clust] += opt_pram[mRNA_all(ii,j)];  
    }
  }
  //Rcout<<cell;
  return cell;
}
