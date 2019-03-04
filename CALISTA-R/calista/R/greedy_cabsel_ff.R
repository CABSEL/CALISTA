#' A calista Function
#'
#' @param as_all,log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,optimize,opt_idx_a,p,sum_prob_tot,population,loops,expected_clusters,algorithm,display
#' @keywords calista
#' @export
#' @examples
#'
#' greedy_cabsel_ff()




greedy_cabsel_ff<-function(as_all,log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,
                        optimize,opt_idx_a,p,sum_prob_tot,population,loops,expected_clusters,algorithm,display){

  if (display && loops>1){
    writeLines('\nCALISTA_clustering is running...\n')
  }

  # library(doParallel)
  # library(Matrix)
  # library(Rcpp)

  ncore=makeCluster(p)
  registerDoParallel(ncore)
  my_result_tep=list()  #here mark is very important to make my_reslut lengt not be 0!!!!
  my_results=list()
  overallsum = matrix(0,1,50)
  optpar = list()
  #sourceCpp('./R/up_cell.cpp')
  t1 = Sys.time()
  #time3 <- system.time({
    #clusterExport(cl, c("n")) # Export max number of iteration to workers
    items = list()
    my_for<-foreach(jj=1:loops, .multicombine = TRUE, .export = "greedy_cabsel_f", .packages = c('Matrix',"Rcpp","tcltk")) %dopar% {             #.export = 'my_results' .combine = 'c',.multicombine = TRUE

      if (optimize==1){
        if(!exists("pb")) pb <- tkProgressBar("CALISTA clustering progress", min=1, max=loops)
        setTkProgressBar(pb, jj)
        Sys.sleep(0.01)
        log2(jj)
      }

      #for (jj in 1:loops){
      #sink("foreach.txt",append = TRUE)
      #sourceCpp('./R/up_cell.cpp')
      opt_when_false=list()
      #opt_when_false$as_all_as=as_all[[jj,]]

      as=as_all[jj,]
      list(
      my_results_f = greedy_cabsel_f(as,log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,optimize,p,
                                     sum_prob_tot,population,1,expected_clusters,algorithm,display,opt_idx_a)
      )
      # opt_idx_a[[jj]]  = my_results_f$opt_idx_a
      # population[jj,]  = my_results_f$population
      # sum_prob_tot[jj]   = my_results_f$overallsum
      #sum_prob_tot[jj] = my_results_f$sum_prob_tot

      #print(paste("this time as_",jj," is:"))
      #print(as)

    }
    #stopImplicitCluster(ncore)
    #stopImplicitCluster()
    stopCluster(ncore)
    t2 = Sys.time()
    t = t2-t1

    for (kk in 1:loops){

      population [kk,] = my_for[[kk]]$my_results_f$population[1,]
      opt_idx_a [[kk]] = my_for[[kk]]$my_results_f$opt_idx_a
      sum_prob_tot [kk] = my_for[[kk]]$my_results_f$overallsum

    }

      #my_results$aaa=my_for
  #   if(optimize==TRUE){
      indx_best=which(sum_prob_tot==max(sum_prob_tot))[1]
      best_as=population[indx_best,]
      my_results$population=population
      my_results$best=best_as
      my_results$overallsum=sum_prob_tot
      my_results$cluster = expected_clusters
      #my_results$runtime_eachloop = my_results_f$runtime_eachloop

  #     #print(indx_best)
      opt_id=opt_idx_a[[indx_best]]
      for (j in 1:nrow(opt_id)) {
        optpar[[j]]=K_new[opt_id[j,],]
      }
      my_results$parameter=optpar

      # browser()

      return(my_results)
  #   }else{
  #     my_results$population=NULL
  #     opt_id=opt_idx_a[[1]]
  #
  #   browser()
  #   if(loops>1){
  #     for (ii in 1:loops) {
  #       opt_idx_a[[ii]]=opt_idx_a_temp[[ii]]$opt_idx_clusters_temp
  #       population[ii,]=opt_idx_a_temp[[ii]]$population_jj
  #       sum_prob_tot[ii]=opt_idx_a_temp[[ii]]$sum_pro_tot_jj
  #     }
  #   }else
  #   {
  #     opt_idx_a[[1]]=opt_idx_a_temp[[1]]$opt_idx_clusters_temp
  #     population=opt_idx_a_temp[[1]]$population_jj
  #     my_results$clusterprobabilities=opt_idx_a_temp[[1]]$clusterprobabilities
  #     my_results$distance=opt_idx_a_temp[[1]]$distance
  #     my_results$cell_prob=opt_idx_a_temp[[1]]$cell_prob
  #     my_results$best=opt_idx_a_temp[[1]]$best
  #
  #   }
  #   optpar=list()
  #    #here ,watch out
  #     #optpar=matrix(0,nrow = nrow(opt_id),ncol = ncol(K_new))
  #     for(j in 1:nrow(opt_id)){
  #       optpar[[j]]=K_new[opt_id[j,],]
  #     }
  #     my_results$parameter=optpar
  #   }
  # })
}



