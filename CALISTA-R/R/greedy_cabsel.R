 greedy_cabsel<-function(as_all,log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,
                            optimize,opt_idx_a,p,sum_prob_tot,population,loops,expected_clusters,algorithm,display){
  
   ####3=test
   #greedy_cabsel(as_all,log_p_mat_ma,K_new,mRNA_all,n_genes,max_iter,nvars,
   #optimize,opt_idx_a,p,sum_prob_tot,in_population,loops,expected_clusters,algorithm,display)
   
  #population=in_population 
  ### some explanation
  if (display && loops>1){
    writeLines('\nCALISTA_clustering is running...\n')
    #writeLines('Progress:')
    #writeLines(paste('\n','..................','\n\n',sep = ''))
  }
  library(doParallel)
  library(Matrix)
  library(Rcpp)
  #ncore=makeCluster(p,outfile="/home/Tao/Pro/AAA-Final_code_reduced_functions/foreach.txt")
  ncore=makeCluster(p)
  registerDoParallel(ncore)
  my_result_tep=list()  #here mark is very important to make my_reslut lengt not be 0!!!!
  my_results=list()
  sourceCpp('./R/up_cell.cpp')
  opt_idx_a_temp<-foreach(jj=1:loops, .multicombine = TRUE,.packages = c('Matrix',"Rcpp")) %dopar% {             #.export = 'my_results' .combine = 'c',.multicombine = TRUE
  #for (jj in 1:loops){
    #sink("foreach.txt",append = TRUE)
    sourceCpp('./R/up_cell.cpp')
    opt_when_false=list()
    #opt_when_false$as_all_as=as_all[[jj,]]
    as=as_all[jj,]
    #print(paste("this time as_",jj," is:"))
    #print(as)
    max_mRNA_counts=nrow(log_p_mat_ma)
    max_mRNA_counts=max_mRNA_counts-1
    num_cells_mRNA=matrix(0,max_mRNA_counts+1,n_genes)
    my_break=FALSE
    
    for (iterations in 1:max_iter){
      if (iterations==max_iter){
        writeLines(' Max number of iterations reached \n')
      }
      sorted_as=sort(as)
      sortedIDX=order(as)
      data_sorted=NULL
      data_sorted=as.matrix(mRNA_all[sortedIDX,])
      data_sorted=unname(data_sorted)
      clusters=unique(sorted_as)
      lastIDX=which((!duplicated(sorted_as,fromLast = TRUE)==TRUE))   #check here
      num_clusters=length(clusters)   #number of clusters would not change here 
      #print('one time here')
      #print(iterations)
      #print(num_clusters)
      bounds=NULL
      if(length(lastIDX)>1){
        bounds=matrix(data = c(1,lastIDX[1:(length(lastIDX)-1)]+1,lastIDX),byrow = TRUE,nrow = 2)
      }else{
        bounds=matrix(data = c(1,lastIDX),byrow = TRUE,nrow = 2)
      }
       #here
      bounds=t(bounds)
      X=t(log_p_mat_ma)
      opt_idx_clusters=matrix(0,nrow = num_clusters,ncol = n_genes)
      for (clust in 1:num_clusters){   #question here
        cells_in_each_cluster=data_sorted[bounds[clust,1]:bounds[clust,2],]
        if (bounds[clust,1]==bounds[clust,2]){
          idx_1=(1:n_genes)*(ncol(num_cells_mRNA))+cells_in_each_cluster+1   ##check here
          idx_1=unname(unlist(idx_1))
          #print(idx)  #less than 597*100
          num_cells_mRNA[idx_1]=1  #here change idx to idx_1
        }else{
          ####problem remains here>>>must center???
          num_cells_mRNA=as.matrix(apply(cells_in_each_cluster,2,function(x) (hist(x,breaks=seq(-0.5,max_mRNA_counts+0.5,1),plot=FALSE)$counts)))      #here
        }
        ###sparse here 
        num_cells_mRNA=as.matrix(num_cells_mRNA)
        Y=Matrix(num_cells_mRNA,sparse = TRUE)
        #Y=num_cells_mRNA
        Y=Y/sum(Y[,1])
        Y=unname(Y)
        X=unname(X)
        Z=X%*%Y  # X * Y : non-conformable arrays ???!!!1 why ? NA values
        idx_max_L=lapply(1:ncol(Z), function(i){
          temp_Z= Z[,i]
          which(temp_Z==max(temp_Z))[1]})#check heres
        #print(idx_max_L)   #597*100times
        opt_idx_clusters[clust,]=unlist(idx_max_L)
      }
      opt_when_false$opt_idx_clusters_temp=opt_idx_clusters
      if(my_break==TRUE){break}
      idx_max_cell_prob=matrix(0,nrow = 1,ncol = nvars)
      cell_prob=matrix(0,nrow = nvars,ncol = num_clusters)
      my_distance=cell_prob
      #----------------
      # for(i in 1:nvars){
      #   # cell_prob[i,]=sapply(1:ncol(cell_prob),function(clust){
      #   #   sum(sapply(1:n_genes,function(j){log_p_mat_ma[,opt_idx_clusters[clust,j]][mRNA_all[i,j]+1]}))
      #   # })
      # cell_i=cell_prob[i,]
      # ii=i-1
      # cell_prob[i,]=up_cell(cell_i,opt_idx_clusters,log_p_mat_ma,mRNA_all,num_clusters,n_genes,ii)
      #   
      # switch (algorithm,
      #     'greedy_cabsel' = {
      #       sum_prob=max(cell_prob[i,])
      #       idx_max_cell_prob[i]=which(cell_prob[i,]==sum_prob)[1]  #alway notice that, not just one max value 
      #       sum_prob_tot[jj]=sum_prob_tot[jj]+sum_prob #this is different as matlab
      #       
      #       maxrel=1
      #     },
      #     'greedy_maxdiff'={
      #       tmp=cell_prob[i,]-sum(cell_prob[i,])
      #       sum_prob=max(tmp)
      #       idx_max_cell_prob[i]=which(tmp==sum_prob)
      #       sum_prob_tot[jj]=sum_prob_tot[jj]+sum_prob
      #       maxrel=1
      #     },
      #     'cabsel_sabec'={
      #       maxrel=0.95
      #       p_temp=cell_prob[i,]-sum(cell_prob[i,])
      #       p_temp=p_temp/sum(p_temp)
      #       p_temp=p_temp^(10*iterations)
      #       my_zeros=which((p_temp<0)==TRUE)
      #       p_temp[my_zeros]=0
      #       if(sum(p_temp)==0){
      #         idx_max_cell_prob[i]=as[i]
      #       }else{
      #         idx_max_cell_prob[i]=sample(expected_clusters,1,TRUE,p_temp)
      #       }
      #     }
      #   )
      # }
      #----------------------------
      temp_mRNA_all=mRNA_all+1
      opt_idx_clusters_list=unlist(lapply(1:nrow(opt_idx_clusters), function(j) opt_idx_clusters[j,]) )
      opt_param_each_gene=array(log_p_mat_ma[,opt_idx_clusters_list],dim=c(max_mRNA_counts+1,n_genes,num_clusters))
      for(j in 1:n_genes){
        cell_prob=cell_prob+opt_param_each_gene[temp_mRNA_all[,j],j,]
      }
      
      switch (algorithm,
          'greedy_cabsel' = {
            sum_prob=apply(cell_prob,1,max)
            idx_max_cell_prob=apply(cell_prob,1,function(x){
              which(x==max(x))[1]
            })  #alway notice that, not just one max value
            #idx_max_cell_prob=matrix(0,nrow = 1,ncol = nvars)
            #idx_max_cell_prob=matrix(idx_max_cell_prob,nrow = 1,ncol = nvars)
            sum_prob_tot[jj]=sum(sum_prob) #this is different as matlab

            maxrel=1
          },
          'greedy_maxdiff'={
            tmp=cell_prob-matrix(rep(rowSums(cell_prob),expected_clusters),ncol = expected_clusters,byrow = FALSE)
            sum_prob=apply(tmp,1,max)
            idx_max_cell_prob=apply(tmp,1,function(x){
              which(x==max(x))[1]
            }) 
            #idx_max_cell_prob=matrix(idx_max_cell_prob,nrow = 1,ncol = nvars)
            sum_prob_tot[jj]=sum_prob_tot[jj]+sum_prob
            maxrel=1
          }
        )
      
      
      opt_when_false$sum_pro_tot_jj=sum_prob_tot[jj]
      comp=(as==idx_max_cell_prob)
      rel=sum(comp)
      #print('ready into maxrel')
      #print(as)
      #print('rel/nvars')
      #print(rel/nvars)
      #print(as)
      #print(idx_max_cell_prob)
      opt_when_false$population_jj=population[jj,]
      if(rel/nvars>=maxrel || optimize==FALSE){
        #print('rel/nvars>=maxrel')
        #print(rel/nvars>=maxrel)
        #print('rel/nars')
        #print(as)
        opt_when_false$population_jj=as
        #print('pulupation[jj,]')
        #print(population[jj,])
        my_break=TRUE
        if(optimize==FALSE){
          for (m in 1:expected_clusters) {
            aaa=which(as==m)
            my_distance[aaa,]=abs(rep(cell_prob[aaa,m],1*expected_clusters)-cell_prob[aaa,])   #check heres
          }
          clusterprobabilities=2^cell_prob;
          clusterprobabilities=clusterprobabilities/rep(rowSums(clusterprobabilities),1*expected_clusters)
          opt_when_false$clusterprobabilities=clusterprobabilities
          opt_when_false$distance=my_distance
          opt_when_false$cell_prob=cell_prob
          opt_when_false$best=as
          # print('my_results.best')
          # print(my_results[[1]]$best)  #it seem can not assigntm in foreach
        }
      }
      as=idx_max_cell_prob
    }
    #print(paste(toString(jj),'. loop finished'))
    #opt_idx_clusters_temp
    #sink()
    opt_when_false
  }
  #stopImplicitCluster(ncore)
  #stopImplicitCluster()
  stopCluster(ncore)
  if(loops>1){
    for (ii in 1:loops) {
      opt_idx_a[[ii]]$run=opt_idx_a_temp[[ii]]$opt_idx_clusters_temp
      population[ii,]=opt_idx_a_temp[[ii]]$population_jj
      sum_prob_tot[ii]=opt_idx_a_temp[[ii]]$sum_pro_tot_jj
    }
  }else
  {
    opt_idx_a[[1]]$run=opt_idx_a_temp[[1]]$opt_idx_clusters_temp
    population=opt_idx_a_temp[[1]]$population_jj
    my_results$clusterprobabilities=opt_idx_a_temp[[1]]$clusterprobabilities
    my_results$distance=opt_idx_a_temp[[1]]$distance
    my_results$cell_prob=opt_idx_a_temp[[1]]$cell_prob
    my_results$best=opt_idx_a_temp[[1]]$best
    
  }
    optpar=list()
    if(optimize==TRUE){
      indx_best=which(sum_prob_tot==max(sum_prob_tot))[1]
      best_as=population[indx_best,]
      my_results$population=population
      my_results$best=best_as
      my_results$overallsum=sum_prob_tot
      #print(indx_best)
      opt_id=opt_idx_a[[indx_best]]$run
      for (j in 1:nrow(opt_id)) {
        optpar[[j]]=K_new[opt_id[j,],]
      }
      my_results$parameter=optpar
    }else{
      my_results$population=NULL
      opt_id=opt_idx_a[[1]]$run  #here ,watch out
      #optpar=matrix(0,nrow = nrow(opt_id),ncol = ncol(K_new))
      for(j in 1:nrow(opt_id)){
        optpar[[j]]=K_new[opt_id[j,],]
      }
      my_results$parameter=optpar
    }
  return(my_results)

 }