#' A calista Function
#'
#' @param DATA,Results,cell_assignments
#' @keywords calista
#' @export
#' @examples
#'
#'
#'  CALISTA_transition_genes_main()


CALISTA_transition_genes_main<-function(DATA,INPUTS,Results,cell_assignments){
    if(nargs()<3){
    stop('Not enough input variables.')
  }

  if (length(Results)==0){
    if (nargs()<4){
      stop("Not enough input arguments")
    }
    Results=jump_clustering(DATA,cell_assignments)
    Results=jump_transition(DATA,Results)

    ###3.b- plot mean expression for each cluster
    Results=Plot_cluster_mean_exp(Results,DATA)

    ###3.c- cell-celll variability analysis
    Results=cell_variability(Results,DATA)
  }

  writeLines('\nCALISTA_transition_genes is running...\n')
  Parameters=DATA$Parameters;
  thr=INPUTS$thr_transition_genes
  my_results_final=Results$clustering_struct
  hh=Results$TRANSITION$final_graph
  nodes_connection2=get.edgelist(hh)
  nodes_connection2=nodes_connection2[order(nodes_connection2[,1]),]
  # if(nrow(nodes_connection2)>3){
  #   numSubplot_tList=numSubplot(nrow(nodes_connection2))
  #   p=numSubplot_tList$p
  # }else{
  #   p=integer()
  #   p[1]=1
  #   p[2]=nrow(nodes_connection2)
  # }
  sorted_gene_prob=list()
  idx_transition_genes=list()
  transition_gene_ranking=list()
  num_transition_genes=list()
  final_transition_genes=list()
  null_LL=list()
  # x11(title = 'CALISTA_transition_genes')
  #par(mfrow=c(p[1],p[2]))
  for(i in 1:nrow(nodes_connection2)){
    count=1
    temp_mRNA_all=NULL
    selected_clusters=NULL
    for(j in 1:ncol(nodes_connection2)){
      idx_final_group_temp=which(Results$final_groups==nodes_connection2[i,j])
      temp_mRNA_all=rbind(temp_mRNA_all,DATA$totDATA[idx_final_group_temp,])
      selected_clusters=cbind(selected_clusters,count*matrix(1,1,length(idx_final_group_temp)))
      count=count+1
    }
    nvar_temp=nrow(temp_mRNA_all)
    prob1=get_prob_transition_genes(selected_clusters,Parameters$Parameters[[3]],temp_mRNA_all,DATA$numGENES,nvar_temp)
    prob2=get_prob_transition_genes(matrix(1,1,nvar_temp),Parameters$Parameters[[3]],temp_mRNA_all,DATA$numGENES,nvar_temp)
    prob_separated_clusters=colSums(prob1)
    prob_all_cells=prob2
    sorted_gene_prob[[i]]=sort(prob_separated_clusters-prob_all_cells,decreasing = TRUE)
    idx_transition_genes[[i]]=order(prob_separated_clusters-prob_all_cells,decreasing = TRUE)
    transition_gene_ranking[[i]]=DATA$genes[idx_transition_genes[[i]]]
    null_LL[[i]]=sum(prob2)
    nodes_connection2[i,]=sort(nodes_connection2[i,])
    cum_LL=cumsum(sorted_gene_prob[[i]])
    idx_thr_tmp=which(cum_LL<=(thr*max(cum_LL)/100))
    idx_thr=idx_thr_tmp[length(idx_thr_tmp)]
    #print(idx_thr)
    num_transition_genes[[i]]=idx_thr
    final_transition_genes[[i]]=transition_gene_ranking[[i]][1:num_transition_genes[[i]]]
    Results$GENES$sorted_prob[[i]]=sorted_gene_prob[[i]][1:num_transition_genes[[i]]]
    #print(sorted_gene_prob[[i]][1:num_transition_genes[[i]]])
    # barplot(sorted_gene_prob[[i]][1:num_transition_genes[[i]]],
    #         main = paste(nodes_connection2[i,1],nodes_connection2[i,2],sep = '-'),
    #         col = 'blue',ylab = 'logP',
    #         names.arg = final_transition_genes[[i]],
    #         #cex.names = 0.8,
    #         las=2
    # )
    # legend('topright',legend = paste(num_transition_genes[[i]],' TGenes'))

  }
  idx_final_transition_genes=list()
  transition_gene_parameters=list()
  for(transition_num in 1:nrow(nodes_connection2)){
    idx_final_transition_genes[[transition_num]]=match(final_transition_genes[[transition_num]],DATA$genes)
    array_xDim=length(my_results_final$all$all$parameter[[1]][1,])
    transition_gene_parameters[[transition_num]]=array(0,c(array_xDim,ncol(nodes_connection2),num_transition_genes[[transition_num]]))
    for(i in 1:num_transition_genes[[transition_num]]){
      for(j in 1:ncol(nodes_connection2)){
        transition_gene_parameters[[transition_num]][,j,i]=t(my_results_final$all$all$parameter[[nodes_connection2[transition_num,j]]][idx_final_transition_genes[[transition_num]][i],])
      }
    }
  }
  idx_tot_transition_genes=sort(unique(unlist(idx_final_transition_genes)))
  tot_transition_genes=DATA$genes[idx_tot_transition_genes]
  Results$GENES$mRNA_tot_transition_genes=DATA$totDATA[,idx_tot_transition_genes]
  Results$GENES$tot_transition_genes=tot_transition_genes
  Results$GENES$thr=thr
  Results$GENES$final_transition_genes=final_transition_genes
  Results$GENES$transition_gene_parameters=transition_gene_parameters

  return(Results)

}
