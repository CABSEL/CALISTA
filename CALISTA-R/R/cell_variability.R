cell_variability<-function(Results,DATA){
  writeLines('\nCalculating cell to cell variability...\n')
  idx_sorted_groups=order(Results$final_groups)
  mRNA_all=DATA$totDATA[idx_sorted_groups,]
  histcounts=hist(Results$final_groups,plot = FALSE)$counts
  cutDIMENSION=histcounts[which(histcounts!=0)]
  DDD=list()
  t_mRNA_all=t(mRNA_all)
  for (i in 1:length(cutDIMENSION)){
    DDD[[i]]=t_mRNA_all[,1:cutDIMENSION[i]]
    t_mRNA_all=t_mRNA_all[,-(1:cutDIMENSION[i])]
    
  }
  min_num_cells=10
  H=entropy_calculation(DDD,min_num_cells)
  x11(title = 'cell_variability')
  par(mfrow=c(2,2))
  boxplot.matrix(t(H),main='Boxplot for the entropy',
                 xlab='Cluster',ylab='Entropy')
  plot(colMeans(t(H)),type = "o",col='red',xlab = "",ylab = "")
  par(new=TRUE)
  plot(colMedians(t(H)),type = "o",col='Green',xlab = "cluster",ylab = "Entropy value")
  title('Entropy Mean and Median')
  legend('topright',legend = c('MeanEntropy','MedianEntropy'),
         col=c('red','green'),lty = 1)
  
  cell_cell_correlation=list()
  mean_cell_cell_correlation=numeric()
  for(i in 1:Results$expected_clusters){
    aaa=which(Results$final_groups==Results$cluster_predicted[i])
    cell_cell_correlation[[i]]=cor(t(DATA$totDATA[aaa,]))
    temp_corr=sort(cell_cell_correlation[[i]],decreasing = TRUE)
    temp_corr=temp_corr[(length(aaa)+1):length(temp_corr)]
    mean_cell_cell_correlation[i]=mean(abs(temp_corr))    
  }
  # plot(mean_cell_cell_correlation,type = 'o',
  #      xlab = 'cluster',ylab = 'correlation',
  #      main = 'Mean cell-cell correlation in each other',col='red')
  # Results$cluster_entropy=H
  # Results$mean_cell_cell_corelation=mean_cell_cell_correlation
  # 
  # plot(1:Results$expected_clusters,Results$mean_prob_in_out_cluster[,1],
  #      type = 'o',col='green',main = 'Mean Probability',
  #      xlab = 'cluster',ylab = 'Log P')

  return(Results)
  
}