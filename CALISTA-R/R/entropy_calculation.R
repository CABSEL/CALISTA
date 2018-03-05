entropy_calculation<-function(DDD,min_num_cells){
  H=matrix(0,length(DDD),nrow(DDD[[1]]))
  for(i in 1:length(DDD)){
    n_genes=nrow(DDD[[i]])
    n_cells=ncol(DDD[[i]])
    for(j in 1:n_genes){
      x=DDD[[i]][j,]
      sorted_x=sort(x)
      in_bin=numeric()
      for(xx in 1:length(sorted_x)){
        if(xx==1){
          num_bin=1
          in_bin[[num_bin]]=0
        }
        in_bin[[num_bin]]=in_bin[num_bin]+1
        if(xx<length(sorted_x) & sorted_x[xx]!=sorted_x[xx+1]&in_bin[[num_bin]]>=min_num_cells){
          num_bin=num_bin+1
          in_bin[[num_bin]]=0
        }
      }
      sorted_x_1=sorted_x+1
      binning=list()
      for( jjj in 1:length(in_bin)){
        binning[[jjj]]=sorted_x_1[1:in_bin[[jjj]]]
        sorted_x_1=sorted_x_1[-(1:in_bin[[jjj]])]
      }
      min_mRNA_in_bin=1
      in_bin2=in_bin
      for( equiBIN in 1:length(in_bin)){
        max_mRNA_in_bin=max(binning[[equiBIN]])
        in_bin2[[equiBIN]]=length(min_mRNA_in_bin:max_mRNA_in_bin)
        min_mRNA_in_bin=max_mRNA_in_bin+1
        
      }
      Frequency=in_bin/n_cells
      bin_size=in_bin2
      H[i,j]=-sum(Frequency*log2(Frequency/bin_size))
    }
  }
  return(H)
}