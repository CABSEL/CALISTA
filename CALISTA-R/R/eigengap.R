eigengap<-function(consensus,max_clust){
  D=diag(colSums(consensus))
  L=D-consensus
  #Dtemp=expm(-1/2*logm(D))
  #Dtemp=D%^%(-1/2)
  E <- eigen(D) 
  V <- E$values 
  Q <- E$vectors 
  Dtemp <- Q%*%diag(1/sqrt(V))%*%t(Q) 
  L=Dtemp%*%L%*%Dtemp
  v=eigen(L)$value
  v=sort(v)
  x11()
  all_diff=get_diff(v[1:max_clust])
  ind_maxdiff=which(all_diff==max(all_diff))[1]
  plot(v[1:max_clust],col="blue",type="p",
       main="Eigengap value",
       xlab="Number of cluster",ylab="")
  points(ind_maxdiff,v[ind_maxdiff],col="red")
  writeLines(paste("Optimal number of cluster according to max. eigenvalue:",ind_maxdiff,"\nif you want to use this value press 0, else provide desired number of cluster:"))
  input1=as.integer(readLines(n=1))
  if(input1==0){
    expected_clusters=ind_maxdiff
  }else{
    expected_clusters=input1
  }
  
  return(expected_clusters)
}