numSubplot<-function(n){
  
  numSubplot_list=list()
  while(is_prime(n) & n>4){
    n=n+1
  }
  p=primeFactors(n)
  if (length(p)==1){
    numSubplot_list$p=p
    numSubplot_list$n=n
    return(c(1,p))
  }
  while(length(p)>2){
    if(length(p)>=4){
      p[1]=p[1]*p[length(p)-1]
      p[2]=p[2]*p[length(p)]
      p=p[1:(length(p)-2)]
    }else{
      p[1]=p[1]*p[2]
      p=p[-2]
    }
    p=sort(p)
  }
  
  while((p[2]/p[1])>2.5){
    N=n+1
    temp_list=numSubplot(N)
    p=temp_list$p
    n=temp_list$n
  }
  
  numSubplot_list$p=p
  numSubplot_list$n=n
  return(numSubplot_list)
}