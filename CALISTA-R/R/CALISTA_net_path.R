CALISTA_net_path<-function(Results,INPUTS){
  # % CALISTA_net_path function infers the co-expression network for each
  # % developmental path detected previously
  # % 
  # % Usage:
  #   % 
  # % Results=CALISTA_net_path_main(Results,INPUTS)
  # % 
  # % Please refer to the README file for further details
  # % 
  # % Created by Nan Papili Gao
  # %            Institute for Chemical and Bioengineering 
  # %            ETH Zurich
  # %            E-mail:  nanp@ethz.ch
  # %
  # % Copyright. June 1, 2017.
  if(nargs()<2){
    stop("Not enought input arguments")
  }
  method=INPUTS$method
  switch(method,
         {
           value_cutoff=0.4
           pvalue_cutoff=0.05
         },
         {
           value_cutoff=0.8
           pvalue_cutoff=0.01
         })
  writeLines("\nCALISTA_net_path is running...\n")
  Theta=list()
  switch(method,
         {
           for(k in 1:length(Results$PATH$CELL_path)){
             pop=Results$PATH$smoothExpr[[k]]
             par_cor=pcor(pop)
             bbb=par_cor$estimate
             pvalue=par_cor$p.value
             Theta[[k]]=bbb*((pvalue<=pvalue_cutoff)&(abs(bbb)>=value_cutoff))
             Theta[[k]][is.na(Theta[[k]])]=0
           }
           str="Partial correlation"
         },
         {
           for(k in 1:length(Results$PATH$CELL_path)){
             pop=Results$PATH$smoothExpr[[k]]
             par_cor=rcorr(pop,type = "pearson")
             bbb=par_cor$r
             pvalue=par_cor$P
             Theta[[k]]=bbb*((pvalue<=pvalue_cutoff)&(abs(bbb)>=value_cutoff))
             Theta[[k]][is.na(Theta[[k]])]=0
           }
           str="Correlation"
         })
  #x11()
  #p1=sqrt(length(Results$PATH$CELL_path))+1
  #par(mfrow=c(p1,p1))
  for(i in 1:length(Results$PATH$CELL_path)){
    yvec=1:length(Results$GENES$tot_transition_genes)
    xvec=yvec
    xLabel=Results$GENES$tot_transition_genes
    yLabel=Results$GENES$tot_transition_genes
    x11()
    colMap <- colorRampPalette(c("red","blue","yellow" ))(nrow(Theta[[i]]*ncol(Theta[[i]])))
    image.plot(x=xvec,
          y=yvec,
          z=Theta[[i]],
          zlim=c(-1,1),
          useRaster=TRUE,
          xlab="",
          ylab="",
          main=paste("Path num", i),
          axes=FALSE,
          col=colMap
          #col=rainbow(nrow(Theta[[i]]*ncol(Theta[[i]])))
          #col=topo.colors(nrow(Theta[[i]]*ncol(Theta[[i]])))
          #col=terrain.colors(nrow(Theta[[i]]*ncol(Theta[[i]])))
          )
    grid (length(xvec),length(yvec), lty = 1, col = "black") 
    axis(1,at=xvec,labels = xLabel,tick = FALSE,las=2,font = 0.5)
    axis(2,at=yvec,labels = yLabel,tick = FALSE,las=1,font = 0.5)
    title(xlab="Target Gene j",ylab="Source Gene i",line = 4.8,cex.lab=1.5)
    #legend(grconvertX(0.5, "device"), grconvertY(1, "device"),
    #       c("0",".5","1"), fill = colMap[c(1, 10, 20)], xpd = NA)
    #abline (v=length(xvec),h=length(yvec), lty = 2, col = "black") 
    #plot_settings(xvec,yvec,Results$GENES$tot_transition_genes,Results$GENES$tot_transition_genes)
  }
  Results$PATH$Theta=Theta
  return(Results)
}