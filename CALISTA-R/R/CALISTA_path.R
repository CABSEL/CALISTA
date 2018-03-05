CALISTA_path<-function(Results,INPUTS){
  plot_fig=INPUTS$plot_fig
  hclustering=INPUTS$hclustering
  CV_plot=INPUTS$CV_plot
  #dev.off(figure102)
  #x11(title = '102')
  #figure102=dev.cur()
  
  if (!"PATH" %in% names(Results)){
  
  hh=Results$TRANSITION$final_graph
  add_paths=TRUE
  count=1
  CELL_path=list()
  x11(title = '')
  figure101=dev.cur()
  while(add_paths){
    writeLines(paste("Path num: ",count))
    writeLines("\n*********************************\n")
    #x11(display = figure101)
    grid.arrange(Results$orderP1,Results$orderP2,Results$orderP3)
    writeLines("Please check the figure 'Cell Ordering.R' for reference and  ")
    writeLines('Key the clusters in the path (e.g. 1 2 3 4):')
    temp_path=strsplit(readLines(n=1), " ")
    temp_path=c(as.integer(temp_path[[1]]))
    CELL_path[[count]]=temp_path
    writeLines('Please press 1 to add another path, 0 otherwise:')
    continue_add=as.integer(readLines(n=1))
    if(continue_add){
      count=count+1
    }else{
      add_paths=FALSE
    }
  }
  } else {CELL_path=Results$PATH$CELL_path}

  n_paths=length(CELL_path)
  Results$num_cells_each_cluster=hist(Results$final_groups,(1:(Results$expected_clusters+1))-0.5,plot = FALSE)$counts
  path_mRNA_all=list()
  smoothExpr=list()
  smoothExpr_std=list()
  CV=list()
  windowCenters=list()
  clustergrams=list()
  path_transition_genes=list()
  gg_datas=list()
  cv_datas=list()
  for(index in 1:length(Results$GENES$tot_transition_genes)){
    gg_datas[[index]]=data.frame()
    cv_datas[[index]]=data.frame()
  }

  for(i in 1:n_paths){
    path_len=length(CELL_path[[i]])
    nodes_connections=cbind(CELL_path[[i]][1:path_len-1],CELL_path[[i]][2:path_len])
    idx_swap_nodes=which(nodes_connections[,1]>nodes_connections[,2])
    nodes_connections[idx_swap_nodes,]=cbind(nodes_connections[idx_swap_nodes,2],nodes_connections[idx_swap_nodes,1])
    #nodes_connections_vec=lapply(1:nrow(nodes_connections),function(i){return(nodes_connections[i,])})
    #graph_vec=lapply(1:nrow(get.edgelist(Results$TRANSITION$final_graph)),function(i){return(get.edgelist(Results$TRANSITION$final_graph)[i,])})
    graph_matrix=get.edgelist(Results$TRANSITION$final_graph)
    idx_actual_edges_in_path=sapply(1:nrow(nodes_connections),function(i){
      return_val=0
      for(j in 1:nrow(graph_matrix)){
       if(all(nodes_connections[i,]==graph_matrix[j,])){
         return_val=j
       }
      }
      return(return_val)
    })
    idx_cells_in_path=t(which(Results$final_groups %in% CELL_path[[i]]==1))
    #define cells in path
    cells_in_path=Results$ORDERING$cell_ordering[idx_cells_in_path,1]
    #define transition genes of the actural path
    actual_transition_genes=Results$GENES$final_transition_genes[idx_actual_edges_in_path]
    actual_transition_genes=unique(unlist(actual_transition_genes))
    idx_actual_transition_genes=match(actual_transition_genes,Results$GENES$tot_transition_genes)
    n_transition_genes=length(idx_actual_transition_genes)
    idx_sorted_cells=sort(cells_in_path)
    idxSORTED_cells_in_path=order(cells_in_path)
    path_mRNA_all[[i]]=Results$GENES$mRNA_tot_transition_genes[idx_cells_in_path,]
    path_mRNA_all[[i]]=path_mRNA_all[[i]][idxSORTED_cells_in_path,]
    radius_width=round(max(length(cells_in_path)*INPUTS$moving_average_window/100,10))
    smoothExpr[[i]]=apply(path_mRNA_all[[i]],2,function(x){
      runMean(x,n=radius_width)
    })
    smoothExpr[[i]]=as.matrix(na.omit(data.frame(smoothExpr[[i]])))
    if(CV_plot){
      smoothExpr_std[[i]]=apply(path_mRNA_all[[i]],2,function(x){
        runSD(x,n=radius_width)
      })
      smoothExpr_std[[i]]=as.matrix(na.omit(data.frame(smoothExpr_std[[i]])))
      CV[[i]]=smoothExpr_std[[i]]/smoothExpr[[i]]
      Results$PATH$CV[[i]]=CV[[i]]
    }
    windowCenters[[i]]=1:nrow(smoothExpr[[i]])
    if(plot_fig){
      #plot one figure with all genes
      p1=sqrt(length(Results$GENES$tot_transition_genes))
      colorMARK_smooth=character()
      colorMARK_smooth_matrix=data.frame()
      cell_in_clusters=Results$num_cells_each_cluster[CELL_path[[i]]]
      tmp_len=length(cell_in_clusters)
      if(tmp_len>=3){
        cell_in_clusters=c(cell_in_clusters[1],round(cell_in_clusters[2:(tmp_len-1)]/2),cell_in_clusters[tmp_len])
      }else{
        cell_in_clusters=c(cell_in_clusters[1],cell_in_clusters[tmp_len])
      }
      counts_cells_in_path_edges=rowSums(t(rbind(cell_in_clusters[1:tmp_len-1],cell_in_clusters[2:tmp_len])))
      for(j in 1:nrow(nodes_connections)){
        #colfunc <- colorRampPalette(c(Results$colorMARK_calista[nodes_connections[j,1]], Results$colorMARK_calista[nodes_connections[j,2]]),
        #                            space="rgb",interpolate="linear",bias=0.01)
        #col2rgb and rgb
        #col2rgb(Results$colorMARK_calista[nodes_connections[j,1]])
        #rgb(0,127/255,255/255)
        color1=Results$colorMARK_calista[nodes_connections[j,1]]
        color2=Results$colorMARK_calista[nodes_connections[j,2]]
        color1_rgb=t(col2rgb(color1))/255
        color2_rgb=t(col2rgb(color2))/255
        colorMARK_smooth_matrix=rbind(colorMARK_smooth_matrix,
                                      data.frame(R=seq(color1_rgb[1],color2_rgb[1],length.out= counts_cells_in_path_edges[[j]]),
                                                 G=seq(color1_rgb[2],color2_rgb[2],length.out= counts_cells_in_path_edges[[j]]),
                                                 B=seq(color1_rgb[3],color2_rgb[3],length.out= counts_cells_in_path_edges[[j]])))
      }
      colorMARK_smooth=rgb(colorMARK_smooth_matrix)
      idx_subsampling=round(seq(from = 1,to=length(cells_in_path),length.out =nrow(smoothExpr[[i]]) ))
      #colorMARK_smooth_matrix=colorMARK_smooth_matrix[idx_subsampling,]
      colorMARK_smooth=colorMARK_smooth[idx_subsampling]
      #x11()
      #plot(1:length(colorMARK_smooth),1:length(colorMARK_smooth),col=colorMARK_smooth)
      #par(mfrow=c(p1,p1))
      #p_ggplot=list()
      for(index in 1:length(Results$GENES$tot_transition_genes)){
        
        gg_datas[[index]]=rbind(gg_datas[[index]],data.frame(x=windowCenters[[i]],
                           y=smoothExpr[[i]][,index],
                           color=colorMARK_smooth))
        # p<-ggplot(gg_data)
        # p<-p+geom_point(mapping = aes(x=x,y=y,color=color),data = gg_data)
        # p<-p+labs(x='Cell ordering',y='Mean Expr')
        # p<-p+theme(legend.position = "none")
        # p<-p+ggtitle(Results$GENES$tot_transition_genes[index])
        # p_ggplot[[index]]=p
        #par(new=TRUE)
      }
      
      ###there is a big problem in smoothExpr function?????
      #plot(p_ggplot[[1]])
      #grid.arrange(grobs=p_ggplot)
      if(CV_plot){
        for(index in 1:length(Results$GENES$tot_transition_genes)){
          cv_datas[[index]]=rbind(cv_datas[[index]],data.frame(x=windowCenters[[i]],
                                                               y=CV[[i]][,index],
                                                               color=colorMARK_smooth))
        }
      }
    }
    if(hclustering){
      #x11=title("heatmap")
      #figureheatmap=dev.cur()
      writeLines("\nHierarchical clustering based on cell ordering...\n")
      color_hierarchical=Results$c[idx_cells_in_path,]
      ColumnLabelsColorValue=list()
      ColumnLabelsColorValue$Labels=as.character(idx_sorted_cells)
      ColumnLabelsColorValue$Colors=color_hierarchical[idxSORTED_cells_in_path]
      transition_expression=Results$GENES$mRNA_tot_transition_genes[idx_cells_in_path,idx_actual_transition_genes]
      transition_expression=transition_expression[idxSORTED_cells_in_path,]
      transition_expression=scale(transition_expression)
      transition_expression=scale(t(transition_expression))
      rownames(transition_expression)=actual_transition_genes
      x11(title = "heatmap")
      figureheatmap=dev.cur()
      #distance= stats::dist(transition_expression, method ="euclidean")    
      #hcluster = hclust(distance, method ="average")
      distance_func= function(x){
        as.dist(1-cor(t(x), method="spearman"))
        }   
      hcluster_func= function(x){
        hclust(x, method ="average")
      }
      #hr <- hclust(as.dist(1-cor(t(transition_expression), method="spearman")), method="average")
      #hc <- hclust(as.dist(1-cor(transition_expression, method="spearman")), method="average") 
      hv=heatmap.2(transition_expression,
                      col = greenred(255),
                      #Rowv = FALSE,
                      Colv = FALSE,
                      labCol = NA,
                      distfun = distance_func,
                      hclustfun = hcluster_func,
                      dendrogram = "row",
                      key=TRUE,
                      density.info = 'none',
                      trace = 'none',
                      ColSideColors = ColumnLabelsColorValue$Colors, 
                      main = paste("Hierarchical clustering for path num:",i))
      
      
      clustergrams[[i]]=hv
    }
    #path_transition_genes[[i]]=actual_transition_genes[clustergrams[[i]]$rowInd]
    path_transition_genes[[i]]=actual_transition_genes
  }
  if(plot_fig){
    p_ggplot=list()
    p_ggplot1=list()
    for(index in 1:length(Results$GENES$tot_transition_genes)){
      p<-ggplot(gg_datas[[index]], aes(x=x,y=y))
      p<-p+geom_point(colour=gg_datas[[index]][,3])
      #p<-p+scale_color_gradient(low = color1,high=color2)
      p<-p+labs(x='Cell ordering',y='Mean Expr')
      p<-p+theme(legend.position = "none")
      p<-p+ggtitle(Results$GENES$tot_transition_genes[index])
      p_ggplot[[index]]=p
      
      if(CV_plot){
        p1<-ggplot(cv_datas[[index]],aes(x=x,y=y))
        p1<-p1+geom_point(colour=cv_datas[[index]][,3])
        p1<-p1+labs(x='Cell ordering',y='Mean Expr')
        p1<-p1+theme(legend.position = "none")
        p1<-p1+ggtitle(Results$GENES$tot_transition_genes[index])
        p_ggplot1[[index]]=p1
      }
    }
    ##there is a big problem in smoothExpr function?????
    #x11(title = "5000")
    #figure5000=dev.cur()
    ml=marrangeGrob(grobs = p_ggplot,ncol = 2,nrow = 2)
    ggsave("Mean expression profiles of transition genes.pdf",ml)
    writeLines('\n You can find "Mean expression profiles of transition genes "in "Mean expression profiles of transition genes.pdf" in current folder  \n')
    #grid.arrange(grobs=p_ggplot)
    
    if(CV_plot){
      # x11(title = "6000")
      # figure6000=dev.cur()
      # grid.arrange(grobs=p_ggplot1)  
      ml1=marrangeGrob(grobs = p_ggplot1,ncol = 2,nrow = 2)
      ggsave("CV_plot.pdf",ml)
      writeLines('\n You can find "CV_plot"in "CV_plot.pdf" in current folder\n')
    }
  }
  
  
  if(exists("clustergrams")){
    Results$PATH$hclustering=clustergrams
  }
  Results$PATH$CELL_path=CELL_path
  Results$PATH$smoothExpr=smoothExpr
  Results$PATH$windowCenters=windowCenters
  Results$PATH$path_transition_genes=path_transition_genes
  return(Results)
  
  
  
}