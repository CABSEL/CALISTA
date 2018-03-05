AddGraph<-function(nodes,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,
                   normed_pseudotime,score3,legendInfo,az,el,...){
  ###test
  #AddGraph(nodes,Results$colorMARK_calista,x_center,y_center,z_center,Results$final_groups,Results$expected_clusters,
  #         Results$cell_cluster_progression,Results$score3,Results$legendInfo_calista,az,el)
  # colorMARK2=Results$colorMARK_calista
  # ClusterGroup2=Results$final_groups
  # expected_clusters=Results$expected_clusters
  # normed_pseudotime=Results$cell_cluster_progression
  # score3=Results$score3
  # legendInfo=Results$legendInfo_calista
  #first form of calling 
  #AddGraph(nodes,Results$colorMARK_caession,Results$score3,Results$legendInfo_calista,az,el)
  varargin=list(...)
  if (length(varargin)>0){   #check out lista,x_center,y_center,z_center,Results$final_groups,Results$expected_clusters,
  #         Results$cell_cluster_progrhere, when coming into further code 
    h=varargin[[1]]
    NumberOfEdges=nrow(get.edgelist(h))
  }else{
    h=NULL
    non=length(unique(as.vector(nodes[,1:2])))
    h=graph.empty(n=non,directed = FALSE)
    h_begin_layout=matrix(c(1:non,rep(1,non)),ncol = 2)
    l=layout.reingold.tilford(h)
    l=h_begin_layout
    NumberOfEdges=1
    h=add_edges(h,c(nodes[1,1],nodes[1,2]),weight=nodes[1,3])
    #plot.igraph(h,layout=l)
  }
  xMinLim=min(normed_pseudotime)-1
  xMaxLim=max(normed_pseudotime)+1
  yMinLim=min(score3[,1])-1
  yMaxLim=max(score3[,1])+1
  zMinLim=min(score3[,2])-1
  zMaxLim=max(score3[,2])+1
  #x11()
  gg_data=data.frame()
  label_color=c()
  label_names=c()
  for(k in 1:expected_clusters){
    gg_data=rbind(gg_data,data.frame(x=as.vector(normed_pseudotime[ClusterGroup2==k]),
                       y=score3[ClusterGroup2==k,1],
                       z=score3[ClusterGroup2==k,2],
                       stage_colors=rep(colorMARK2[[k]],length(score3[ClusterGroup2==k,1]))))
                       #stages=rep(legendInfo[[k]],length(score3[ClusterGroup2==k,1]))
                       
    label_color=c(label_color,colorMARK2[[k]])
    label_names=c(label_names,legendInfo[[k]])
  }
  p<-ggplot(gg_data)
  #stage_color_sca=gg_data[,"stage_colors"]
  #names(stage_color_sca)=gg_data[,"stages"]
  p<-p+geom_point(mapping = aes(x=y,y=z),data = gg_data,color=gg_data[,"stage_colors"])
  p<-p+labs(x='COMP1',y='COMP2')
  #legend('topright',legend=legendInfo,pch = 16,col = colorMARK2)
  V(h)$size=1
  V(h)$color=colorMARK2
  E(h)$width=2.5
  h_after_layout=data.frame(x_center,y_center,z_center,label_color,label_names,stringsAsFactors = FALSE)
  colnames(h_after_layout)=c('x_center','y_center','z_center','label_color',"label_names")
  txt_lab=as.character(c(1:length(label_names)))
  p<-p+geom_point(mapping =aes(x=y_center,y=z_center),data = h_after_layout,size=10,color=h_after_layout[,"label_color"]) #+geom_text(label=txt_lab,data = h_after_layout)
  h_edge=get.edgelist(h)
   segment_data=data.frame(x=h_after_layout[h_edge[,1],2],
                           y=h_after_layout[h_edge[,1],3],
                           xend=h_after_layout[h_edge[,2],2],
                           yend=h_after_layout[h_edge[,2],3])
   p<-p+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                     size=1)
  plot(p)
  AddGraph_list=list()
  AddGraph_list$h=h
  AddGraph_list$i=NumberOfEdges
  AddGraph_list$p=p
  return(AddGraph_list)

}

