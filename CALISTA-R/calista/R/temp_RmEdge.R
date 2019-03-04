#' A calista Function
#'
#' @param h,node,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,az,el
#' @keywords calista
#' @export
#' @examples
#' temp_RmEdge()


temp_RmEdge<-function(h,node,colorMARK2,x_center,y_center,z_center,ClusterGroup2,expected_clusters,normed_pseudotime,score3,legendInfo,az,el){

  gg_data=data.frame()
  label_color=c()
  label_names=c()
  for(k in 1:expected_clusters){
    gg_data=rbind(gg_data,data.frame(x=as.vector(normed_pseudotime[ClusterGroup2==k]),
                                     y=score3[ClusterGroup2==k,1],
                                     z=score3[ClusterGroup2==k,2],
                                     stage_colors=rep(colorMARK2[[k]],length(score3[ClusterGroup2==k,1]))
                                     #stages=rep(legendInfo[[k]],length(score3[ClusterGroup2==k,1]))
                                     ))
    label_names=c(label_names,legendInfo[[k]])
    label_color=c(label_color,colorMARK2[[k]])
  }
  h_after_layout=data.frame(x_center,y_center,z_center,label_color,stringsAsFactors = FALSE)
  colnames(h_after_layout)=c('x','y','z','stage_colors')
  txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
  gg_data=rbind(gg_data,h_after_layout)
  p<-ggplot(gg_data,mapping = aes(x=y,y=z))
  p<-p+geom_point(data = gg_data,color=gg_data[,"stage_colors"])+geom_text(aes(label=txt_lab),hjust=0,vjust=0,size=20)
  p<-p+labs(x='COMP1',y='COMP2')
  #legend('topright',legend=legendInfo,pch = 16,col = colorMARK2)
  h<-delete.edges(h,get.edge.ids(h,c(node[1],node[2])))
  V(h)$size=1
  V(h)$color=colorMARK2
  E(h)$width=2.5
  p<-p+geom_point(mapping =aes(x=y_center,y=z_center),data = h_after_layout,size=10,color=h_after_layout[,"stage_colors"])
  h_edge=get.edgelist(h)
  segment_data=data.frame(x=h_after_layout[h_edge[,1],2],
                          y=h_after_layout[h_edge[,1],3],
                          xend=h_after_layout[h_edge[,2],2],
                          yend=h_after_layout[h_edge[,2],3])
  p<-p+geom_segment(data = segment_data,aes(x=x,y=y,xend=xend,yend=yend),
                    size=1)
  plot(p)
  return(h)














}
