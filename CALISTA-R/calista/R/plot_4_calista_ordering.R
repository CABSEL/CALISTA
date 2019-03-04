#' A calista Function
#'
#' @param groups,colorMARK,group_names,Results
#' @keywords calista
#' @export
#' @examples
#' plot_4_calista_ordering()




plot_4_calista_ordering<-function(groups,colorMARK, group_names, Results){


  if (!is.null(nrow(groups)) && nrow(groups)<ncol(groups)){
    groups=t(groups)
  }

  unique_groups=unique(groups)
  Xcoordinate=Results$ORDERING$Xcoordinate
  Ycoordinate=Results$ORDERING$Ycoordinate
  Zcoordinate=Results$ORDERING$Zcoordinate
  x_center=Results$ORDERING$x_center
  y_center=Results$TRANSITION$y_center
  z_center=Results$TRANSITION$z_center
  h=Results$TRANSITION$final_graph
  h_edge=get.edgelist(h)

  if (length(colorMARK)==length(groups)){

    p=plot_ly(x=Xcoordinate,y=Ycoordinate,z=Zcoordinate,  color= colorMARK) %>%
      add_trace(x = ~Xcoordinate, y = ~Ycoordinate, z = ~Zcoordinate,  type = 'scatter3d',
                mode = 'markers',size=1, alpha=0.8) %>%
      layout(title = 'Cell Ordering',
             scene = list(camera = list(eye = list(x = -1.7, y = -1.7, z = 0.5)),
                          xaxis = list(title = 'Pseudotime',
                                       gridcolor = 'rgb(255, 255, 255)',
                                       #range = c(xMinLim,xMaxLim),
                                       #type = 'log',
                                       #zerolinewidth = 1,
                                       #ticklen = 1,
                                       nticks=10,
                                       gridwidth = 3),

                          yaxis = list(title = 'COMP1',
                                       gridcolor = 'rgb(255, 255, 255)',
                                       #range = c(yMinLim,yMaxLim),
                                       #zerolinewidth = 1,
                                       #ticklen = 1,
                                       nticks=10,
                                       gridwith = 3),
                          zaxis = list(title = 'COMP2',
                                       gridcolor = 'rgb(255, 255, 255)',
                                       #type = 'log',
                                       #range = c(zMinLim,zMaxLim),
                                       #zerolinewidth = 1,
                                       #ticklen = 1,
                                       nticks=10,
                                       gridwith = 3)),
             paper_bgcolor = 'rgb(243, 243, 243)',
             plot_bgcolor = 'rgb(243, 243, 243)')

  }else{


    gg_data=data.frame()
    label_color=c()
    label_names=c()
    for(k in 1:length(unique_groups)){
      gg_data=rbind(gg_data,data.frame(x=Xcoordinate[groups==unique_groups[k]],
                                       y=Ycoordinate[groups==unique_groups[k]],
                                       z=Zcoordinate[groups==unique_groups[k]],
                                       stage_colors=rep(colorMARK[k],length(Xcoordinate[groups==unique_groups[k]])),
                                       #stages=rep(Results$legendInfo_calista[[k]],length(Results$score3[groups==unique_groups[k],1])),
                                       gg_names=rep(paste(group_names,k),length(Xcoordinate[groups==unique_groups[k]]))
      ))
      label_color=c(label_color,colorMARK[k])
      label_names=c(label_names,paste(group_names,groups[k]))
    }
    names(label_color)=label_names
    h_after_layout=data.frame(x_center,y_center,z_center,stringsAsFactors = FALSE)
    colnames(h_after_layout)=c('x','y','z')
    txt_lab=as.character(c(rep("",nrow(gg_data)),1:length(label_names)))
    # gg_data=rbind(gg_data,h_after_layout)
    X=matrix(0,2,length(h_after_layout[h_edge[,1],1]))
    Y=X
    Z=X
    X[1,]=h_after_layout[h_edge[,1],1]
    X[2,]=h_after_layout[h_edge[,2],1]
    Y[1,]=h_after_layout[h_edge[,1],2]
    Y[2,]=h_after_layout[h_edge[,2],2]
    Z[1,]=h_after_layout[h_edge[,1],3]
    Z[2,]=h_after_layout[h_edge[,2],3]


    segment_data=data.frame(x=c(X),
                            y=c(Y),
                            z=c(Z))

    splitting=c(t(replicate(2, 1:(length(segment_data$x)/2))))
    ###########
    p <- plot_ly(gg_data) %>%
      add_trace(x = ~x, y = ~y, z = ~z, color = ~gg_names, colors = colorMARK, type = 'scatter3d',
                mode = 'markers',size=1, alpha=0.8) %>%
      add_trace(x = segment_data$x, y = segment_data$y,
                z=segment_data$z, split=splitting,name='Graph',type = 'scatter3d',
                mode = 'lines', opacity = 1, line = list(color = 'rgb(22, 96, 167)',width = 3,alpha=0.8),
                showlegend = FALSE)  %>%
      layout(title = 'Cell Ordering',
             scene = list(camera = list(eye = list(x = -0.1, y = -2.5, z = 0.9)),
                          xaxis = list(title = 'Pseudotime',
                                       gridcolor = 'rgb(255, 255, 255)',
                                       #range = c(xMinLim,xMaxLim),
                                       #type = 'log',
                                       #zerolinewidth = 1,
                                       #ticklen = 1,
                                       nticks=10,
                                       gridwidth = 3),

                          yaxis = list(title = 'COMP1',
                                       gridcolor = 'rgb(255, 255, 255)',
                                       #range = c(yMinLim,yMaxLim),
                                       #zerolinewidth = 1,
                                       #ticklen = 1,
                                       nticks=10,
                                       gridwith = 3),
                          zaxis = list(title = 'COMP2',
                                       gridcolor = 'rgb(255, 255, 255)',
                                       #type = 'log',
                                       #range = c(zMinLim,zMaxLim),
                                       #zerolinewidth = 1,
                                       #ticklen = 1,
                                       nticks=10,
                                       gridwith = 3)),
             paper_bgcolor = 'rgb(243, 243, 243)',
             plot_bgcolor = 'rgb(243, 243, 243)')

  }
  return (p)

}
