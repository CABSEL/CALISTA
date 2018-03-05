# find_progression2(final_groups,nt,method,my_stages);
find_progression2<-function(Results,DATA,method){
  
  final_groups=Results$final_groups
  nvars=length(final_groups)
  final_num_clusters=Results$expected_clusters
  cluster_progression=numeric(length = final_num_clusters)
  timeline=DATA$timeline
  final_groups_progression=list()
  legendInfo=list()
  if(length(unique(DATA$timeline))==1 && unique(DATA$timeline==0)){
    if(is.null(Results$starting_method)||Results$starting_method!=1||Results$starting_method!=2){
      writeLines('\nNo time info found. Please enter the starting cell or the marker gene whenever available\n')
      writeLines('Press 1 to enter the starting cell, 2 to enter the marker gene, 3 otherwise: ')
      Results$starting_method=readLines(n=1)
      Results$starting_method=as.integer(Results$starting_method)
    }
    switch (Results$starting_method,
            {
              if(is.null(Results$starting_cell)){
                writeLines('\n Key the index of the starting cell(e.g. 1):')
                Results$starting_cell=readLines(n=1)
                Results$starting_cell=as.integer(Results$starting_cell)
              }
              if(Results$starting_cell>DATA$nvars || Results$starting_cell<1){
                Results$starting_cell=NULL
                stop("Invaidl index, Please re-run the section and key a correct staring cell.")
              }else{
                final_groups2=t(final_groups)
                if(final_groups[Results$starting_cell]!=1){
                  final_groups2[which(final_groups==final_groups[Results$starting_cell])]=1
                  final_groups2[which(final_groups==1)]=final_groups[Results$starting_cell]
                }
              }
            },
            {
              if(is.null(Results$markerGENE)){
                writeLines('\n Enter the marker gene chosen(e.g. Erg, case insensitive):')
                Results$markerGENE=readLines(n=1)
              }
              if(is.null(Results$UPorDOWN)){
                writeLines('\npress 1 if it is upregulated, 2 if it is downregulated (e.g.2): ')
                Results$UPorDOWN=readLines(n=1)
                Results$UPorDOWN=as.integer(Results$UPorDOWN)
              }
              tf=tolower(Results$markerGENE) %in% tolower(DATA$genes)
              idx_markerGENE=match(tolower(Results$markerGENE), tolower(DATA$genes))
              if(tf==0){
                Results$markerGENE=NULL
                stop("Marker gene not found, please re-run the section and enter a correct gene name")
              }else{
                mean_markerGENE=numeric()
                for(i in 1:Results$expected_clusters){
                  mean_markerGENE[i]=mean(DATA$totDATA[which(Results$final_groups==i),idx_markerGENE])
                }
                if(Results$UPorDOWN==1){
                  initial_cluster=which(mean_markerGENE==min(mean_markerGENE))[1]
                }else if(Results$UPorDOWN){
                  initial_cluster=which(mean_markerGENE==max(mean_markerGENE))[1]
                }else{
                  Results$UPorDown=NULL
                  stop("Please re-run the section and enter 1 for up or 2 for down regulated marked gene")
                }
                final_groups2=t(final_groups)
                if(initial_cluster!=1){
                  final_groups2[which(final_groups==initial_cluster)]=1
                  final_groups2[which(final_groups==1)]=initial_cluster
                }
              }
            },
            {
              writeLines("\nSkip relabeling step\n")
              final_groups2=t(final_groups)
            }
    )
    stages=1:final_num_clusters
    Results$cluster_progression=(stages-min(stages))/(max(stages)-min(stages)) # normalized
    stt1='Cluster  '
    stt2='pseudo-stage'
    for (clust in 1:final_num_clusters) {
      idx_final_groups=which(final_groups2==clust)
      final_groups_progression[idx_final_groups]=stages[clust]
      legendInfo[clust]=paste(stt1,clust)
    }
  }else{
    MEAN_progression=matrix(0,1,final_num_clusters)
    SD_progression=matrix(0,1,final_num_clusters)
    for (clust in 1:final_num_clusters) {
      idx_final_groups=which(final_groups==clust)
      switch (method,
              {cluster_progression[clust]=median(timeline[idx_final_groups])},
              {cluster_progression[clust]=mean(timeline[idx_final_groups])},
              {cluster_progression[clust]=Mode(timeline[idx_final_groups])}
      )
      MEAN_progression[clust]=mean(timeline[idx_final_groups])
      SD_progression[clust]=sd(timeline[idx_final_groups])  #watch out here
    }
    stages=sort(cluster_progression)
    Results$cluster_progression=(stages-min(stages))/(max(stages)-min(stages)) # normalized
    idx_progression=order(cluster_progression)
    final_groups2=matrix(0,nvars,1)
    stt1='Cluster:'
    stt2='Cluster pseudotime:'
    for (clust in 1:final_num_clusters) {
      new_clust_order=which(idx_progression==clust)
      idx_final_groups=which(final_groups==clust)
      final_groups2[idx_final_groups,]=new_clust_order
      final_groups_progression[idx_final_groups]=stages[clust]
      legendInfo[clust]=paste(stt1,clust,stt2,format(round(Results$cluster_progression[clust], 2), nsmall = 2))
    }
  }

  Results$final_groups=t(final_groups2)
  Results$legendInfo_calista=legendInfo
  Results$final_groups_progression=final_groups_progression
  Results$cluster_predicted=sort(unique(Results$final_groups[1,]))
  Results$cell_cluster_progression=matrix(0,length(Results$final_groups),1)
  DATA$my_stages=t(final_groups_progression)
  for (i in 1:Results$expected_clusters){
    Results$cell_cluster_progression[which(Results$final_groups==Results$cluster_predicted[i])]=Results$cluster_progression[i]
  }
  #Results$cell_cluster_progression=(Results$cell_cluster_progression-min(Results$cell_cluster_progression))/max(Results$cell_cluster_progression-min(Results$cell_cluster_progression))
  find_progression2_list=list()
  find_progression2_list$Results=Results
  find_progression2_list$DATA=DATA
  return(find_progression2_list)
}
























