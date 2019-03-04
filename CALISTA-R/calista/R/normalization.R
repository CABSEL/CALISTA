#' A calista Function
#'
#' @param data_type,format_data,perczeros_genes,perczeros_cells,cut_variable_genes,cells_2_cut
#' @keywords calista
#' @export
#' @examples
#' normalization()


normalization <-function(data_type,format_data,perczeros_genes,perczeros_cells,cut_variable_genes,cells_2_cut){
  filename=file.choose()
  # choose what kind of file? /home/Tao/Pro/AAA-CALISTA/RT-qPCR/6-data_WITH_VPA_0_24_48_72_CT_NORM.xlsx
  #imported_temp=read.xlsx(filename,1)  # in matlab code, imported_data.data is data, imported_data.text data is first line
  filetype=file_ext(filename)
  if(filetype=='xlsx'){
    imported_temp=read.xlsx(filename,1,rowNames = FALSE)
    imported_temp=na.omit(imported_temp)
  }else
  {
    if(format_data!=4){
      imported_temp=fread(filename,header = TRUE,
                          stringsAsFactors = TRUE,
                          data.table = FALSE,blank.lines.skip = TRUE)
      imported_temp=na.omit(imported_temp)
    }else{
      imported_temp=fread(filename,header = FALSE,
                          stringsAsFactors = TRUE,
                          data.table = FALSE,blank.lines.skip = TRUE)
      imported_temp=na.omit(imported_temp)
    }
  }
  #   if(filetype=='csv'){
  #   imported_temp=read.csv(filename,header=FALSE,as.is = FALSE,stringsAsFactors = FALSE)
  #   imported_temp=na.omit(imported_temp)
  # }else{
  #   imported_temp=read.table(filename,header = FALSE,stringsAsFactors = FALSE)
  #   imported_temp=na.omit(imported_temp)
  # }
  writeLines('\nData loading...\n')
  imported_data=list()
  DATA=list()
  DATA$FileName=filename
  switch (format_data,
          {
            #1
            #imported_data$data=as.matrix(unname(imported_temp[-1,]))
            #imported_data$textdata=imported_temp[1,]
            imported_data$data=as.matrix(unname(imported_temp))
            imported_data$textdata=colnames(imported_temp)
            NUM=imported_data$data
            NUM=NUM[!(is.nan(NUM[,1])),]
            totDATA=NUM[,1:ncol(NUM)-1]
            timeline=NUM[,ncol(NUM)]
            TXT=imported_data$textdata
            DATA$genes=TXT[1:length(TXT)-1]
            class(totDATA)='numeric'
            DATA$genes=unname(unlist(DATA$genes))
            class(timeline)='numeric'
            timeline=unname(timeline)
          },
          {
            #2

            imported_data$data=imported_temp[,-1]
            imported_data$textdata=imported_temp[,1]
            #imported_data$data=as.matrix(unname(imported_temp))
            #imported_data$textdata=colnames(imported_temp)
            time=as.numeric(colnames(imported_temp)[-1])
            #NUM=cbind(t(imported_data$data),time)
            #NUM=NUM[!(is.nan(NUM[,1])),]
            totDATA=t(imported_data$data)
            timeline=time
            TXT=imported_data$textdata
            DATA$genes=TXT
            class(totDATA)='numeric'
            DATA$genes=unname(unlist(DATA$genes))
            class(timeline)='numeric'
            timeline=unname(timeline)
          },
          {
            #3
            # imported_data$data=as.matrix(unname(imported_temp[-1,]))
            # imported_data$textdata=imported_temp[1,]
            imported_data$data=as.matrix(unname(imported_temp))
            imported_data$textdata=colnames(imported_temp)
            NUM=imported_data$data
            NUM=NUM[!(is.nan(NUM[,1])),]
            totDATA=NUM
            timeline=matrix(0,nrow(NUM),1)
            TXT=imported_data$textdata
            if(length(TXT)>ncol(totDATA)){
              DATA$genes=TXT[-1]
            }else{
              DATA$genes=TXT
            }
            class(totDATA)='numeric'
            DATA$genes=unname(unlist(DATA$genes))
            class(timeline)='numeric'
            timeline=unname(timeline)
          },
          {
            #4
            imported_data$data=t(as.matrix(unname(imported_temp[,-1])))
            imported_data$textdata=imported_temp[,1]
            #imported_data$data=as.matrix(unname(imported_temp))
            #imported_data$textdata=colnames(imported_temp)
            NUM=t(imported_data$data)
            NUM=NUM[!(is.nan(NUM[,1])),]
            totDATA=NUM
            timeline=matrix(0,nrow(NUM),1)
            TXT=imported_data$textdata
            if(length(TXT)>ncol(totDATA)){
              DATA$genes=TXT[-1]
            }else{
              DATA$genes=TXT
            }
            class(totDATA)='numeric'
            DATA$genes=unname(unlist(DATA$genes))
            class(timeline)='numeric'
            timeline=unname(timeline)
          },
          {
            #5
            writeLines("IMPORTED DATA preview:\n\n")
            print(imported_temp[1:5,1:10])
            writeLines(paste("Dimension of IMPORTED DATA:\n",
                             "Rows=",nrow(imported_temp),"\n",
                             "Columns=",ncol(imported_temp),"\n"))
            #head(imported_temp)
            writeLines("Press 1 if columns=genes, 0 otherwise (e.g. 1):")
            genes_vector=as.integer(readLines(n=1))
            if(genes_vector){
              writeLines("Key the starting column for expression data (e.g. 6):")
              start_genes=as.integer(readLines(n=1))
              writeLines("Selected EXPRESSION DATA preview:\n\n")
              print(imported_temp[1:nrow(imported_temp),start_genes:ncol(imported_temp)][1:5,1:10])
              #print(dim(imported_temp[1:nrow(imported_temp),start_genes:ncol(imported_temp)]))
              writeLines(paste("Dimension of selected EXPRESSION DATA:\n",
                               "Rows=",nrow(imported_temp[1:nrow(imported_temp),start_genes:ncol(imported_temp)]),"\n",
                               "Columns=",ncol(imported_temp[1:nrow(imported_temp),start_genes:ncol(imported_temp)]),"\n"))
              NUM=imported_temp[1:nrow(imported_temp),start_genes:ncol(imported_temp)]
              NUM=unname(as.matrix(NUM))
              TXT=colnames(imported_temp)
            }else{
              writeLines("key the starting row for expression data (e.g. 1):")
              start_genes=as.integer(readLines(n=1))
              print(imported_temp[start_genes:nrow(imported_temp),2:ncol(imported_temp)][1:5,1:10])
              print(dim(imported_temp[start_genes:nrow(imported_temp),2:ncol(imported_temp)]))
              NUM=imported_temp[start_genes:nrow(imported_temp),2:ncol(imported_temp)]
              TXT=imported_temp[,1]

            }
            writeLines("* Key starting and ending rows for the expression values (e.g. 1 405):*")
            rows=strsplit(readLines(n=1), " ")
            rows=c(as.integer(rows[[1]][1]),as.integer(rows[[1]][2]))
            writeLines("* Key starting and ending columns for the expression values (e.g. 1 22524):*")
            cols=strsplit(readLines(n=1), " ")
            cols=c(as.integer(cols[[1]][1]),as.integer(cols[[1]][2]))
            totDATA=NUM[rows[1]:rows[2],cols[1]:cols[2]]
            if(genes_vector){
              writeLines("*Key the starting and ending columns for gene names(e.g. 6 22529)*")
              cols=strsplit(readLines(n=1), " ")
              cols=c(as.integer(cols[[1]][1]),as.integer(cols[[1]][2]))
              DATA$genes=TXT[cols[1]:cols[2]]
            }else{
              writeLines("*Key the starting and ending rows for gene names(e.g. 6 22529)*")
              rows=strsplit(readLines(n=1), " ")
              rows=c(as.integer(rows[[1]][1]),as.integer(rows[[1]][2]))
              totDATA=t(totDATA)
              DATA$genes=TXT[rows[1]:rows[2]]
            }
            writeLines("*Add time info (1-Yes,0-No):*")
            time_info=as.integer(readLines(n=1))
            if (time_info){
              writeLines("*Key the column or row (e.g 5) vector in the raw data set with time/cell stage info:*")
              time_vector=as.integer(readLines(n=1))
              if(genes_vector){
                timeline=imported_temp[(rows[1]):(rows[2]),time_vector]
              }else{
                timeline=imported_temp[time_vector,(cols[1]):(cols[2])]
              }
            }else{
              timeline=matrix(0,nrow(totDATA),1)
            }
            class(totDATA)='numeric'
            DATA$genes=unname(unlist(DATA$genes))
            class(timeline)='numeric'
            timeline=unname(timeline)
          }
  )

  DATA$genes=as.character(DATA$genes)

  if(cells_2_cut==1){
    writeLines("***Select the csv file containing cell indeces to remove***\n\n")
    idxfilename=file.choose()
    cells_2_cut=read.table(idxfilename)  ## check the data format later???
    cells_2_cut=as.integer(unlist(cells_2_cut))
    timeline=timeline[-(cells_2_cut)]
    totDATA=totDATA[-(cells_2_cut),]
  }

  if (data_type==1){
    totDATA[totDATA>28]=28
    zeros_genes=colSums(totDATA==max(max(totDATA)))*100/nrow(totDATA)
    zeros_cells=rowSums(totDATA==max(max(totDATA)))*100/ncol(totDATA)
  }else{
    zeros_genes=colSums(totDATA==0)*100/nrow(totDATA)
    zeros_cells=rowSums(totDATA==0)*100/ncol(totDATA)}
  ### cutting genes
  idx2cut=which(zeros_genes>=perczeros_genes)
  if(length(idx2cut)!=0){
    totDATA=totDATA[,-idx2cut]
    DATA$genes=DATA$genes[-idx2cut]
    DATA$gene_removed=DATA$genes[idx2cut]
  }
  DATA$cut_sort$idx2cutGENES=idx2cut;
  ### cutting cells
  idx2cut=which(zeros_cells>=perczeros_cells)
  if(length(idx2cut)!=0){
    totDATA=totDATA[-idx2cut,]     # why same opereation for cuttint cell
    timeline=timeline[-idx2cut]
  }
  DATA$cut_sort$idx2cutCELL=idx2cut;



  ###Most variable genes
  if (data_type>=1){
    SD=apply(totDATA,2,sd)
    MEAN=colSums(totDATA)/colSums(totDATA!=0)
    CV=(SD^2)/MEAN
    nanCV=which(is.nan(CV))
    genes=DATA$genes
    if(length(nanCV)>=1){
      CV=CV[-(nanCV)]
      MEAN=MEAN[-(nanCV)]
      SD=SD[-(nanCV)]
      totDATA=totDATA[,-nanCV]
      genes=genes[-nanCV]
    }
    DATA$genes=genes
    nbin=21
    mean_hist=hist(x=MEAN,breaks = seq(min(MEAN),max(MEAN),length.out = nbin),right = FALSE,plot = FALSE)
    BIN=findInterval(MEAN,mean_hist$breaks,all.inside = TRUE)
    z_scoredCV=NULL
    most_variable_genes_idx=numeric()
    zscores_most_variable_genes=numeric()
    for(i in 1:(nbin-1)){
      if(length(which(BIN==i))>=0){
        idx_genes_each_bin=which(BIN==i)
        z_scored_CV=scale(CV[idx_genes_each_bin])
        # temp_genes=idx_genes_each_bin[which(abs(z_scored_CV)>=cut_variable_genes)]
        most_variable_genes_idx=c(most_variable_genes_idx,idx_genes_each_bin)
        # temp_zscores=z_scored_CV[which(abs(z_scored_CV)>=cut_variable_genes)]
        zscores_most_variable_genes=c(zscores_most_variable_genes,z_scored_CV)
      }
    }
    idx=order(abs(zscores_most_variable_genes),decreasing = TRUE)
    zscores_most_variable_genes=zscores_most_variable_genes[idx]
    most_variable_genes=genes[most_variable_genes_idx]
    most_variable_genes=most_variable_genes[idx]
    totDATA=totDATA[,most_variable_genes_idx]
    totDATA=totDATA[,idx]
    DATA$genes=most_variable_genes
    DATA$zscores_most_variable_genes=zscores_most_variable_genes

  }


  switch(data_type,
         {
         log_max_mRNA=log2(200)
         if(min(totDATA)<0){
           totDATA=totDATA-min(totDATA)
         }
         totDATA[totDATA>28]=28
         ctmax=max(totDATA)
         log2Ex=ctmax-totDATA
         base=2^(log_max_mRNA/max(log2Ex))
         totDATA=round(base^log2Ex)-1
         },
         {
           totDATA=log2(totDATA+1)  # here there is  a problem, r can't cal log(-int)
           log_max_mRNA=log2(200)
           totDATA[totDATA>28]=28
           if (min(totDATA)<0){
             totDATA=totDATA-min(totDATA)
           }
           ctmax=max(totDATA)
           log2Ex=totDATA
           tmp_mat=matrix(rep(apply(log2Ex, 2, max),nrow(log2Ex)),nrow = nrow(log2Ex),byrow = TRUE)
           exponent=log_max_mRNA*(log2Ex/tmp_mat)
           totDATA=round(2^exponent)-1
         },
         {
           totDATA=log2(totDATA+1)  # here there is  a problem, r can't cal log(-int)
           log_max_mRNA=log2(200)
           totDATA[totDATA>28]=28
           if (min(totDATA)<0){
             totDATA=totDATA-min(totDATA)
           }
           ctmax=max(totDATA)
           log2Ex=totDATA
           tmp_mat=matrix(rep(apply(log2Ex, 2, max),nrow(log2Ex)),nrow = nrow(log2Ex),byrow = TRUE)
           exponent=log_max_mRNA*(log2Ex/tmp_mat)
           totDATA=round(2^exponent)-1
         },
         {
           totDATA=log2(log2(totDATA+1)+1)  # here there is  a problem, r can't cal log(-int)
           log_max_mRNA=log2(200)
           totDATA[totDATA>28]=28
           if (min(totDATA)<0){
             totDATA=totDATA-min(totDATA)
           }
           ctmax=max(totDATA)
           log2Ex=totDATA
           tmp_mat=matrix(rep(apply(log2Ex, 2, max),nrow(log2Ex)),nrow = nrow(log2Ex),byrow = TRUE)
           exponent=log_max_mRNA*(log2Ex/tmp_mat)
           totDATA=round(2^exponent)-1
         }

  )

  DATA$timeline=timeline
  DATA$time=unique(timeline)
  DATA$num_time_points=length(DATA$time)
  sortTOTdata=numeric()
  sortTIMELINE=numeric()
  idx_sorted_cells=integer()
  cutDIMENSION=list()
  for (k in 1:DATA$num_time_points) {
    I=which(timeline==DATA$time[k])
    idx_sorted_cells=c(idx_sorted_cells,I)
    cutDIMENSION[k]=length(I)
    sortTOTdata=rbind(sortTOTdata,totDATA[I,])
    sortTIMELINE=c(sortTIMELINE,timeline[I])
  }
  cutDIMENSION=unlist(cutDIMENSION)
  DATA$cut_sort$idx_sorted_cells=idx_sorted_cells
  DATA$totDATA=sortTOTdata
  DATA$timeline=sortTIMELINE
  #####
  DATA$totDATA[is.nan(totDATA)]=0
  ###
  DATA$nvars=nrow(DATA$totDATA)
  DATA$numGENES=ncol(DATA$totDATA)


  data=list()
  data[[1]]=t(DATA$totDATA)[,1:cutDIMENSION[1]]
  if(length(cutDIMENSION)>1){
    for (i in 2:length(cutDIMENSION)) {
      data[[i]]=t(DATA$totDATA)[,(sum(cutDIMENSION[1:(i-1)])+1):sum(cutDIMENSION[1:i])]
    }
  }
  DATA$singleCELLdata=data
  DATA$imported_data=imported_temp
  return(DATA)
}














