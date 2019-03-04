#' A calista Function
#'
#' @param INPUTS
#' @keywords calista
#' @export
#' @examples
#'
#'
#' import_data()



import_data <-function(INPUTS){
  print('**** Please upload normalized data. File formats accepted: .txt , .xlxs , .csv ****')
  if(nargs()<1){
    stop("Not enough input arguments")

  }

  ## CHECK INPUT ARGUMENTS
  if (!exists('data_type', where=INPUTS)){
    stop('Please specify the data type in INPUTS.data_type')
  }

  if  (!exists('format_data', where=INPUTS)){
    stop('Please specify the format of the data in INPUTS$format_data')
  }

  # if  (!exists('cluster_time', where=INPUTS)){
  # INPUTS$cluster_time=0
  # }

  if  (!exists('cells_2_cut', where=INPUTS)){
    INPUTS$cells_2_cut=0;
  }

  if  (!exists('cut_variable_genes', where=INPUTS)){
    INPUTS$cut_variable_genes=1000
  }

  if(!exists('data_selection', where=INPUTS)){
    INPUTS$data_selection=integer()
  }

  #   if (INPUTS$cluster_time~=0 && ~INPUTS$data_selection==0){
  #   writelines('\nCALISTA Time CLustering is active. INPUTS$data_selection set as 0. All time points are considered for the analysis\n')
  #   INPUTS$data_selection=0
  # }

  if  (!exists('perczeros_genes', where=INPUTS)){
    INPUTS$perczeros_genes=100
  }

  if  (!exists('perczeros_cells', where=INPUTS)){
    INPUTS$perczeros_cells=100

  }

  if  (!exists('cells_2_cut', where=INPUTS)){
    INPUTS$cells_2_cut=0

  }

  if  (!exists('optimize', where=INPUTS)){
    INPUTS$optimize=0

  }

  if  (!exists('parallel', where=INPUTS)){
    INPUTS$parallel=1

  }

  if  (!exists('runs', where=INPUTS)){
    INPUTS$runs=50

  }

  if  (!exists('max_iter', where=INPUTS)){
    INPUTS$max_iter=100

  }


  if  (!exists('Cluster', where=INPUTS)){
    INPUTS$Cluster='kmedoids'

  }

  if  (!exists('thr_transition_genes', where=INPUTS)){
    INPUTS$thr_transition_genes=50

  }

  data_type=INPUTS$data_type
  format_data=INPUTS$format_data
  data_selection=INPUTS$data_selection
  perczeros_genes=INPUTS$perczeros_genes
  perczeros_cells=INPUTS$perczeros_cells
  cut_variable_genes=0
  cells_2_cut=INPUTS$cells_2_cut


  #uploading
  if(data_type>=1){

    # here normalization == old code.selection.load data at the beging(DATA.timelime,.time,.unm_time,.totDATA,.genes,.singeCELLDATA...)
    DATA=normalization(data_type,format_data,perczeros_genes,perczeros_cells,cut_variable_genes,cells_2_cut)
  }else{
    DATA=normalization(data_type,format_data,perczeros_genes,perczeros_cells,NULL,cells_2_cut=cells_2_cut)
  }

  ###check input arguments
  if(length(data_selection)==0){
    data_selection=c(1:length(DATA$singleCELLdata))
  }
  nc=length(data_selection)

  ###initialization
  timeline=numeric()
  genes=DATA$genes
  nt=list()
  for (i in 1:nc) {
    mRNA=t(DATA$singleCELLdata[[data_selection[i]]])
    nt[[i]]=nrow(mRNA)
    timeline=c(timeline,DATA$timeline[DATA$timeline==DATA$time[data_selection[[i]]]])
    if (i==1){
      mRNA_all=mRNA
    }else{
      mRNA_all=rbind(mRNA_all,mRNA)
    }
    if (length(data_selection)!=DATA$num_time_points){
      DATA$singleCELLdata[[i]]=t(mRNA)
    }
  }
  labels=as.character(timeline)
  my_stages=timeline
  DATA$totDATA=mRNA_all
  DATA$timeline=timeline
  if(length(data_selection)!=DATA$num_time_points){
    DATA$time=DATA$time[data_selection]
    DATA$num_time_points=length(DATA$time)
  }
  DATA$nc=nc
  DATA$nt=nt
  DATA$labels=labels
  DATA$my_stages
  DATA$nvars=nrow(DATA$totDATA)
  #selectt top genes
  if(data_type>=1){
    perc_top_genes=INPUTS$perc_top_genes
    n_top_genes=min(200,round(DATA$nvars*perc_top_genes/100),length(DATA$genes))
    if (DATA$numGENES>n_top_genes){
      DATA$genes=DATA$genes[1:n_top_genes]
      DATA$totDATA=DATA$totDATA[,1:n_top_genes]
      DATA$numGENES=length(DATA$genes)
      for(i in 1:DATA$num_time_points){
        DATA$singleCELLdata[[i]]=DATA$singleCELLdata[[i]][1:n_top_genes,]
      }
    }
  }
  ###parameters loading
  #Parameters=readMat('./Two-state model parameters/Parameters.mat')

  DATA$Parameters=calista:::Parameters
  return(DATA)
}






