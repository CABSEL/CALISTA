#   %IMPORT_DATA upload single-cell expression data and perform preprocessing
#   % CALISTA accepts single-cell RTqPCR and RNA sequencing data. 
#   % Besides the expression data matrix, users can also provide capture time 
#   % or cell stage information. 
#   %
#   % Usage:
#   % DATA=import_data(INPUTS)
#   % Upload and preprocess data using user-defined specifications
#   % 
#   % Input:
#   % INPUTS - a structure containing data description and data preprocessing 
#   % settings
#   %
#   % ** INPUTS$data_type **
#   % 1- for scRT-qPCR  (CT values)
#   % 2- for scRT-qPCR  (Expression values - Num. mRNA molecules)
#   % 3- for scRNA-seq  (Expression values - e.g log(TPM+1) or log(RPKM+1))
#   % 4- for scRNA-seq  (Expression values - e.g TPM+1 or RPKM)
#   %
#   % ** INPUTS$format_data **
#   % Five data formats are accepted.
#   %
#   % 1- Rows= cells and Columns= genes with time/stage info in the last column  
#   % ------------------------------------------
#   % Gene1  Gene2  Gene3  ... Genej  Time/Stage
#   %  27     80     56    ...  69        0
#   %  73     20     90    ...  45        0
#   %   .      .      .    ...   .        .
#   %   .      .      .    ...   .        .
#   %   .      .      .    ...   .        .
#   % ------------------------------------------
#   %
#   % 2- Rows= genes and Columns= cells with time/stage info in the first row
#   % ----------------------------------------------
#   % Time/Stages   0      0      1   ...  16     16
#   %   Gene1      27     80     56   ...  69      0
#   %   Gene2      73     20     90   ...  45      0
#   %   Gene3       .      .      .   ...   .      .
#   %   Gene4       .      .      .   ...   .      .
#   %   Genej       .      .      .   ...   .      .
#   % ----------------------------------------------
#   % 
#   % 3- Rows= cells and Columns= genes without time info
#   % ------------------------------
#   % Gene1  Gene2  Gene3  ... Genej  
#   %  27     80     56    ...  69       
#   %  73     20     90    ...  45       
#   %   .      .      .    ...   .       
#   %   .      .      .    ...   .       
#   %   .      .      .    ...   .       
#   % ------------------------------
#   % 
#   % 4- Rows= genes and Columns= cells without time info
#   % --------------------------------------------
#   % Gene1      27     80     56   ...  69      0
#   % Gene2      73     20     90   ...  45      0
#   % Gene3       .      .      .   ...   .      .
#   % Gene4       .      .      .   ...   .      .
#   % Genej       .      .      .   ...   .      .
#   % --------------------------------------------
#   %
#   % 5- Manual selection from data table
#   %
#   % ** INPUTS$data_selection **
#   % When cells come from different capture times or stages, users can select 
#   % cells from specific time or stage for further analysis by CALISTA. 
#   % For instance, considering a dataset with cells taken from 4 time points, 
#   % snapshots can be selected, as follows:
#   % INPUTS$data_selection = integer() or c(1:4) for all time points/cells 
#   % INPUTS$data_selection = 1          for cells the first time point
#   % INPUTS$data_selection = c(1 3 4)    for cells 1st, 3rd and 4th time points
#   % 
#   % ** INPUTS$perczeros_genes **
#   % Users can exclude genes with a certain percentage of zeros. 
#   % INPUTS$perczeros_genes = 100 (recommended)
#   %
#   % ** INPUTS$perczeros_cells **
#   % Users can exclude cells with more than a certain percentage of zeros.
#   % INPUTS$perczeros_cells = 100 (recommended)
#   % 
#   % ** INPUTS$cells_2_cut **
#   % Users can exclude cells from further analysis. 
#   % 1- Remove cells based on their indices in the expression matrix. 
#   %    Indices need to be uploaded as a separate csv file. 
#   % 0- No cells deletion
#   % 
#   % ** INPUTS$perc_top_genes **
#   % Users can specify a certain percentage of genes to be used for CALISTA. 
#   % For computational efficiency, the number of genes is set to the
#   % smallest between the number of genes in the percentage above and 200. 
#   % INPUTS$perc_top_genes = 10 
#   %
#   % Output:
#   % DATA - a structure containing preprocessed expression values for further
#   % analysis in CALISTA
#   % 
#   % Created by Nan Papili Gao (R version implemented by Tao Fang)
#   %            Institute for Chemical and Bioengineering 
#   %            ETH Zurich
#   %            E-mail:  nanp@ethz.ch
#   %
#   % Copyright. June 1, 2017.
    

import_data <-function(INPUTS){
  print('**** Please upload normalized data. File formats accepted: .txt , .xlxs , .csv ****')
  if(nargs()<1){ 
    stop("Not enough input arguments")
    # INPUTS$data_type=1
    # INPUTS$format_data=1
    # INPUTS$cells_2_cut=0
    # INPUTS$data_selection=integer()
    # INPUTS$perczeros_genes=90
    # INPUTS$perczeros_cells=100
    # INPUTS$perc_top_genes=50
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
  Parameters=readMat('./R/Two-state model parameters/Parameters.mat')
  DATA$Parameters=Parameters
  return(DATA)
}






