function [DATA,INPUTS]=import_data(INPUTS)
%IMPORT_DATA upload single-cell expression data and perform preprocessing
% CALISTA accepts single-cell RTqPCR, RNA-seq, Drop-seq and DropNc-seq data. 
% Besides the expression data matrix, users can also provide capture time 
% or cell stage information. 
%
% Usage:
% [DATA,INPUTS]=import_data(INPUTS)
% Upload and preprocess data using user-defined specifications
% 
% Input:
% INPUTS - a structure containing data description and data preprocessing 
% settings
%
% ** INPUTS.data_type **
% 1- for scRT-qPCR  (CT values)
% 2- for scRT-qPCR  (Expression values - Num. mRNA molecules)
% 3- for scRNA-seq  (Expression values - e.g log(TPM+1) or log(RPKM+1))
% 4- for scRNA-seq  (Expression values - e.g TPM or RPKM)
% 5- for Drop-seq and DropNc-seq (Expression values - e.g UMI or log(UMI+1))
%
% ** INPUTS.format_data **
% 1- Rows= cells and Columns= genes with time/stage info in the last column  
% ------------------------------------------
% Gene1  Gene2  Gene3  ... Genej  Time/Stage
%  27     80     56    ...  69        0
%  73     20     90    ...  45        0
%   .      .      .    ...   .        .
%   .      .      .    ...   .        .
%   .      .      .    ...   .        .
% ------------------------------------------
%
% 2- Rows= genes and Columns= cells with time/stage info in the first row
% ----------------------------------------------
% Time/Stages   0      0      1   ...  16     16
%   Gene1      27     80     56   ...  69      0
%   Gene2      73     20     90   ...  45      0
%   Gene3       .      .      .   ...   .      .
%   Gene4       .      .      .   ...   .      .
%   Genej       .      .      .   ...   .      .
% ----------------------------------------------
% 
% 3- Rows= cells and Columns= genes without time info
% ------------------------------
% Gene1  Gene2  Gene3  ... Genej  
%  27     80     56    ...  69       
%  73     20     90    ...  45       
%   .      .      .    ...   .       
%   .      .      .    ...   .       
%   .      .      .    ...   .       
% ------------------------------
% 
% 4- Rows= genes and Columns= cells without time info
% --------------------------------------------
% Gene1      27     80     56   ...  69      0
% Gene2      73     20     90   ...  45      0
% Gene3       .      .      .   ...   .      .
% Gene4       .      .      .   ...   .      .
% Genej       .      .      .   ...   .      .
% --------------------------------------------
%
% 5- Manual selection from data table
%
% 6- (RECOMMENDED FOR Droplet-based datasets) File .mat containing the following variables:
%    'expMatrix': NxG expression matrix with rows=cells and columns=genes;
%    'cellID': cell array/character vector containing cell IDs (not required)
%    'geneNAMES': cell array containing gene names (not required)
%    'timeline': 1xN numerical vector containing the time information (not required)
%
% 7- (RECOMMENDED FOR Droplet-based datasets) File .mat containing the following variables:
%    'expMatrix': NxG Log-expression matrix with rows=cells and columns=genes;
%    'cellID': cell array/character vector containing cell IDs (not required)
%    'geneNAMES': cell array containing gene names (not required)
%    'timeline': 1xN numerical vector containing the time information (not required)
%
% ** INPUTS.data_selection **
% When cells come from different capture times or stages, users can select 
% cells from specific time or stage for further analysis by CALISTA. 
% For instance, considering a dataset with cells taken from 4 time points, 
% snapshots can be selected, as follows:
% INPUTS.data_selection = [ ] or 1:4 for all time points/cells 
% INPUTS.data_selection = 1          for cells the first time point
% INPUTS.data_selection = [1 3 4]    for cells 1st, 3rd and 4th time points
% By default: INPUTS.data_selection=[]
% 
% ** INPUTS.zeros_genes **
% Gene filtering.
% If INPUTS.zeros_genes is a value in the interval [0,1]=> Remove genes
% with more than certain % of zeros (e.g. INPUTS.zeros_genes=0.9 -> 90%)
% Otherwise Remove genes with more than INPUTS.zeros_genes zeros
% By default:
% if INPUTS.data_type<5 -> INPUTS.zeros_genes=1;
% otherwise -> INPUTS.zeros_genes=3;
%
% ** INPUTS.zeros_cells **
% Cell filtering.
% If INPUTS.zeros_cells is a value in the interval [0,1]=> Remove cells
% with more than certain % of zeros (e.g. INPUTS.zeros_cells=1 -> 100%)
% Otherwise Remove cells with more than INPUTS.zeros_cells zeros
% By default:
% if INPUTS.data_type<5 -> INPUTS.zeros_cells=1;
% otherwise -> INPUTS.zeros_cells=200;
%
% ** INPUTS.cut_variable_genes **
% Most variable gene selection.
% Select X most variable genes with X=min(total number of genes,INPUTS.cut_variable_genes)
% By default: INPUTS.cut_variable_genes = 1000 
%
% ** INPUTS.top_genes **
% Top gene selection based on CALISTA.
% If INPUTS.top_genes is a value in the interval [0,1]=> Select top Y of genes
% with Y=round(INPUTS.top_genes*INPUTS.cut_variable_genes)(e.g. INPUTS.top_genes=0.9 -> 90% of INPUTS.cut_variable_genes)
% Otherwise Select top Y genes with Y=min(INPUTS.top_genes,INPUTS.cut_variable_genes)
% By default: INPUTS.top_genes = min(300, tot number of genes)
% 
% ** INPUTS.cells_2_cut **
% Users can exclude cells from further analysis. 
% 1- Remove cells based on their indices in the expression matrix. 
%    Indices need to be uploaded as a separate csv file. 
% 0- No cells deletion (by default)
% 
% ** INPUTS.use_drop_prob_in_clustering **
% 1- Consider the dropout probability in CALISTA
% 0- Otherwise
% By default:
% If INPUTS.data_type<5 -> INPUTS.use_drop_prob_in_clustering=0;
% Otherwise INPUTS.use_drop_prob_in_clustering=1;
%
% **  INPUTS.cluster_time **
% 1- Run CALISTA clustering for each time point separately
% 0- Run CALISTA clustering on the entire dataset (by default)
%
% ** INPUTS.plot **
% 1- Plot additional figures (by default)
% 0- Otherwise
%
% Output:
% DATA - a structure containing preprocessed expression values for further
% analysis in CALISTA
% 
% INPUTS - a structure containing data description and data preprocessing 
% settings
%
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright.  October 1, 2018.

fprintf('\n\n**** Please upload normalized data. File formats accepted: .txt , .xlxs , .csv, .mat, .loom ****\n\n')

%% CHECK INPUT ARGUMENTS
if nargin <1 %Set to default variables
    error('Not enough input arguments')
end

if ~isfield(INPUTS,'data_type')
     error('Please specify the data type in INPUTS.data_type')
end

if ~isfield(INPUTS,'format_data')
     error('Please specify the format of the data in INPUTS.format_data')
end

if ~isfield(INPUTS,'cluster_time')
     INPUTS.cluster_time=0; 
end

if ~isfield(INPUTS,'cells_2_cut')
     INPUTS.cells_2_cut=0; 
end

if ~isfield(INPUTS,'cut_variable_genes')
     INPUTS.cut_variable_genes=1000;
end

if ~isfield(INPUTS,'data_selection')
     INPUTS.data_selection=[];
end

if INPUTS.cluster_time~=0 && ~isempty(INPUTS.data_selection)
    fprintf('\nCALISTA Time CLustering is active. INPUTS.data_selection set as []. All time points are considered for the analysis\n')
    INPUTS.data_selection=[];
end

if ~isfield(INPUTS,'zeros_genes')
    if INPUTS.data_type<5
        INPUTS.zeros_genes=1;
    else
        INPUTS.zeros_genes=3;
    end
end

if ~isfield(INPUTS,'zeros_cells')
    if INPUTS.data_type<5
        INPUTS.zeros_cells=1;
    else
        INPUTS.zeros_cells=200;
    end
end

if ~isfield(INPUTS,'top_genes')
     INPUTS.top_genes=300;
end

if ~isfield(INPUTS,'plot')
        INPUTS.plot=1; 
     
end

if ~isfield(INPUTS,'use_drop_prob_in_clustering')
    if INPUTS.data_type<5
        INPUTS.use_drop_prob_in_clustering=0;
    else
        INPUTS.use_drop_prob_in_clustering=1;
    end
end

if ~isfield(INPUTS,'remove_batch_effects')
    INPUTS.remove_batch_effects=0;
end

if INPUTS.format_data==8
    if ~isfield(INPUTS.loom_opts,'time_info')
        INPUTS.loom_opts.time_info='Timepoint';
    end
    if ~isfield(INPUTS.loom_opts,'cell_id')
        INPUTS.loom_opts.cell_id='CellID';
    end
    if ~isfield(INPUTS.loom_opts,'gene_name')
        INPUTS.loom_opts.gene_name='Gene';
    end
    
end

%% UPLOADING

[DATA]=new_normalization(INPUTS);

%     cut_variable_genes=1.7;
%     binning_method=2;
%     [DATA]=normalization_large_dataset(data_type,format_data,perczeros_genes,perczeros_cells,cut_variable_genes,cells_2_cut,binning_method);

fprintf('\nDONE.\n');
end