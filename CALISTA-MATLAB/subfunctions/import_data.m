function [DATA]=import_data(INPUTS)
%IMPORT_DATA upload single-cell expression data and perform preprocessing
% CALISTA accepts single-cell RTqPCR and RNA sequencing data. 
% Besides the expression data matrix, users can also provide capture time 
% or cell stage information. 
%
% Usage:
% [DATA]=import_data(INPUTS)
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
% 4- for scRNA-seq  (Expression values - e.g TPM+1 or RPKM)
%
% ** INPUTS.format_data **
% Five data formats are accepted.
%
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
% ** INPUTS.data_selection **
% When cells come from different capture times or stages, users can select 
% cells from specific time or stage for further analysis by CALISTA. 
% For instance, considering a dataset with cells taken from 4 time points, 
% snapshots can be selected, as follows:
% INPUTS.data_selection = [ ] or 1:4 for all time points/cells 
% INPUTS.data_selection = 1          for cells the first time point
% INPUTS.data_selection = [1 3 4]    for cells 1st, 3rd and 4th time points
% 
% ** INPUTS.perczeros_genes **
% Users can exclude genes with a certain percentage of zeros. 
% INPUTS.perczeros_genes = 90 (recommended)
%
% ** INPUTS.perczeros_cells **
% Users can exclude cells with more than a certain percentage of zeros.
% INPUTS.perczeros_cells = 100 (recommended)
% 
% ** INPUTS.cells_2_cut **
% Users can exclude cells from further analysis. 
% 1- Remove cells based on their indices in the expression matrix. 
%    Indices need to be uploaded as a separate csv file. 
% 0- No cells deletion
% 
% ** INPUTS.perc_top_genes **
% Users can specify a certain percentage of genes to be used for CALISTA. 
% For computational efficiency, the number of genes is set to 
% min(200, INPUTS.perc_top_genes * num of cells/100, num of genes)
% INPUTS.perc_top_genes = 100
%
% Output:
% DATA - a structure containing preprocessed expression values for further
% analysis in CALISTA
% 
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. June 1, 2017.

fprintf('\n\n**** Please upload normalized data. File formats accepted: .txt , .xlxs , .csv ****\n\n')

%% CHECK INPUT ARGUMENTS
if nargin <1 %Set to default variables
    error('Not enough input arguments')
end

data_type=INPUTS.data_type;
format_data=INPUTS.format_data;
data_selection=INPUTS.data_selection;
perczeros_genes=INPUTS.perczeros_genes;
perczeros_cells=INPUTS.perczeros_cells;
cut_variable_genes=0;
cells_2_cut=INPUTS.cells_2_cut;

%% UPLOADING
if data_type>=1
    
    [DATA]=normalization(data_type,format_data,perczeros_genes,perczeros_cells,cut_variable_genes,cells_2_cut);
else
    [DATA]=normalization(data_type,format_data,perczeros_genes,perczeros_cells,[],cells_2_cut);
end
%% CHECK INPUT ARGUMENTS

if isempty(data_selection)
    data_selection=(1:size(DATA.singleCELLdata,2)); % Default all the dataset
end

nc=length(data_selection);


%% INITIALIZATION
timeline=[];
genes=DATA.genes;
for i=1:nc
    mRNA=DATA.singleCELLdata{data_selection(i)}';
    nt(i)=size(mRNA,1);
    timeline=[timeline; DATA.timeline(DATA.timeline==DATA.time(data_selection(i)))];
    if i==1
        mRNA_all=mRNA;
    else
        mRNA_all=[mRNA_all;mRNA];
    end
    if length(data_selection)~=DATA.num_time_points
        DATA.singleCELLdata{i}=mRNA';
    end
end

labels=cellstr(num2str((timeline)));
my_stages=str2num(cell2mat(labels));
DATA.totDATA=mRNA_all;

DATA.timeline=timeline;
if length(data_selection)~=DATA.num_time_points
    DATA.time=DATA.time(data_selection);
    DATA.num_time_points=length(DATA.time);
end

DATA.nc=nc;
DATA.nt=nt;
DATA.labels=labels;
DATA.my_stages=my_stages;
DATA.nvars=size(DATA.totDATA,1);


% Select top genes for RNA-seq data
if data_type>=1
    perc_top_genes=INPUTS.perc_top_genes;%mean(DATA.zscores_most_variable_genes)+std(DATA.zscores_most_variable_genes);%INPUTS.perc_top_genes;
    n_top_genes=min([200;round(DATA.nvars*perc_top_genes/100);length(DATA.genes)]);%length(DATA.genes);%min(200,round(DATA.nvars/5)); %round(length(DATA.genes)/2));
    %find(abs(DATA.zscores_most_variable_genes)>=perc_top_genes,1,'last');%
    if DATA.numGENES > n_top_genes
        DATA.genes=DATA.genes(1:n_top_genes);
        DATA.totDATA=DATA.totDATA(:,1:n_top_genes);
        DATA.numGENES=length(DATA.genes);
        for i=1: DATA.num_time_points
            DATA.singleCELLdata{i}=DATA.singleCELLdata{i}(1:n_top_genes,:);
        end
    end
end 

% *** Parameters loading ***
load('Parameters.mat')
DATA.Parameters=Parameters;

end