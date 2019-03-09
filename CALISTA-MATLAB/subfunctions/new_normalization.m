function [DATA]=new_normalization(INPUTS)
% *** Parameters loading ***
load('Parameters.mat')
DATA.Parameters=Parameters;



data_type=INPUTS.data_type;
format_data=INPUTS.format_data;
zeros_genes=INPUTS.zeros_genes;
zeros_cells=INPUTS.zeros_cells;
cells_2_cut=INPUTS.cells_2_cut;
cut_variable_genes=INPUTS.cut_variable_genes;



[FileName,PathName,FilterIndex] = uigetfile('*.*');
filename=strcat(PathName, FileName);


DATA.FileName=filename;
switch format_data
    case 1
        imported_data=importdata(filename);
        NUM=imported_data.data;
        NUM(isnan(NUM(:,1)),:)=[]; % REMOVE ROW WITH AT LEAST ONE NaN
        totDATA=NUM(:,1:end-1);
        timeline=NUM(:,end);
        %         [totDATA,timeline,outlier_idx,outliers]=remove_outliers(totDATA,timeline);
        TXT=imported_data.textdata;
        DATA.genes=TXT(1,1:end-1);
        cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
    case 2
        imported_data=importdata(filename);
        NUM=[imported_data.data(2:end,:)' imported_data.data(1,:)'];
        NUM(isnan(NUM(:,1)),:)=[]; % REMOVE ROW WITH AT LEAST ONE NaN
        totDATA=NUM(:,1:end-1);
        timeline=NUM(:,end);
        TXT=imported_data.textdata';
        DATA.genes=TXT(1,2:end);
        cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
    case 3
        imported_data=importdata(filename);
        NUM=imported_data.data;
        NUM(isnan(NUM(:,1)),:)=[]; % REMOVE ROW WITH AT LEAST ONE NaN
        totDATA=NUM;
        timeline=zeros(size(NUM,1),1);
        TXT=imported_data.textdata;
        if length(TXT(1,:))>size(totDATA,2)
            DATA.genes=TXT(1,2:end);
        else
            DATA.genes=TXT(1,:);
        end
        cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
    case 4
        imported_data=importdata(filename);
        NUM=imported_data.data';
        NUM(isnan(NUM(:,1)),:)=[]; % REMOVE ROW WITH AT LEAST ONE NaN
        totDATA=NUM;
        timeline=zeros(size(NUM,1),1);
        TXT=imported_data.textdata';
        if length(TXT(1,:))>size(totDATA,2)
            DATA.genes=TXT(1,2:end);
        else
            DATA.genes=TXT(1,:);
        end
        cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
    case 5
        imported_data=importdata(filename);
        fprintf(' Text data extracted preview: \n\n')
        disp(imported_data.textdata(1:5,1:7))
        [row,col]=size(imported_data.textdata)
        
        fprintf(' Expression data extracted preview: \n\n')
        disp(imported_data.data(1:5,1:10))
        [row,col]=size(imported_data.data)
        
        rows= input(' * Key starting and ending rows for the expression values (e.g. [1 405]): ');
        cols= input(' * Key starting and ending columns for the expression values (e.g. [2 22525]): ');
        genes_vector=input(' * Press 1 if columns=genes, 0 otherwise: ');
        
        NUM=imported_data.data;
        TXT=imported_data.textdata;
        
        
        %         NUM(isnan(NUM(:,1)),:)=[]; % REMOVE ROW WITH AT LEAST ONE NaN
        totDATA=NUM(rows(1):rows(2),cols(1):cols(2));
        if genes_vector
            cols= input(' * Key starting and ending columns for gene names (e.g. [6 22529]): ');
            DATA.genes=TXT(1,cols(1):cols(2));
        else
            rows= input(' * Key starting and ending rows for gene names (e.g. [6 22529]): ');
            totDATA=totDATA';
            DATA.genes=TXT(rows(1):rows(2),1);
        end
        
        time_info=input(' * Add time info (1-Yes, 0-No): ');
        
        if time_info
            time_vector=input(' * Key the column or row vector (e.g 1) of the EXPRESSION DATA MATRIX with time / cell stage info: ');
            if genes_vector
                timeline=NUM(rows(1):rows(2),time_vector);
            else
                timeline=NUM(time_vector,cols(1):cols(2));
            end
        else
            timeline=zeros(size(totDATA,1),1);
        end
        %         [totDATA,timeline,outlier_idx,outliers]=remove_outliers(totDATA,timeline);
         cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
    
    case 6
        imported_data=importdata(filename);
        if isstruct(imported_data)
            totDATA=imported_data.expMatrix;
        else
            totDATA=imported_data;
        end
        
        if isfield(imported_data,'timeline')
            if size(imported_data.timeline,2)==1
                timeline=imported_data.timeline;
            else
                timeline=imported_data.timeline';
            end
        else
            timeline=zeros(size(totDATA,1),1);
        end
        if isfield(imported_data,'geneNAMES')
            DATA.genes=imported_data.geneNAMES;
        else
            DATA.genes=cellstr(num2str([1:size(totDATA,2)]'));
        end
        
        DATA.imported_genes=DATA.genes;
        % Check if cellID is a cell array otherwise convert into it
        if isfield(imported_data,'cellID')
            if ~iscell(imported_data.cellID)
                cell_ID=cellstr(imported_data.cellID);
            else
                cell_ID=imported_data.cellID;
            end
        else
             cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
        end
    case 7
        imported_data=importdata(filename);
        if isstruct(imported_data)
            totDATA=imported_data.expMatrix;
        else
            totDATA=imported_data;
        end
        
        totDATA=2.^(totDATA)-1;  % transform back first!!
        
        if isfield(imported_data,'timeline')
            if size(imported_data.timeline,2)==1
                timeline=imported_data.timeline;
            else
                timeline=imported_data.timeline';
            end
        else
            timeline=zeros(size(totDATA,1),1);
        end
        if isfield(imported_data,'geneNAMES')
            DATA.genes=imported_data.geneNAMES;
        else
            DATA.genes=num2cell(1:size(totDATA,2))';
        end
        
        DATA.imported_genes=DATA.genes;
        % Check if cellID is a cell array otherwise convert into it
        if isfield(imported_data,'cellID')
            if ~iscell(imported_data.cellID)
                cell_ID=cellstr(imported_data.cellID);
            else
                cell_ID=imported_data.cellID;
            end
        else
            cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
        end
    case 8
        info = h5info(filename);
        totDATA = h5read(filename,'/matrix');
        totDATA=sparse(double(totDATA));
        DATA.genes={};
        timeline=[];
        cell_ID={};
        
        for ii=1:size(info.Groups,1)
            match_col_attrs=ismember(info.Groups(ii).Name,'/col_attrs');
            match_row_attrs=ismember(info.Groups(ii).Name,'/row_attrs');
            if sum(match_col_attrs)==length(match_col_attrs)
                for iii=1:size(info.Groups(ii).Datasets,1)
                     match_cellID=ismember(info.Groups(ii).Datasets(iii).Name,INPUTS.loom_opts.cell_id);
                      if sum(match_cellID)==length(match_cellID)
                          cell_ID=h5read(filename,fullfile('/col_attrs',INPUTS.loom_opts.cell_id));
                      end
                      
                       match_time_info=ismember(info.Groups(ii).Datasets(iii).Name,INPUTS.loom_opts.time_info);
                      if sum(match_time_info)==length(match_time_info)
                          timeline=double(h5read(filename,fullfile('/col_attrs',INPUTS.loom_opts.time_info)));
                          if size(timeline,1)<size(timeline,2)
                              timeline=timeline';
                          end
                      end
                end
                
            end
            if sum(match_row_attrs)==length(match_row_attrs)
                for iii=1:size(info.Groups(ii).Datasets,1)
                     match_gene=ismember(info.Groups(ii).Datasets(iii).Name,INPUTS.loom_opts.gene_name);
                      if sum(match_gene)==length(match_gene)
                          DATA.genes = h5read(filename,fullfile('/row_attrs',INPUTS.loom_opts.gene_name));
                      end
                end
            end
        end
        
        
        imported_data.info=info;
        
        if length(cell_ID)~=size(totDATA,1)
            cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
        end
        if isempty(timeline)
            timeline=zeros(size(totDATA,1),1);
        end
        if isempty(DATA.genes)
            DATA.genes=cellstr(num2str([1:size(totDATA,2)]'));
        end
        
        DATA.imported_genes=DATA.genes;
        % Check if cellID is a cell array otherwise convert into it
        if isempty(cell_ID)
             cell_ID=cellstr(num2str([1:size(totDATA,1)]'));
        end
end


if INPUTS.remove_batch_effects==0
    
    if INPUTS.cluster_time==0
        DATA.imported_data=imported_data;
        [DATA]=new_normalization_subfun(INPUTS,DATA,totDATA,cell_ID,timeline,data_type,zeros_genes,zeros_cells,cells_2_cut,cut_variable_genes,Parameters);
    else
        time=unique(timeline);
        for t=1:length(time)
            fprintf( '\n\n%s %3.2f\n\n', 'IMPORTING DATA AT TIME:  ', time(t))
            temp_totDATA=totDATA(timeline==time(t),:);
            temp_cell_ID=cell_ID(timeline==time(t));
            temp_timeline=timeline(timeline==time(t));
            [TIME_DATA{t}]=new_normalization_subfun(INPUTS,DATA,temp_totDATA,temp_cell_ID,temp_timeline,data_type,zeros_genes,zeros_cells,cells_2_cut,cut_variable_genes,Parameters);
            %         TIME_DATA{t}=DATA;
            if t==1
                TIME_DATA{t}.ORIGINAL_DATA.totDATA=totDATA;
                TIME_DATA{t}.ORIGINAL_DATA.genes=DATA.genes;
                TIME_DATA{t}.ORIGINAL_DATA.cell_ID=cell_ID;
                TIME_DATA{t}.ORIGINAL_DATA.timeline=timeline;
                TIME_DATA{t}.ORIGINAL_DATA.Parameters=Parameters;
            end
        end
        DATA=TIME_DATA;
    end
else
    DATA.imported_data=imported_data;
    [DATA]=new_normalization_subfun_batch_correction(INPUTS,DATA,totDATA,cell_ID,timeline,data_type,zeros_genes,zeros_cells,cells_2_cut,cut_variable_genes,Parameters);
%     fprintf( '\nPlease upload the batch info.\n')
%     [FileName,PathName,FilterIndex] = uigetfile('*.*');
%     filename=strcat(PathName, FileName);
%     batch_num=importdata(filename);
%     
%     batch=unique(batch_num);
%     for t=1:length(batch)
%         fprintf( '\n\n%s %3.2f\n\n', 'IMPORTING DATA AT TIME:  ', batch(t))
%         idx_to_select=find(batch_num==batch(t));
%         temp_totDATA=totDATA(idx_to_select,:);
%         temp_cell_ID=cell_ID(idx_to_select);
%         temp_timeline=timeline(idx_to_select);
%         [BATCH_DATA{t}]=new_normalization_subfun(INPUTS,DATA,temp_totDATA,temp_cell_ID,temp_timeline,data_type,zeros_genes,zeros_cells,cells_2_cut,cut_variable_genes,Parameters);
%         %         TIME_DATA{t}=DATA;
%         if t==1
%             BATCH_DATA{t}.ORIGINAL_DATA.totDATA=totDATA;
%             BATCH_DATA{t}.ORIGINAL_DATA.genes=DATA.genes;
%             BATCH_DATA{t}.ORIGINAL_DATA.cell_ID=cell_ID;
%             BATCH_DATA{t}.ORIGINAL_DATA.timeline=timeline;
%             BATCH_DATA{t}.ORIGINAL_DATA.Parameters=Parameters;
%         end
%     end
    
    % Unify data    
%     if length(BATCH_DATA)==1
%         intersect_GENES=BATCH_DATA{1}.genes;
%     else
%         intersect_GENES=BATCH_DATA{1}.genes;
%         for t=2:length(BATCH_DATA)
%             [intersect_GENES]=intersect(intersect_GENES,BATCH_DATA{t}.genes);
%         end
%     end
%     
%     for t=1:length(BATCH_DATA)
%         [~,idx_intersect_GENES]=ismember(intersect_GENES,BATCH_DATA{t}.genes);
%         BATCH_DATA{t}.genes=intersect_GENES;
%         BATCH_DATA{t}.totDATA=BATCH_DATA{t}.totDATA(:,idx_intersect_GENES);
%         BATCH_DATA{t}.numGENES=length(BATCH_DATA{t}.genes);
%         %     [ my_results_final] = CALISTA_clustering(BATCH_DATA{t},INPUTS,Results_batch{t},'parallel',parallel,'get_k',final_groups{t},'algorithm',algorithm);
%         
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     totDATA=[];
%     totCELL_ID=[];
%     timeline=[];
%     shift=0;
%     Results.final_groups=[];
%     Results.final_groups_batch=[];
%     for t=1:length(BATCH_DATA)
%         %     totGENES=vertcat(totGENES,BATCH_DATA{t}.genes);
%         totCELL_ID=vertcat(totCELL_ID,BATCH_DATA{t}.cell_ID);
%         timeline=[timeline; BATCH_DATA{t}.timeline];
%         temp=final_groups{t}+shift;
%         Results.final_groups=[Results.final_groups temp];
%         Results.final_groups_batch=[Results.final_groups_batch final_groups{t}];
%         shift=max(Results.final_groups);
%         totDATA=[totDATA;  BATCH_DATA{t}.totDATA];
%     end
%     DATA=BATCH_DATA{1}.ORIGINAL_DATA;
%     DATA=BATCH_DATA;
    
end
fprintf('\nFinishing setup of the CALISTA data structure...\n');
end