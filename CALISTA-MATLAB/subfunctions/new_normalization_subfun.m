function [DATA]=new_normalization_subfun(INPUTS,DATA,totDATA,cell_ID,timeline,data_type,zeros_genes,zeros_cells,cells_2_cut,cut_variable_genes,Parameters)


% cell_ID=num2cell(1:size(totDATA,1))';

if cells_2_cut==1
    fprintf('\n\n**** Select the csv file containing cell indices to remove ****\n\n')
    [FileName,PathName,FilterIndex] = uigetfile('*.*');
    filename2=strcat(PathName, FileName);
    cells_2_cut = csvread(filename2);
    timeline(cells_2_cut)=[];
    totDATA(cells_2_cut,:)=[];
    cell_ID(cells_2_cut)=[];
    DATA.cut_sort.cell_2_cut=cells_2_cut;
else
    DATA.cut_sort.cell_2_cut=[];
end

abs_indices=1:size(totDATA,1);
%% Cut first genes and cells with high % of zeros

if zeros_genes<=1
    zeros_genes=round(size(totDATA,1)*zeros_genes);
else
    zeros_genes=size(totDATA,1)-zeros_genes;
    
end

if zeros_cells<=1
    zeros_cells=round(size(totDATA,2)*zeros_cells);
else
    zeros_cells=size(totDATA,2)-zeros_cells;
end


if data_type==1
    totDATA(totDATA>28)=28; %ct max
    zeros_genes_in_data=sum(totDATA==max(max(totDATA)));
    zeros_cells_in_data=sum(totDATA==max(max(totDATA)),2);
    
else
    zeros_genes_in_data=sum(totDATA==0);
    zeros_cells_in_data=sum(totDATA==0,2);
    
end


% CUTTING genes
fprintf('\nFiltering genes by number of cells...\n');
idx2cut=find(zeros_genes_in_data>=zeros_genes);
totDATA(:,idx2cut)=[];
DATA.genes(idx2cut)=[];
DATA.cut_sort.idx2cutGENES=idx2cut;
% CUTTING cells
fprintf('\nFiltering cells by number of genes...\n');
idx2cut=find(zeros_cells_in_data>=zeros_cells);
abs_indices(idx2cut)=[];
totDATA(idx2cut,:)=[];
timeline(idx2cut,:)=[];
DATA.cut_sort.idx2cutCELL=idx2cut;
cell_ID(idx2cut,:)=[];

tot_perc_zeros=sum(sum(totDATA==0))*100/(size(totDATA,1)*size(totDATA,2));
fprintf( '\n%s %3.2f\n', '% of zeros in the data ', full(tot_perc_zeros))


% totDATA=totDATA(6574:end,:);

DATA.numGENES_raw_reduced=size(totDATA,2);
DATA.totDATA_raw_reduced=totDATA;
DATA.genes_raw_reduced=DATA.genes;
DATA.cell_raw_reduced=cell_ID;


%% Process Atac-seq data

if data_type==6
    DATA.totDATA_ori_full=totDATA;
    fprintf('\nPerforming NMF...\n')
    eigenpeaks = nnmf(totDATA,INPUTS.top_genes);
    totDATA=eigenpeaks(:,2:end); % Remove first component (technical noise)
    DATA.genes=cellstr(num2str([1:size(totDATA,2)]'));
end

%% Scaling the expression data
if data_type==5 || data_type==6
    fprintf('\nScaling the data...\n');
    DATA.original_totDATA=totDATA;
    DATA.original_genes=DATA.genes;
    % totDATA=log(totDATA+1);
    tot_mRNA_each_cell=sum(totDATA,2);
    new_tot_mRNA_each_cell=median(tot_mRNA_each_cell);
    %     scaling=repmat(new_tot_mRNA_each_cell./tot_mRNA_each_cell,1,size(totDATA,2));
    %     scaled_totDATA=totDATA.*scaling;
    scaling=new_tot_mRNA_each_cell./tot_mRNA_each_cell;
    scaled_totDATA=bsxfun(@times,scaling,totDATA);
    totDATA=scaled_totDATA;
end

% %% Select cell training set
% expressed_genes_per_cell=totDATA~=0;
% idx_training_set=find(sum(expressed_genes_per_cell,2)>900);
% % projection_set=totDATA;
% % projection_set(idx_training_set,:)=[];
% training_set=totDATA(idx_training_set,:);
% temp_totDATA=totDATA;
% totDATA=training_set;


%% Most variable genes

% if data_type>=1
%     fprintf('\nSelecting most variable genes...\n');
%     %Remove all genes that have less than 25 molecules in total over all cells
%     df.top_genes=min(length(DATA.genes),cut_variable_genes);
%     idx_genes2remove=find(sum(totDATA)<=25);
%     totDATA(:,idx_genes2remove)=[];
%     DATA.genes(idx_genes2remove)=[];
%     df.sd=std(totDATA);
%     df.mean=mean(totDATA);%sum(totDATA, 1) ./ sum(totDATA~=0, 1);%mean(nonzeros(totDATA));%nanmean(totDATA);
%     df.cv=df.sd./df.mean;
%     y=log2(df.cv');
%     x=df.mean';
%     f = fit(x,y, 'log2(x.^a+b)','Lower',[-1,0],'Upper',[0,1]);
%     fitted_y=log2(x.^f.a+f.b);
%     df.dispersion=abs(y-fitted_y);
%     df.a=f.a;
%     df.b=f.b;
%     df.fitted_y=fitted_y;
%     [df.dispersion_norm_sorted,df.dispersion_norm_sortedIDX]=sort(df.dispersion,'descend');
%     df.disp_cut_off=df.dispersion_norm_sorted(df.top_genes);
%     df.idx_top_variable_genes=find(df.dispersion>=df.disp_cut_off);
%
%
%     figure;
% %     plot( f, x, y )
%     plot( log2(x), y ,'k*')
%     hold on
%     plot(log2(x(df.idx_top_variable_genes)),y(df.idx_top_variable_genes),'r*')
%     x4fit=linspace(min(x),max(x),1000);
%     plot(log2(x4fit),log2(x4fit.^f.a+f.b),'LineWidth',2)
%     xlabel('log2(mean)')
%     ylabel('log2(CV)')
%     grid on
%     legend('Detected genes','Fitted curve')
%     figure
%     bar(df.dispersion_norm_sorted)
%     grid on
%     title('Distance between empirical and expected log2CV values')
%
% %     totDATA=temp_totDATA;
%     DATA.genes=DATA.genes(df.idx_top_variable_genes);
%     totDATA=totDATA(:,df.idx_top_variable_genes);
%     DATA.df=df;
%
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if data_type==5 
    fprintf('\nSelecting most variable genes...\n');
    df.top_genes=min(length(DATA.genes),cut_variable_genes);
    df.sd=std(totDATA);
    df.mean=mean(totDATA);%sum(totDATA, 1) ./ sum(totDATA~=0, 1);%mean(nonzeros(totDATA));%nanmean(totDATA);
    df.cv=df.sd./df.mean;
    df.var=df.sd.^2;
    df.dispersion=df.var./df.mean;
    df.quantile=[-inf quantile(df.mean,[0.1:0.05:1])];
    scale_factor=1;%1.4826; %from R function
    var_by_bin.mean_bin= discretize(df.mean,df.quantile);
    for i=1:max(unique(var_by_bin.mean_bin))
        actual_idx=find(var_by_bin.mean_bin==i);
        s=df.dispersion(actual_idx);
        var_by_bin.bin_median(i)=median(s);
        var_by_bin.bin_mad(i)=scale_factor*median(abs(s - median(s)));%sum(abs(s-median(s)))/length(s); %median absolute deviation
    end
    df.bin_disp_median=var_by_bin.bin_median(var_by_bin.mean_bin);
    df.bin_disp_mad=var_by_bin.bin_mad(var_by_bin.mean_bin);
    df.dispersion_norm= abs(df.dispersion-df.bin_disp_median)./df.bin_disp_mad;
    [df.dispersion_norm_sorted,df.dispersion_norm_sortedIDX]=sort(df.dispersion_norm,'descend');
    
    
    [nanDISPERSION] = find(isnan(df.dispersion_norm_sorted),1,'last');%
    if isempty(nanDISPERSION)
        nanDISPERSION=0;
    end
    
    df.disp_cut_off_idx=df.dispersion_norm_sortedIDX(1+nanDISPERSION:df.top_genes+nanDISPERSION);
    
    DATA.numGENES_ori_full=size(totDATA,2);
    DATA.totDATA_ori_full=totDATA(:,df.dispersion_norm_sortedIDX);
    DATA.genes_ori_full=DATA.genes(df.dispersion_norm_sortedIDX);
    
    % Sort totDATA based on most variable genes
    totDATA=totDATA(:,df.disp_cut_off_idx);
    DATA.genes=DATA.genes(df.disp_cut_off_idx);
    
    if INPUTS.plot
        % plot dispersion vs. mean for the genes
        figure
        plot(log10(df.mean),log10(df.dispersion),'k*')
        hold on
        grid on
        plot(log10(df.mean(df.disp_cut_off_idx)),log10(df.dispersion(df.disp_cut_off_idx)),'r*')
        xlabel('Log mean expression')
        ylabel('Log dispersion')
        %     legend(' Variable genes')
        %
        %     figure
        %     plot(log(df.mean+1),log(df.cv+1),'r*')
    end
    DATA.df=df;
    DATA.var_by_bin=var_by_bin;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if data_type<5
    fprintf('\nSelecting most variable genes...\n');
    SD=std(totDATA);
    MEAN=sum(totDATA, 1) ./ sum(totDATA~=0, 1);%mean(nonzeros(totDATA));%nanmean(totDATA);
    CV=(SD.^2)./MEAN;
    
    % Remove genes with CV=NaN
    [nanCV] = find(isnan(CV));
    
    CV(nanCV)=[];
    MEAN(nanCV)=[];
    SD(nanCV)=[];
    totDATA(:,nanCV)=[];
    genes=DATA.genes;
    genes(nanCV)=[];
    DATA.genes=genes;
    nbin=20;
    %     if format_data>=10  % equal frequency
    %         chuckCell = {};
    %         vectLength = numel(v);
    %         splitsize = 1/n*vectLength;
    %
    %         for i = 1:n
    %             %newVector(end + 1) =
    %             idxs = [floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1]+1;
    %             chuckCell{end + 1} = v(idxs);
    %         end
    %
    %         [~,~,BIN]=histcounts(MEAN,cat(dim,-Inf, edges, Inf));%,'BinMethod','fd');
    %     else
    [~,~,BIN]=histcounts(MEAN,nbin);%,'BinMethod','fd');
    %     end
    %
    
    z_scoredCV=[];
    most_variable_genes_idx=[];
    zscores_most_variable_genes=[];
    for i=1:nbin
        if length(find(BIN==i))>0
            idx_genes_each_bin{i}=find(BIN==i);
            z_scored_CV{i}=zscore(CV(idx_genes_each_bin{i}));
            %             temp_genes=idx_genes_each_bin{i}(find(abs(z_scored_CV{i})>=cut_variable_genes));
            most_variable_genes_idx=[most_variable_genes_idx idx_genes_each_bin{i}];
            %             temp_zscores=z_scored_CV{i}(find(abs(z_scored_CV{i})>=cut_variable_genes));
            zscores_most_variable_genes=[zscores_most_variable_genes z_scored_CV{i}];
        end
    end
    [aaa,idx] = sort(abs(zscores_most_variable_genes),'descend');
    zscores_most_variable_genes = zscores_most_variable_genes(idx);
    %     cutting_idx=min(length(DATA.genes),cut_variable_genes);%find(abs(zscores_most_variable_genes)>=cut_variable_genes,1,'last');
    most_variable_genes=genes(most_variable_genes_idx);
    most_variable_genes=most_variable_genes(idx);
    %     most_variable_genes=most_variable_genes(1:cutting_idx);
    totDATA=totDATA(:,most_variable_genes_idx);
    totDATA=totDATA(:,idx);
    %     totDATA=totDATA(:,1:cutting_idx);
    DATA.genes=most_variable_genes;
    DATA.zscores_most_variable_genes=zscores_most_variable_genes;%(1:cutting_idx);
    
    DATA.totDATA_ori_full=totDATA;
    DATA.numGENES_ori_full=size(totDATA,2);
    DATA.genes_ori_full=DATA.genes;
end
%% Normalization to max mRNA=200
% Save before normalising

fprintf('\nNormalizing the data...\n');
[DATA, totDATA] = internal_norm_calista(DATA,data_type,totDATA);

%% plots
if INPUTS.use_drop_prob_in_clustering || data_type==5
    fprintf('\nDropout estimation...\n');
    DD=totDATA;
    dropoutRATE=sum(DD==0)/size(DD,1);
    mu=mean(DD);
    
    % Fit dropout probability
    
    f = fit( mu',dropoutRATE', 'exp(-lambda*x)','Lower',0,'Upper',inf);
    opt_lambda=f.lambda;
    DATA.opt_lambda=opt_lambda;
    nnz_totDATA=totDATA~=0;
    nnz_totDATA_cell=sum(nnz_totDATA,2);
    num_molecules_expressed_cell=sum(totDATA,2)./nnz_totDATA_cell;
    
    if INPUTS.plot
        figure;
        subplot(2,2,1)
        % % plot(mu,dropoutRATE,'o','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','c','MarkerFaceColor','c')
        % % %     str= 'Results for cluster ';
        % % %     title( sprintf( '%s %d', str, Results.batches(i) ) );
        plot( f, mu, dropoutRATE )
        xlabel('Mean expression after normalization')
        ylabel('Frequency 0 values')
        grid on
        
        
        
        
        subplot(2,2,2)
        plot(nnz_totDATA_cell,num_molecules_expressed_cell,'o','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','m','MarkerFaceColor','m')
        %     str= 'Results for cluster ';
        %     title( sprintf( '%s %d', str, Results.batches(i) ) );
        xlabel('Genes detected in each cell')
        ylabel('Average molecules detected per cell')
        grid on
        
        subplot(2,2,3)
        hist(nnz_totDATA_cell,100)
        xlabel('Genes detected in each cell')
        ylabel('Num of cells')
        grid on
        
        subplot(2,2,4)
        totDATA_zeros_cells=sum(totDATA==0,2)*100/size(totDATA,2);
        hist(totDATA_zeros_cells,100)
        xlabel('% zero values in each cell')
        ylabel('Num of cells')
        grid on
        pause(1)
    end
    %     DATA.opt_lambda=0.01;
    %     opt_lambda=0.01;
    %
    %% CALISTA top gene selection
    if data_type==5
        fprintf('\nTop genes selection...\n');
        [nvars_temp,numGENES] =size(totDATA);
        p_mat_ma=Parameters.P;
        opt_lambda=DATA.opt_lambda;
        [prob_transition_genes]=get_prob_transition_genes_with_drop_prob( ones(1,nvars_temp),p_mat_ma,totDATA,numGENES,nvars_temp,opt_lambda);
        [sorted_gene_prob,idx_transition_genes]=sort(prob_transition_genes,'ascend');
        transition_gene_ranking=DATA.genes(idx_transition_genes);
        if INPUTS.plot
            figure
            bar(sorted_gene_prob)
            totDATA=totDATA(:,idx_transition_genes);
            DATA.genes=transition_gene_ranking;
            DATA.sorted_gene_prob=sorted_gene_prob;
            title('CALISTA gene ranking')
            xlabel('Num of genes')
            ylabel('\Sigma Log P')
        end
    end
end


% %% PCA top genes
% totDATA=totDATA(:,1:300);
% DATA.genes=DATA.genes(1:300);
% [coeff, score3,latent] = pca(totDATA);
% Results.coeff=coeff;
% explainedVar=latent*100/sum(latent);
% num_pcs2keep=min(50, find(cumsum(explainedVar)>=50,1,'first'));
% [~,idx_top_PC_genes]=sort(abs(coeff(:,1:num_pcs2keep)),'descend');
% idx_top_PC_genes=idx_top_PC_genes(1:3,:);
% unique_idx_top_PC_genes=unique(idx_top_PC_genes(:));
%
% totDATA=totDATA(:,unique_idx_top_PC_genes);
% DATA.genes=DATA.genes(unique_idx_top_PC_genes);
%% Final initialization

DATA.timeline=timeline;
DATA.time=unique(timeline);
DATA.num_time_points=length(DATA.time);
sortTOTdata=[];
sortTIMELINE=[];
idx_sorted_cells=[];
for k=1:DATA.num_time_points
    I=find(timeline==DATA.time(k));
    idx_sorted_cells=[idx_sorted_cells; I];
    cutDIMENSION(k)=length(I);
    sortTOTdata=[sortTOTdata; totDATA(I,:)];
    sortTIMELINE=[sortTIMELINE; timeline(I)];
end
DATA.cutDIMENSION=cutDIMENSION;
DATA.cut_sort.idx_sorted_cells=idx_sorted_cells;
abs_indices=abs_indices(idx_sorted_cells);
DATA.totDATA=sortTOTdata;
DATA.timeline=sortTIMELINE;
DATA.original_data=DATA.original_data(idx_sorted_cells,:);
DATA.totDATA_ori_full=DATA.totDATA_ori_full(idx_sorted_cells,:);
DATA.cell_ID=cell_ID(idx_sorted_cells);
% if data_type>1
%     DATA.max_totDATA=DATA.max_totDATA(idx_sorted_cells,:);
% end
%%%%%%%%
DATA.totDATA(isnan(totDATA)) = 0 ;
%%%%%%%%
[DATA.nvars,DATA.numGENES]=size(DATA.totDATA);

data= mat2cell(DATA.totDATA',DATA.numGENES,cutDIMENSION);% now rows=genes and columns=cells
DATA.singleCELLdata=data;
%% CHECK INPUT ARGUMENTS
data_selection=INPUTS.data_selection;
if isempty(data_selection)
    data_selection=(1:size(DATA.singleCELLdata,2)); % Default all the dataset
end

nc=length(data_selection);


%% INITIALIZATION
timeline=[];
idx_selected_nc=[];
for i=1:nc
    mRNA=DATA.singleCELLdata{data_selection(i)}';
    nt(i)=size(mRNA,1);
    I=find(DATA.timeline==DATA.time(data_selection(i)));
    idx_selected_nc=[idx_selected_nc; I];
    timeline=[timeline; DATA.timeline(I)];
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
DATA.cell_ID=DATA.cell_ID(idx_selected_nc);
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
DATA.cut_sort.abs_indices=abs_indices;

% Select top genes for RNA-seq data
% Select top genes for RNA-seq data
top_genes=INPUTS.top_genes;
if top_genes<=1
    top_genes=round(length(DATA.genes)*top_genes);
end

n_top_genes=min([top_genes;length(DATA.genes)]);

DATA.genes_ori=DATA.genes;
DATA.totDATA_ori=DATA.totDATA;
DATA.numGENES_ori=DATA.numGENES;
if DATA.numGENES > n_top_genes
    DATA.genes=DATA.genes(1:n_top_genes);
    DATA.totDATA=DATA.totDATA(:,1:n_top_genes);
    DATA.numGENES=length(DATA.genes);
    for i=1: DATA.num_time_points
        DATA.singleCELLdata{i}=DATA.singleCELLdata{i}(1:n_top_genes,:);
    end
end

end