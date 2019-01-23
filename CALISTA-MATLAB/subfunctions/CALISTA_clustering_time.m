function [Results, DATA, INPUTS]=CALISTA_clustering_time(DATA,INPUTS)

Parameters=DATA{1}.Parameters;
optimize=INPUTS.optimize;
algorithm=INPUTS.algorithm;
Cluster=INPUTS.cluster;
parallel=INPUTS.parallel;
loops=INPUTS.runs;
max_iter=INPUTS.max_iter;
if ~isfield(INPUTS,'max_clusters')
    % *** 2.b-Cell clustering ***
    Results.expected_clusters=50;
else
    Results.expected_clusters=INPUTS.max_clusters;
end

TIME_DATA=DATA;
num_time_points=length(DATA);
min_clust=1;
for t=num_time_points:-1:1
    INPUTS.step=1; % for internal code
    DATA=TIME_DATA{t};
    fprintf( '\n\n%s %3.2f\n\n', 'CLUSTERING DATA AT TIME:  ', DATA.time)
    if Results.expected_clusters==1
        optimize=3;
    end
    switch optimize
        case 1
            fprintf('\nDetecting the optimal number of clusters with Eigengap heuristics...\n');
            if DATA.nvars>20000
                fprintf('\n***More than 20kcells detected. Please await a few minutes or re-run CALISTA with INPUTS.optimize=2 ***\n');
            end
            [ my_results] = CALISTA_clustering(DATA,INPUTS,Results,'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops);
            % Eigengap heuristic
            consensus=my_results.all.all.consensus;
            if t==num_time_points
                [Results.expected_clusters]=eigengap(consensus,Results.expected_clusters,'min_clust',min_clust);
            else
                [Results.expected_clusters]=eigengap(consensus,Results.expected_clusters,'min_clust',min_clust,'auto',1); %automatically evaluate opt number of clusters
            end
        case 2
            fprintf('\nDetecting the optimal number of clusters with BIC criterion...\n');
            if t==num_time_points
                [Results.expected_clusters] = opt_numb_of_clusts(DATA,INPUTS,Results );
            else
                [Results.expected_clusters] = opt_numb_of_clusts(DATA,INPUTS,Results,'auto',1 );%automatically evaluate opt number of clusters
            end
        case 0
            expected_clusters=input('\nNumber of clusters expected: ');
            if t==num_time_points
                Results.expected_clusters=expected_clusters;
            else
                
                while(Results.expected_clusters<expected_clusters)
                    fprintf('\nNumber of clusters selceted exceeds the previous number of cluster predicted. Please try again.\n');
                    expected_clusters=input('\nNumber of clusters expected: ');
                end
                Results.expected_clusters=expected_clusters;
            end
        case 3
           Results.expected_clusters=1;  
    end
    if t<=3
        min_clust=1;
    else
        min_clust=ceil(Results.expected_clusters/3);
    end
    INPUTS.step=2;    
    tic
    % *** 2.b-Cell clustering ***
    [ my_results_c] = CALISTA_clustering(DATA,INPUTS,Results,'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops);    
    runtime(t,:)=my_results_c.all.all.runtime_eachloop;
    final_groups{t}=my_results_c.all.all.idx';
    toc
    fprintf('\n******************************************************\n');
end

totGENES=[];
totCELL_ID=[];
timeline=[];
shift=0;
Results.final_groups=[];
for t=1:length(TIME_DATA)
    totGENES=vertcat(totGENES,TIME_DATA{t}.genes);
    totCELL_ID=vertcat(totCELL_ID,TIME_DATA{t}.cell_ID);
    timeline=[timeline; TIME_DATA{t}.timeline];
    temp=final_groups{t}+shift;
    Results.final_groups=[Results.final_groups temp];
    shift=max(Results.final_groups);
end
totGENES=unique(totGENES);
Results.expected_clusters=max(Results.final_groups);
check=length(totCELL_ID)==length(unique(totCELL_ID));
if ~check
    error('*** Cell IDs must be unique ***')
end

%% Now unify the data
DATA=[];
data_type=INPUTS.data_type;
fprintf('\nUnifying the data...\n');
totDATA=TIME_DATA{1}.ORIGINAL_DATA.totDATA;
genes=TIME_DATA{1}.ORIGINAL_DATA.genes;
cell_ID=TIME_DATA{1}.ORIGINAL_DATA.cell_ID;

TIME_DATA{1}.ORIGINAL_DATA.numGENES_raw_reduced=size(totDATA,2);
TIME_DATA{1}.ORIGINAL_DATA.totDATA_raw_reduced=totDATA;
TIME_DATA{1}.ORIGINAL_DATA.genes_raw_reduced=genes;
TIME_DATA{1}.ORIGINAL_DATA.cell_ID_raw_reduced=cell_ID;



[~,idx_genes]=ismember(totGENES,genes);
[~,idx_cell_ID]=ismember(totCELL_ID,cell_ID);
totDATA=totDATA(idx_cell_ID,idx_genes);
DATA.genes=totGENES;
DATA.cut_sort.idx2cutGENES=[];
DATA.cut_sort.idx2cutCELL=[];
tot_perc_zeros=sum(sum(totDATA==0))*100/(size(totDATA,1)*size(totDATA,2));
fprintf( '\n%s %3.2f\n', '% of zeros in the data ', full(tot_perc_zeros))
% Scaling the expression data
if data_type==5
    fprintf('\nScaling the data...\n');
    DATA.original_totDATA=totDATA;
    DATA.original_genes=DATA.genes;
    % totDATA=log(totDATA+1);
    tot_mRNA_each_cell=sum(totDATA,2);
    new_tot_mRNA_each_cell=median(tot_mRNA_each_cell);
    scaling=repmat(new_tot_mRNA_each_cell./tot_mRNA_each_cell,1,size(totDATA,2));
    scaled_totDATA=totDATA.*scaling;
    totDATA=scaled_totDATA;
end

DATA.original_data=totDATA;
DATA.min_totDATA=min(min(totDATA));
DATA.log_max_mRNA=log2(200);

fprintf('\nNormalizing the data...\n');
[DATA, totDATA] = internal_norm_calista(DATA,data_type,totDATA);

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
    
end



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
DATA.totDATA=sortTOTdata;
DATA.timeline=sortTIMELINE;
DATA.original_data=DATA.original_data(idx_sorted_cells,:);
DATA.cell_ID=cell_ID(idx_sorted_cells);
DATA.totDATA(isnan(totDATA)) = 0 ;
%%%%%%%%
[DATA.nvars,DATA.numGENES]=size(DATA.totDATA);

data= mat2cell(DATA.totDATA',DATA.numGENES,cutDIMENSION);% now rows=genes and columns=cells
DATA.singleCELLdata=data;
if INPUTS.cluster_time==0
    DATA.imported_data=imported_data;
end
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



% Select top genes for RNA-seq data
top_genes=length(DATA.genes);
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

DATA.TIME_DATA=TIME_DATA;
DATA.Parameters=TIME_DATA{1}.ORIGINAL_DATA.Parameters;
%% CLUSTERING

% *** 2.c-Relabelling based on time/cell stage info ***
method=3;
[Results,DATA]= find_progression2(Results,DATA,method);

% *** 2.c-Optimal parameter estimation for the final cluster assignment ***
[ my_results_final] = CALISTA_clustering(DATA,INPUTS,Results,'parallel',parallel,'get_k',Results.final_groups,'algorithm',algorithm);

% *** 2.e-Cluster visualization ***
% reduction=2;
[Results]=visualization(INPUTS,DATA,Results);
if isfield(my_results_c.all.all,'consensus')
    my_results_final.all.all.consensus=my_results_c.all.all.consensus;
end
my_results_final.all.all.population=my_results_c.all.all.population;
my_results_final.all.all.runtime_eachloop=my_results_c.all.all.runtime_eachloop';
Results.clustering_struct=my_results_final;
Results.runtime=runtime;
end