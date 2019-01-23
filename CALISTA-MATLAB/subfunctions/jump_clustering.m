function Results=jump_clustering(DATA,INPUTS,cell_assignments)

if ~isfield(INPUTS,'use_drop_prob_in_clustering')
    if INPUTS.data_type<5
        INPUTS.use_drop_prob_in_clustering=0;
    else
        INPUTS.use_drop_prob_in_clustering=1;
    end
end

if ~isfield(INPUTS,'optimize')
    if INPUTS.time==0
        INPUTS.optimize=0;
    else
        INPUTS.optimize=1;
    end
end

if ~isfield(INPUTS,'parallel')
     INPUTS.parallel=1; 
end

if ~isfield(INPUTS,'runs')
    if INPUTS.data_type<5
        INPUTS.runs=50; 
    else
        INPUTS.runs=20; 
    end
     
end

if ~isfield(INPUTS,'max_iter')
     INPUTS.max_iter=200;
end

if ~isfield(INPUTS,'algorithm')
     INPUTS.algorithm='greedy_cabsel';
end

if ~isfield(INPUTS,'cluster')
     INPUTS.cluster='kmedoids';
end

if ~isfield(INPUTS,'plot_tsne')
    INPUTS.plot_tsne=1;
end

if  INPUTS.plot_tsne==1
    if ~isfield(INPUTS,'opts')
    INPUTS.tsne_opts.perplexity = 30;
    INPUTS.tsne_opts.input_data = 0;
    else
        if ~isfield(INPUTS.tsne_opts,'perplexity')
           INPUTS.tsne_opts.perplexity=30;
        end
         if ~isfield(INPUTS.tsne_opts,'input_data')
           INPUTS.tsne_opts.input_data=0;
         end
    end
end



Parameters=DATA.Parameters;
algorithm='greedy_cabsel';  
parallel=0;


if size(cell_assignments,1)>size(cell_assignments,2)
    cell_assignments=cell_assignments';
end

% Create new Results structure
Results.cluster_predicted=unique(cell_assignments);
Results.expected_clusters=length(Results.cluster_predicted);
Results.final_groups=cell_assignments;

% *** 2.c-Relabelling based on time/cell stage info ***
method=3;
[Results,DATA]= find_progression2(Results,DATA,method);


% *** 2.c-Optimal parameter estimation for the final cluster assignment ***
[ my_results_final] = CALISTA_clustering(DATA,INPUTS,Results,'parallel',parallel,'get_k',Results.final_groups,'algorithm',algorithm);

% *** 2.e-Cluster visualization ***
% reduction=2;
[Results]=visualization(INPUTS,DATA,Results);

Results.clustering_struct=my_results_final;