function Results=jump_clustering(DATA,cell_assignments)

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
[ my_results_final] = CALISTA_clustering(DATA.totDATA,Parameters.logP,Parameters.sets,Results.expected_clusters,'parallel',parallel,'get_k',Results.final_groups,'algorithm',algorithm);

% *** 2.e-Cluster visualization ***
reduction=2;
[Results]=visualization(reduction,DATA,Results);

Results.clustering_struct=my_results_final;