function [Results, DATA, INPUTS]=CALISTA_clustering_classic(DATA,INPUTS)

Parameters=DATA.Parameters;
optimize=INPUTS.optimize;
algorithm=INPUTS.algorithm;
Cluster=INPUTS.cluster;
parallel=INPUTS.parallel;
loops=INPUTS.runs;
max_iter=INPUTS.max_iter;
INPUTS.step=1; % for internal code

%% CLUSTERING

switch optimize
    case 1
        fprintf('\nDetecting the optimal number of clusters with Eigengap heuristics...\n');
        if DATA.nvars>20000
            fprintf('\n***More than 20kcells detected. Please await a few minutes or re-run CALISTA with INPUTS.optimize=2 ***\n');
        end
        
        if INPUTS.data_type==5
            if ~isfield(INPUTS,'max_clusters')
                % *** 2.b-Cell clustering ***
                Results.expected_clusters=50;
            else
                Results.expected_clusters=INPUTS.max_clusters;
            end
            
        else
            if ~isfield(INPUTS,'max_clusters')
                % *** 2.b-Cell clustering ***
                Results.expected_clusters=12;
            else
                Results.expected_clusters=INPUTS.max_clusters;
            end
            
        end
        
        [ my_results] = CALISTA_clustering(DATA,INPUTS,Results,'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops);
        
        % Eigengap heuristic
        consensus=my_results.all.all.consensus;
        [Results.expected_clusters]=eigengap(consensus,Results.expected_clusters);
    case 2
        fprintf('\nDetecting the optimal number of clusters with Akaike information criterion...\n');
        if INPUTS.data_type==5
            if ~isfield(INPUTS,'max_clusters')
                % *** 2.b-Cell clustering ***
                Results.expected_clusters=50;
            else
                Results.expected_clusters=INPUTS.max_clusters;
            end
            
            
        else
            if ~isfield(INPUTS,'max_clusters')
                % *** 2.b-Cell clustering ***
                Results.expected_clusters=12;
            else
                Results.expected_clusters=INPUTS.max_clusters;
            end
            
        end
        
        [Results.expected_clusters] = opt_numb_of_clusts(DATA,INPUTS,Results );
        
        
    case 0
        Results.expected_clusters=input('\nNumber of clusters expected: ');
end

INPUTS.step=2;

tic
% *** 2.b-Cell clustering ***
[ my_results_c] = CALISTA_clustering(DATA,INPUTS,Results,'parallel',parallel,'cluster',Cluster,'algorithm',algorithm,'max_iterations',max_iter,'pop_size',loops);
toc

Results.final_groups=my_results_c.all.all.idx';

% *** 2.c-Relabelling based on time/cell stage info ***
method=3;
[Results,DATA]= find_progression2(Results,DATA,method);

% *** 2.c-Optimal parameter estimation for the final cluster assignment ***
[ my_results_final] = CALISTA_clustering(DATA,INPUTS,Results,'parallel',parallel,'get_k',Results.final_groups,'algorithm',algorithm);

% *** 2.e-Cluster visualization ***
reduction=2;
[Results]=visualization(reduction,DATA,Results);
if isfield(my_results_c.all.all,'consensus')
    my_results_final.all.all.consensus=my_results_c.all.all.consensus;
end
my_results_final.all.all.population=my_results_c.all.all.population;
my_results_final.all.all.runtime_eachloop=my_results_c.all.all.runtime_eachloop';
Results.clustering_struct=my_results_final;

end