function [final_results] = CALISTA_clustering(DATA,INPUTS,Results,varargin )


% CALISTA is a novel computational method for Clustering And 
% Lineage Inference through Single Cell Transcriptomics Analysis

if nargin <3
    error('Not enough input variables.')
end

mRNA_all=DATA.totDATA;
k_new=DATA.Parameters.sets;
expected_clusters=Results.expected_clusters;

if length(expected_clusters)>1
    error('Please input only one value for the expected number of clusters')
end

if INPUTS.use_drop_prob_in_clustering
    p_mat_ma=DATA.Parameters.P;
    opt_lambda=DATA.opt_lambda;
    % invert the probability matrix  %NOT LOG AS IN TRANSITION GENES!
    X=p_mat_ma';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define dropout probability
    epsilon=min(min(X(X>0)));
    mRNA_counts=[0 1:200];
    P_temp=exp(-opt_lambda*mRNA_counts)';
    P=diag([1; 1-P_temp(2:end)]);
    P(:,1)=P_temp;
    Z_temp=X*P;
    Z_temp(Z_temp==0)=epsilon;  %%%%%%% SET ZEROS TO A SMALL VALUE
    Z_temp=log(Z_temp)';
else
    Z_temp=DATA.Parameters.logP;
end


%% Default values
% run optimiziation
optimize=true;
% initial population
in_pop=false;
% number of runs
loops=50;
% run in parallel
run_parallel=false;
% cluster algorithm ('greedy_cabsel',greedy_maxdiff', 'cabsel_sabec')
algorithm='greedy_cabsel';
% number of maximum iterations
max_iter=100;
% number of cells and number of genes
[nvars, n_genes]=size(mRNA_all);
nVarargs=length(varargin);
% plot
my_plot=true;
% clustering
my_cluster='kmedoids';

for k=1:2:nVarargs
    
    
    % check if user only wants to get the optimal rate constants for a
    % given assignment
    if strcmp(varargin{k},'get_k')
        % 1 assignement to calculate optimal rate constants
        as_get_k=varargin{k+1};
        % do not optimize
        optimize=false;
        % therefore number of runs =1
        loops=1;
        % if more than one population throw error
        if size(as_get_k,1)>1
            error('Please provide only 1 initial population in order to calculate the rate constants wihout optimization')
        end
        
    end
    
    % check if user has userdefined initial populations
    if strcmp(varargin{k},'set_initial_pop')
        in_pop=true;
        as_in=varargin{k+1};
        
    end
    
    
    % number of runs (default =50)
    if strcmp(varargin{k},'pop_size')
        loops=varargin{k+1};
    end
    
    % number of iterations for each run (default =100)
    if strcmp(varargin{k},'max_iterations')
        max_iter=varargin{k+1};
    end
    
    
    % run function in parallel (default false)
    if strcmp(varargin{k},'parallel')
        run_parallel=varargin{k+1};
    end
    
    % choose algorithm
    if strcmp(varargin{k},'algorithm')
        algorithm=varargin{k+1};
    end
    
    % choose algorithm
    if strcmp(varargin{k},'cluster')
        my_cluster=varargin{k+1};
    end
    
    if strcmp(varargin{k},'use_imputed_data')
        mRNA_all=varargin{k+1};
        [nvars, n_genes]=size(mRNA_all);
        
    end
end

% if no optimization just run for 1 population
if optimize==false
    as_all=as_get_k;
    % if expected_clusters is empty get clusters from the population
    if isempty(expected_clusters)
        expected_clusters=length(unique(as_all));
    end
elseif length(expected_clusters)==1
    % if optimization is true sample random numbers between 1 and expected clusters
    if INPUTS.data_type==5
        seed=round(nvars/2000*expected_clusters);
    else
        seed=round(nvars/2*expected_clusters);
    end
    rng(seed) % For reproducibility
    as_all=randi([1 expected_clusters],loops,nvars);
    
end

% if user provides initial poplation and at the same time wants to
% calculate the optimal constants, throw warning
if in_pop==true
    
    if optimize==false
        warning('No optimization is carried out, since option "get_k" was set')
        
    else
        
        % if length of initial population is smaller than number of runs defined, throw warning
        if size(as_in,1)<loops
            
            loops=size(as_in,1);
            warning('Number of iteratios was set according to size of the initial population')
        end
        
        as_all=as_in;
        
    end
    
end

% define variables used in the function greedy_cabsel
in_population=zeros(loops,nvars);
sum_prob_tot=zeros(1,loops);
opt_idx=struct('run',[]);
opt_idx_a(1:loops)=opt_idx;

% check option parallel was set
if run_parallel==false
    p=0;
else
    
    aaa=gcp('nocreate');
    % if no parpool was started
    if isempty(aaa)
        parpool;
        aaa=gcp('nocreate');
    end
    % number of workes in parfor
    p=aaa.NumWorkers;
    p=p-1;  % number of CPUs used for parallel execution (minus 1 for system stability)
end
%% RUN Greedy Algorithm
display=1;
[my_results ] = greedy_cabsel( as_all,Z_temp,k_new,mRNA_all,n_genes,max_iter,nvars,...
    optimize,opt_idx_a,p,sum_prob_tot,in_population,loops,expected_clusters,algorithm,display);

my_results.cluster=expected_clusters;
all_my_results.all=my_results;

if optimize==true
    
    switch INPUTS.step
        case 1
            population=my_results.population;
            % consensus
            fprintf('\nConstructing consensus matrix...\n');
            consensus=get_consensus(population,p);
        case 2
            
            
            if INPUTS.data_type==5
                fprintf('\nBest cell assignment...\n');
                idx=my_results.best';
            else
                fprintf('\nPerforming clustering based on consensus matrix...\n');
                population=my_results.population;
                % consensus
                fprintf('\nConstructing consensus matrix...\n');
                consensus=get_consensus(population,p);
                data_2_cluster=consensus;
                
                switch my_cluster
                    case 'kmedoids'
                        fprintf('\nKmedois...\n');
                        stream = RandStream('mrg32k3a');  % Random number stream
                        options = statset('UseParallel',run_parallel,'UseSubstreams',run_parallel,...
                            'Streams',stream);
                        idx=kmedoids(data_2_cluster,expected_clusters,'Algorithm','pam','Replicates',5,'Options',options);
                    case 'hierarchical'
                        fprintf('\nHierarchical clustering...\n');
                        Z = linkage(data_2_cluster,'complete');%'average','correlation'); %'complete);
                        idx=cluster(Z,'maxclust',expected_clusters);
                end
            end
   
    
    display=0;
    [my_results_1 ] = greedy_cabsel(idx',Z_temp,k_new,mRNA_all,n_genes,max_iter,nvars,...
        false,opt_idx_a,p,sum_prob_tot,in_population,1,expected_clusters,algorithm,display);
    all_my_results.all.clusterprobabilities=my_results_1.clusterprobabilities;
    all_my_results.all.idx=idx;
    end
    if exist('consensus','var')
    all_my_results.all.consensus=consensus;
    end
end

final_results.all=all_my_results;
end

