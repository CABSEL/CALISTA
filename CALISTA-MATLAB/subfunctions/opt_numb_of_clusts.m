function [expected_clusters] = opt_numb_of_clusts(DATA,INPUTS,Results ,varargin)

auto=0;
nVarargs=length(varargin);
for k=1:2:nVarargs    
    if strcmp(varargin{k},'auto')
        auto=varargin{k+1};
        
    end
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
loops=expected_clusters;
% run in parallel
run_parallel=INPUTS.parallel;
% cluster algorithm ('greedy_cabsel',greedy_maxdiff', 'cabsel_sabec')
algorithm='greedy_cabsel';
% number of maximum iterations
max_iter=100;
% number of cells and number of genes
[nvars, n_genes]=size(mRNA_all);
% plot
my_plot=true;
% clustering
my_cluster='kmedoids';
seed=round(nvars/2*expected_clusters);
rng(seed) % For reproducibility
for i=1:expected_clusters
as_all(i,:)=randi([1 i],1,nvars);
end
% define variables used in the function greedy_cabsel
in_population=zeros(loops,nvars);
sum_prob_tot=zeros(1,loops);
opt_idx=struct('run',[]);
opt_idx_a(1:loops)=opt_idx;
display=1;
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

[my_results ] = greedy_cabsel( as_all,Z_temp,k_new,mRNA_all,n_genes,max_iter,nvars,...
    optimize,opt_idx_a,p,sum_prob_tot,in_population,loops,expected_clusters,algorithm,display);
figure()
title('LogP')
xlabel('Number of clusters')
plot(my_results.sum_prob_tot,'*');
figure
title('BIC')
xlabel('Number of clusters')
BIC=-2*my_results.sum_prob_tot+log(nvars)*3*(1:expected_clusters); 
plot(BIC,'*');
hold on
v=-diff(BIC);
% fom=v(1:end-1)./v(2:end);
% [min_AIC,idx_min_AIC]=max(fom);
% % [min_AIC,idx_min_AIC]=min(AIC);
% plot(idx_min_AIC+1,AIC(idx_min_AIC+1),'r*');


Var=BIC(1:end-1)-BIC(2:end); %calculate %variance explained
PC=cumsum(Var)/(BIC(1)-BIC(end));
Cutoff=0.90;
[~,idx_min_BIC]=find(PC>Cutoff,1,'first'); %find the best index
plot(idx_min_BIC+1,BIC(idx_min_BIC+1),'r*');

promt=sprintf('\nOptimal number of cluster according to min BIC %d:  \nIf you want to use this value press enter, \nelse provide desired number of cluster: ',idx_min_BIC+1);

if auto==0
    input1=input(promt);
    
    if isempty(input1)
        expected_clusters=idx_min_BIC+1;
    else
        expected_clusters=input1;
    end
else
    expected_clusters=idx_min_BIC+1;
end
    
