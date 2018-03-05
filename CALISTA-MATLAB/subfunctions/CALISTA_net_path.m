function Results=CALISTA_net_path(Results,INPUTS)

% CALISTA_net_path_main function infers the co-expression network for each
% developmental path detected previously
% 
% Usage:
% 
% Results=CALISTA_net_path_main(Results,INPUTS)
% 
% Please refer to the README file for further details
% 
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. June 1, 2017.

if nargin <2 
    error('Not enough input arguments')
end

method=INPUTS.method;
switch method
    case 1
        value_cutoff=0.4;
        pvalue_cutoff=0.05;
    case 2
        value_cutoff=0.8;
        pvalue_cutoff=0.01;
end

fprintf('\nCALISTA_net_path is running...\n')
switch method
    case 1
        for k=1:length(Results.PATH.CELL_path)
            pop=Results.PATH.smoothExpr{k};     
            [bbb,pvalue]=partialcorr(pop);
            Theta{k}=bbb.*((pvalue<=pvalue_cutoff)& (abs(bbb)>=value_cutoff));
            Theta{k}(isnan(Theta{k}))=0;
        end
        str='Partial Correlation';
    case 2
        for k=1:length(Results.PATH.CELL_path)
            pop=Results.PATH.smoothExpr{k};            
            [bbb,pvalue]=corr(pop);
            Theta{k}=bbb.*((pvalue<=pvalue_cutoff)& (abs(bbb)>=value_cutoff));
            Theta{k}(isnan(Theta{k}))=0;
        end
        str='Correlation';
end

figure
[p,~]=numSubplots(length(Results.PATH.CELL_path));
for i=1:length(Results.PATH.CELL_path)
    subplot(p(1),p(2),i)
    imagesc(Theta{i})
    colormap('jet');   % set colormap
    caxis([-1 1])
    colorbar;          % show color scale
    %Ylabel
    yvec=1:length(Results.GENES.tot_transition_genes);
    %Xlabel
    xvec=yvec;
    plot_settings(xvec,yvec,Results.GENES.tot_transition_genes,Results.GENES.tot_transition_genes)
    xlabel('Target Gene j','FontWeight','Bold')
    ylabel('Source Gene i','FontWeight','Bold')
    title([' Path num ' num2str(i)])
%     suptitle(str)
    
end

Results.PATH.Theta=Theta;

pause(3)