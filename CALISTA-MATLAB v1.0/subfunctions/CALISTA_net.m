function Results=CALISTA_net(Results,method)

fprintf('\nCALISTA_net is running...\n')

switch method
   
    case 1
        for k=1:length(Results.TRANSITION.nodes_connection)
            pop=Results.ORDERING.transition_genes_Expr{k};
%             pop=Results.ORDERING.edge_mRNA_all{k};
            [bbb,pvalue]=partialcorr(pop);
            Theta{k}=bbb.*((pvalue<=0.05)& (abs(bbb)>=0.4));
            Theta{k}(isnan(Theta{k}))=0;
        end
    case 2
        for k=1:length(Results.TRANSITION.nodes_connection)
            pop=Results.ORDERING.transition_genes_Expr{k};
%             pop=Results.ORDERING.edge_mRNA_all{k};
            [bbb,pvalue]=corr(pop);
            Theta{k}=bbb.*((pvalue<=0.01)& (abs(bbb)>=0.9));
            Theta{k}(isnan(Theta{k}))=0;
        end
%     case 3
%         parfor k=1:length(Results.TRANSITION.nodes_connection)
%             pop=Results.ORDERING.transition_genes_Expr{k};
%             %              pop=Results.ORDERING.edge_mRNA_all{k};
%             S = cov(pop);
%             tmp = max(max(abs(S)));
%             lambdaList = (tmp ./ 10):(tmp ./ 10):tmp;
%             pathLength = length(lambdaList);
%             rho=lambdaList(1);
%             [Theta{k} W{k}] = graphicalLasso(S, rho);
%             Theta{k}=abs(Theta{k})>0;
%         end
end

figure
[p,~]=numSubplots(length(Results.TRANSITION.nodes_connection));
for i=1:length(Results.TRANSITION.nodes_connection)
    subplot(p(1),p(2),i)
    imagesc(Theta{i})
    colormap('jet');   % set colormap
    caxis([-1 1])
    colorbar;          % show color scale
    %Ylabel
    yvec=1:length(Results.GENES.actual_transition_genes{i});
    %Xlabel
    xvec=yvec;
    plot_settings(xvec,yvec,Results.GENES.actual_transition_genes{i},Results.GENES.actual_transition_genes{i})
    xlabel('Target Gene j','FontWeight','Bold')
    ylabel('Source Gene i','FontWeight','Bold')
     title([' edge ' num2str(i)])
    
end

Results.NET.Theta=Theta;

pause(3)
