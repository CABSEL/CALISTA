function Results=plot_cluster_mean_exp(DATA,Results)

fprintf('\nPlotting mean gene expressions...\n')
p=numSubplots(DATA.numGENES);
for i=1:Results.expected_clusters
    D{i}=DATA.totDATA(find(Results.final_groups==i),:)';
    MEAN(:,i)=mean(D{i},2);
    %        MEDIAN(:,i)=median(singleCELLdata{i},2);
%     SD(:,i)=std(D{i},[],2);
end


figure('position', [500, 500, 800, 800]) 
for j=1:DATA.numGENES
    subplot(p(1),p(2),j)
    
    for i=1:Results.expected_clusters
        plot(Results.cluster_progression(i),MEAN(j,i),'.','Color',Results.colorMARK_calista(i,:),'MarkerSize',30)
        hold on
        grid on
        
        %         errorbar(cluster_progression(i),SD(j,i))
    end
    
    pp=plot(Results.TRANSITION.final_graph);
    pp.MarkerSize=15;
    pp.EdgeAlpha=1;
    pp.NodeColor=Results.colorMARK_calista;
    pp.LineWidth=0.5;
    pp.YData=MEAN(j,:);
    pp.XData=Results.cluster_progression;
    pp.NodeLabel=[];
    st=DATA.genes{j};
    title(st)
    xlim([min(Results.stages_name)-0.5 max(Results.stages_name)+0.5])
    
end

Results.singleCELLclusterDATA=D;
pause(3)