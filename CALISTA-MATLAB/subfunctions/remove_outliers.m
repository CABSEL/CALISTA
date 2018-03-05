function [totDATA,timeline,outlier_idx,outliers]=remove_outliers(totDATA,timeline)
n_cells=size(totDATA,1);
[cells4bin,bin_center]=hist(totDATA(:));
perc=5;
idx_discard=cells4bin<round(n_cells/100*5);

while(remove==1)
    



end
%     
% function [totDATA,timeline,outlier_idx,outliers]=remove_outliers(totDATA,timeline)
% outlier_idx=[];
% outliers=[];
% for i=1:1%size(totDATA,2)
%     x=totDATA(:,i);
%     LB=mean(x)-2*std(x);
%     UB=mean(x)+2*std(x);
%     temp_idx=find(x>UB | x<LB);
%     outlier_idx = [outlier_idx; temp_idx];
%     outliers=[outliers; totDATA(temp_idx,i)];
%     
%     % %     percntiles = prctile(x,[5 95]);
%     % %     outlierIndex = find(x < percntiles(1) | x > percntiles(2));
%     % %     outlier_idx = [outlier_idx; outlierIndex];
%     % %     outliers=[outliers; totDATA(outlierIndex,i)];
% end
% if length(outlier_idx)~=0
%     outlier_idx=unique(outlier_idx);
%     totDATA(outlier_idx,:)=[];
%     timeline(outlier_idx)=[];
% end
% end

% function [totDATA,timeline,outlier_idx,outliers]=remove_outliers(totDATA,timeline)
% outlier_idx=[];
% outliers=[];
% for i=1:size(totDATA,2)
%     x=totDATA(:,i);
%     outlier_idx = [outlier_idx; find(abs(x - mean(x)) > 4*std(x))];
%     outliers=[outliers; totDATA(outlier_idx,i)];
% end
% outlier_idx=unique(outlier_idx);
% totDATA(outlier_idx,:)=[];
% timeline(outlier_idx)=[];
% end