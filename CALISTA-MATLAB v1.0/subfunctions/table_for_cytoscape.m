function []=table_for_cytoscape(connectivityMATRIX,genes,fileNAME,masterREGULATORS)

% gene i=target; gene j=source;
numGENES=length(genes);

if size(genes,2)>size(genes,1)
    genes=genes';
end

%Interactions
interactions=connectivityMATRIX(:);

% %Cutting
% idxZEROS=find(interactions~=0);
% interactions=interactions(idxZEROS);

%Edges
edges={};
for i=1:length(interactions)
    if interactions(i)<0;
        edges=[edges; 'repression'];
    elseif interactions(i)==0;
        edges=[edges; 'nointeraction'];
    else interactions(i)>0;
        edges=[edges; 'activation'];
    end
end


% interactions=num2cell(abs(interactions));
interactions=num2cell(interactions);

%Column of sources
sourceGENES=repmat(genes',1,numGENES)';
sourceGENES=sourceGENES(:);
% sourceGENES=sourceGENES((idxZEROS));
%Column of targets
targetGENES=repmat(genes,1,numGENES)';
targetGENES=targetGENES(:);
% targetGENES=targetGENES((idxZEROS));


if nargin < 4 || isempty(masterREGULATORS)
    excelTABLE=[sourceGENES targetGENES interactions edges];
    titleCOLUMNS={'SourceGENES' 'TargetGENES' 'Interaction' 'Edges'};
    
else
    IN_OUT_degree=repmat(masterREGULATORS',1,numGENES)';
    IN_OUT_degree=IN_OUT_degree(:);
    IN_OUT_degree=IN_OUT_degree((idxZEROS));
    IN_OUT_degree=num2cell(IN_OUT_degree);
    excelTABLE=[sourceGENES targetGENES interactions edges IN_OUT_degree];
    titleCOLUMNS={'SourceGENES' 'TargetGENES' 'Interaction' 'Edges' 'MasterRegulators'};
    
end

table=cell2table(excelTABLE,'VariableNames',titleCOLUMNS);
if exist(fileNAME, 'file')
    delete(fileNAME);
end
writetable(table,fileNAME);
end