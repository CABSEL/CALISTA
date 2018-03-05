function [NoAdding]=CheckNumberOfEdges(i,MaxNumberOfEdges)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



if i>MaxNumberOfEdges || i==0
    NoAdding=true;
    
else
    NoAdding=false;
end



end

