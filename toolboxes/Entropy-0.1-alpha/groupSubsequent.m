function [groupRows] = groupSubsequent(inputVector, order)
% This function will find segment a vector into overlapping blocks of
% length n, consisted of adjacent vector elements. Then find the frequency
% of occurrence of each unique block 
% 
% Input: 
% 
% inputVec - input vector (1xn)
% order        - number of adjacent vector elements to group together 
% 
% Output: 
%
% p - frequency of occurrence of each unique block, this is a vector 

% TO-DO:
% - make sure it has ways to handle wierd order inputs

groupRows = zeros(length(inputVector) - order + 1, order);

nCol = order;
for iCol = 1:nCol
    groupRows(:, iCol) = inputVector(iCol:(end - order + iCol));
end 


end 