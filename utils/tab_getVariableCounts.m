% Takes in a table structure "T" and counts the number of "true" values in
% a TF column of data with variable name "varStrCell" (cell with character
% array). This function supports multiple columnts in one call, if you make
% varStrCell a 1 x n cell array of char vectors. 

function [rowCounts] = tab_getVariableCounts(T, varStrCell)
% varStrCell must be a cell vector of chars

nVars = length(varStrCell);
rowCounts = zeros(1, nVars);

for iVar = 1:nVars
    varStr = varStrCell{iVar};
    iVarTF = T{:, varStr};
    rowCounts(iVar) = sum(iVarTF);
    
end

totCount = sum(rowCounts);

if ~(totCount == height(T))
    error('make sure to include only TF values groups that together sum up to total num of table rows');

end

end
