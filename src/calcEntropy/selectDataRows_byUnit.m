function Extracted = selectDataRows_byUnit(inputTable, SortedUnits, neuTypeFilter)
% find which rows in inputTable align with the neuron Unit Filter 



N = SortedUnits;
isNeuType = strcmp(N.NeuronType, neuTypeFilter);
matchNeurons = N.objectID(isNeuType);




nRows = size(inputTable,1);
extractRows = false(nRows,1);
for iRow = 1:nRows
    unitID = inputTable.Unit_objectID{iRow};
    
    if any(strcmp(matchNeurons, unitID))
        extractRows(iRow) = true;
    end 
    
end

Extracted = inputTable(extractRows,:);


end
