function Extracted = extractSUs(inputTable, SortedUnits)
% find which rows in table are Single Units (SU)



N = SortedUnits;
isSU = strcmp(N.NeuronType,'SU');
SU_neurons = N.objectID(isSU);




% for each row of the "Analyze" table, check if its "Trial_objectID" field
% is included in the list of trials for Contact Sweep
nRows = size(inputTable,1);
extractRows = false(nRows,1);
for iRow = 1:nRows
    unitID = inputTable.Unit_objectID{iRow};
    
    if any(strcmp(SU_neurons, unitID))
        extractRows(iRow) = true;
    end 
    
end

Extracted = inputTable(extractRows,:);


end
