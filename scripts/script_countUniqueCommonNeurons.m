% Code for counting number of unique neurons that were present for all
% specified conditions

neuron = unique(Tdisp.Unit_objectID);
nUnique = numel(neuron)

condition = {'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};
nConditions = numel(condition);

hasCondition = false(nUnique, nConditions);

% C0



for iNeur = 1:nUnique
    T_iNeur = Tdisp(strcmp(Tdisp.Unit_objectID, neuron{iNeur}),:);

    for iCon = 1:nConditions
        
        if any(strcmp(T_iNeur.dbsContact, condition(iCon)))
            hasCondition(iNeur,iCon) = true;

        end

    end
    
end

hasAll = prod(hasCondition, 2);
nUniqueConsistent = sum(hasAll)

