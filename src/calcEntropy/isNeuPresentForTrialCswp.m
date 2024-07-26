function [T] = isNeuPresentForTrialCswp(R_ISI_Cswp)
% go thru each row of original table, and tally which dbs trials each
% unique neuron shows up for, with a "true" value in the table


% for Csweep:
cTab = R_ISI_Cswp;

uniqueNeu = unique(cTab.Unit_objectID);
nUniqueNeu = numel(uniqueNeu);
ContLabels = unique(cTab.dbsElectrode)';
nC = numel(ContLabels);
neuTrialPresence = false(nUniqueNeu, nC);
T = array2table(neuTrialPresence);
T.Properties.VariableNames = ContLabels;
T.Properties.RowNames = uniqueNeu;


nRows = size(cTab, 1);

for iRow = 1:nRows
    % for the current row get the neuronID and dbsTrial
    neuID = cTab.Unit_objectID(iRow);
    dbsTr = cTab.dbsElectrode(iRow);
    
    [~, rowInT] = ismember(neuID, T.Properties.RowNames);
    [~, colInT] = ismember(dbsTr, T.Properties.VariableNames);
    
    T{rowInT,colInT} = true;
    
end



end



% present = table2array(T);
% hasNtrials = sum(present, 2);
% nOrMoreTrials = 7:-1:1;
% 
% for i = 1:nC
%     nNeu(i,1) = sum(hasNtrials >= nC - (i - 1));
%     
% end

