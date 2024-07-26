function [T] = isNeuPresentForTrialFswp(R_ISI_Fswp)
% go thru each row of original table, and tally which dbs trials each
% unique neuron shows up for, with a "true" value in the table


% annoying step where we must transform the "dbsFrequency" column from
% numerical values to a cell-array of strings. Why the bloody hell didn't I
% keep all labels as strings...
col = find(strcmpi(R_ISI_Fswp.Properties.VariableNames, 'dbsFrequency'));
lHalf = R_ISI_Fswp(:,1:(col - 1));
rHalf = R_ISI_Fswp(:,(col + 1):end);
tabFreqNums = R_ISI_Fswp{:,col};
dbsFrequency = cell(length(tabFreqNums), 1);
for i = 1:length(tabFreqNums)
    dbsFrequency{i,1} = ['hz', num2str(tabFreqNums(i))];
    
end
fTab = [lHalf, table(dbsFrequency), rHalf];


% for Fsweep:
uniqueNeu = unique(fTab.Unit_objectID);
nUniqueNeu = numel(uniqueNeu);
freqLabels = unique(fTab.dbsFrequency)';
nF = numel(freqLabels);
neuTrialPresence = false(nUniqueNeu, nF);
T = array2table(neuTrialPresence);
T.Properties.VariableNames = freqLabels;
T.Properties.RowNames = uniqueNeu;


nRows = size(fTab, 1);

for iRow = 1:nRows
    % for the current row get the neuronID and dbsTrial
    neuID = fTab.Unit_objectID(iRow);
    dbsTr = fTab.dbsFrequency(iRow);
    
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

