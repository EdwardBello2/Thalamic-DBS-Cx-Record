function [labels, isRowLabel] = getAlldbsElectrodeLabels(R_Cswp)
% For a table with "R" rows, Get a 1xL string array of all "L" Contact 
% labels, and an RxN boolean matrix of which row R pertains to which 
% label L.

% Get positions of each label in Table
labels = unique(R_Cswp.dbsElectrode)';
labelTable = R_Cswp.dbsElectrode(:);


% Get total number of trials for each label
nTableRows = size(R_Cswp, 1);
nLabs = numel(labels);

isRowLabel = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isRowLabel(:,iLab) = strcmp(labels(iLab), labelTable);

end



end % END function