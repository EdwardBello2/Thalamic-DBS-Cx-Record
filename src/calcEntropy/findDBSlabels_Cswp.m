function [labels, isLab] = findDBSlabels_Cswp(R_Cswp);



% Get positions of each label in Table
labels = unique(R_Cswp.dbsElectrode);
labelTable = R_Cswp.dbsElectrode(:);


% Get total number of trials for each label
nTableRows = size(R_Cswp, 1);
nLabs = numel(labels);

isLab = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLab(:,iLab) = strcmp(labels(iLab), labelTable);

end




end