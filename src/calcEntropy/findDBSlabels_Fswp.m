function [labels, isLab] = findDBSlabels_Fswp(R_Fswp);



% Get positions of each label in Table
labels = unique(R_Fswp.dbsFrequency);
labelTable = R_Fswp.dbsFrequency(:);


% Get total number of trials for each label
nTableRows = size(R_Fswp, 1);
nLabs = numel(labels);

isLab = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLab(:,iLab) = (labels(iLab) == labelTable);

end



end