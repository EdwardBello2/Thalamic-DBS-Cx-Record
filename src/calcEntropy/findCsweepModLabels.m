function [labels, isLab, isLabMod] = findCsweepModLabels(R_Cswp, alpha)



% Get idx's of modulated (phase-locked) trials
isModCswp = R_Cswp.pVal(:) <= alpha;


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



% Get total number of modulated trials for each label
isLabMod = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLabMod(:,iLab) = (isLab(:,iLab) & isModCswp);
    
end



end % END function