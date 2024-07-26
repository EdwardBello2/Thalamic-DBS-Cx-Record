function [labels, isLab, isLabMod] = findFsweepModLabels(R_Fswp, alpha)



% Get idx's of modulated (phase-locked) trials
isModFswp = R_Fswp.pVal(:) <= alpha;


% Get positions of each label in Table
labels = unique(R_Fswp.dbsFrequency);
labelTable = R_Fswp.dbsFrequency(:);


% Get total number of trials for each label
nTableRows = size(R_Fswp, 1);
nLabs = numel(labels);

isLab = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLab(:,iLab) = (labelTable == labels(iLab));

end



% Get total number of modulated trials for each label
isLabMod = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLabMod(:,iLab) = (isLab(:,iLab) & isModFswp);
    
end










end