function [labelStrings, isRowLabel] = getAllFrequencyLabels(R_Fswp)
% For a table with "R" rows, Get a 1xL string array of all "L" Frequency 
% labels, and an RxN boolean matrix of which row R pertains to which 
% label L.




% Get positions of each label in Table
labelNums = unique(R_Fswp.dbsFrequency);
labelTable = R_Fswp.dbsFrequency(:);


% Get total number of trials for each label
nTableRows = size(R_Fswp, 1);
nLabs = numel(labelNums);

isRowLabel = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isRowLabel(:,iLab) = (labelTable == labelNums(iLab));

end

labelStrings = cell(1, nLabs);
for iLab = 1:nLabs
    labelStrings{1,iLab} = num2str(labelNums(iLab));
    
end


end % END function