function [neurons, isNeuron] = getAllNeurons(Table)
% For a table with "R" rows, Get a 1xN string array of all "N" unique
% neurons in the table, and an RxN boolean matrix of which row R pertains 
% to which label N.


% Get positions of each neuron in Table
neurons = unique(Table.Unit_objectID(:))';
neuronTable = Table.Unit_objectID(:);


% Get total number of trials for each label
nTableRows = size(Table, 1);
nLabs = numel(neurons);

isNeuron = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isNeuron(:,iLab) = strcmp(neurons(iLab), neuronTable);

end



end % END function