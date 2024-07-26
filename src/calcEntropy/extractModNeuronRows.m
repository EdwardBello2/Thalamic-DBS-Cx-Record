function Table_ModNeurons = extractModNeuronRows(Table, alpha)
% Outputs a table with rows that pertain to all trials of a particular
% neuron that showed at least one trial as modulated, according to
% Table.pVal < alpha. 

% alpha = pipeParams.pValAlpha;

% Get all unique neurons in Cswp, and idx for which rows pertain to each
[neurons_Cswp, isRowNeuron_Cswp] = getAllNeurons(Table);

% Get all modulated trial idx's in Cswp, and select entire Neurons with at
% least one pattern-modulation
isMod_Cswp = Table.pVal(:) < alpha;
isRowNeuronMod_Cswp = isRowNeuron_Cswp & isMod_Cswp;
isModNeuron_Cswp = any(isRowNeuronMod_Cswp, 1);

% Choose only those neurons that were modulated at least once
modNeurons_Cswp = neurons_Cswp(isModNeuron_Cswp);
isRowModNeuron_Cswp = isRowNeuron_Cswp(:,isModNeuron_Cswp);

% Choose only those table rows pertaining to these "modulated" neurons
isAllModNeurons = any(isRowModNeuron_Cswp, 2);
Table_ModNeurons = Table(isAllModNeurons,:);




end