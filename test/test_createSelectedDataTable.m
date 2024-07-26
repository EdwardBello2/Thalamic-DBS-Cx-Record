% try basic inputs



subjID = 'Uva';
% subjID = 'Kramer';

pn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';

% load Nexprocfiles table, which holds spike-sorted spk times
fn = 'NEXprocfiles_';
load([pn, '\', fn, subjID, '.mat']);

% load Sorted_Units table
fn = 'SortedUnits_';
load([pn, '\', fn, subjID, '.mat']);



Selected = createSelectedDataTable(NEXprocfiles, SortedUnits, subjID, 'SUs', pipeParams);
