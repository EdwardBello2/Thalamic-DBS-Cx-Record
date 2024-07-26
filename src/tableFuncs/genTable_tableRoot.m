% script for creating tableRoot

loadpn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';
savepn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';

% LOAD NEXprocfiles
tableName = 'NEXprocfiles';
load([loadpn, '\', tableName]);
N = NEXprocfiles;

% LOAD SortedUnits
tableName = 'SortedUnits';
load([loadpn, '\', tableName]);
S = SortedUnits;

% LOAD Trial data 
tableName = 'SweepAnalysisTrials4Paper2018';
load([loadpn, '\', tableName]);
T = SweepAnalysisTrials4Paper2018;

% COMBINE the three
R1 = join(N, T, 'LeftKeys', 'Trial_objectID', 'RightKeys', 'objectID', ...
    'KeepOneCopy', 'subjID');

R2 = join(R1, S, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'objectID', ...
    'KeepOneCopy', 'subjID');

% PRUNE the table
R3 = R2;
R3(:,[6, 17]) = [];

% RENAME certain columns
R3.Properties.VariableNames{4} = 'nexFile';
R3.Properties.VariableNames{5} = 'nexFileFolder';
R3.nexFileFolder(:) = {'DataProcessing\NEXprocfiles'};

tableRoot = R3;

save([savepn, '\', 'tableRoot'], 'tableRoot');




