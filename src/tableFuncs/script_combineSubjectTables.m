% quick script to combine tables across subjects...


%% Combine Subject NEXprocfiles tables

loadpn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';
savepn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';

tableName = 'NEXprocfiles';

load([loadpn, '\', tableName, '_Kramer']);
NEXprocfiles_Kramer = NEXprocfiles;
nRecs = size(NEXprocfiles_Kramer, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Kramer'};
NEXprocfiles_Kramer = [NEXprocfiles_Kramer, table(subjID)];

load([loadpn, '\', tableName, '_Uva']);
NEXprocfiles_Uva = NEXprocfiles;
nRecs = size(NEXprocfiles_Uva, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Uva'};
NEXprocfiles_Uva = [NEXprocfiles_Uva, table(subjID)];


NEXprocfiles = [NEXprocfiles_Kramer; NEXprocfiles_Uva];

save([savepn, '\', tableName], 'NEXprocfiles');



%% Combine subject SortedUnits tables

loadpn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';
savepn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';

tableName = 'SortedUnits';

load([loadpn, '\', tableName, '_Kramer']);
SortedUnits_Kramer = SortedUnits;
nRecs = size(SortedUnits_Kramer, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Kramer'};
SortedUnits_Kramer = [SortedUnits_Kramer, table(subjID)];

load([loadpn, '\', tableName, '_Uva']);
SortedUnits_Uva = SortedUnits;
nRecs = size(SortedUnits_Uva, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Uva'};
SortedUnits_Uva = [SortedUnits_Uva, table(subjID)];


SortedUnits = [SortedUnits_Kramer; SortedUnits_Uva];

save([savepn, '\', tableName], 'SortedUnits');



%% Combine subject SweepAnalysisTrials4Paper2018 tables

loadpn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';
savepn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';

tableName = 'SweepAnalysisTrials4Paper2018';

load([loadpn, '\', tableName, '_Kramer']);
T_Kramer = SweepAnalysisTrials4Paper2018;
nRecs = size(T_Kramer, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Kramer'};
T_Kramer = [T_Kramer, table(subjID)];

load([loadpn, '\', tableName, '_Uva']);
T_Uva = SweepAnalysisTrials4Paper2018;
nRecs = size(T_Uva, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Uva'};
T_Uva = [T_Uva, table(subjID)];


SweepAnalysisTrials4Paper2018 = [T_Kramer; T_Uva];

save([savepn, '\', tableName], 'SweepAnalysisTrials4Paper2018');



%% Combine subject RateBins tables

loadpn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';
savepn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';

tableName = 'RateBins_60sec_1secBins';

Tstruct = load([loadpn, '\', tableName, '_Kramer']);
Tcell = struct2cell(Tstruct);
T_Kramer = Tcell{1};
nRecs = size(T_Kramer, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Kramer'};
T_Kramer = [T_Kramer, table(subjID)];

Tstruct = load([loadpn, '\', tableName, '_Uva']);
Tcell = struct2cell(Tstruct);
T_Uva = Tcell{1};
nRecs = size(T_Uva, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Uva'};
T_Uva = [T_Uva, table(subjID)];


RateBins = [T_Kramer; T_Uva];

save([savepn, '\', tableName], 'RateBins');



%% Combine subject EntropyPsthLetter_SU_2HzThresh_60pre_60dbs_10000boots tables

loadpn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';
savepn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';

tabN1 = 'EntropyPsthLetter';
tabN2 = '_SU_2HzThresh_60pre_60dbs_10000boots';

Tstruct = load([loadpn, '\', tabN1, '_Kramer', tabN2]);
Tcell = struct2cell(Tstruct);
T_Kramer = Tcell{1};
nRecs = size(T_Kramer, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Kramer'};
T_Kramer = [T_Kramer, table(subjID)];

Tstruct = load([loadpn, '\', tabN1, '_Uva', tabN2]);
Tcell = struct2cell(Tstruct);
T_Uva = Tcell{1};
nRecs = size(T_Uva, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Uva'};
T_Uva = [T_Uva, table(subjID)];


H_letterResults = [T_Kramer; T_Uva];

save([savepn, '\', tabN1, tabN2], 'H_letterResults');



%% Combine subject EntropyPsthLetter_SU_2HzThresh_60pre_60dbs_10000boots tables

loadpn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';
savepn = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';

tabN1 = 'EntropyDirectISI';
tabN2 = '_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots';

Tstruct = load([loadpn, '\', tabN1, '_Kramer', tabN2]);
Tcell = struct2cell(Tstruct);
T_Kramer = Tcell{1};
nRecs = size(T_Kramer, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Kramer'};
T_Kramer = [T_Kramer, table(subjID)];

Tstruct = load([loadpn, '\', tabN1, '_Uva', tabN2]);
Tcell = struct2cell(Tstruct);
T_Uva = Tcell{1};
nRecs = size(T_Uva, 1);
subjID = cell(nRecs, 1);
subjID(:) = {'Uva'};
T_Uva = [T_Uva, table(subjID)];


H_DirectISIResults = [T_Kramer; T_Uva];

save([savepn, '\', tabN1, tabN2], 'H_DirectISIResults');
