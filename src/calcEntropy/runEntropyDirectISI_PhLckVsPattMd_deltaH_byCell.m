function runEntropyDirectISI_PhLckVsPattMd_deltaH_byCell(pipeParams)
% Compare Entropies for all data using ISI-estimation and PSTH Entropy 
% to look at deltaH Pipeline
%
% This pipeline assumes that the following tables exist and are correct:
%
% NEXprocfiles_subjID.mat
% SortedUnits_sibjID.mat
% SweepAnalysisTrials4Paper2018_subjID.mat
%
% where "subjID" is the nhp name (i.e. 'Uva')

% TO-DO




%% DEFAULT PARAMETERS

DEFAULT.subjID        = 'XXX';
DEFAULT.tablepn       = '\';
DEFAULT.neuTypeFilter = 'SU';
DEFAULT.predbsTime    = 60; % seconds
DEFAULT.dbsTime       = 60; % seconds
DEFAULT.hzThresh      = 2; % Hz
DEFAULT.ordH          = 2;
DEFAULT.binsPD        = 20;



%% SET DEFAULT PARAMETERS IF USER HAS NOT SET THEM

% Check for specific User-defined pipeline parameter inputs:
if ~isfield(pipeParams, 'subjID'), pipeParams.subjID               = DEFAULT.subjID; end
if ~isfield(pipeParams, 'tablepn'), pipeParams.tablepn             = DEFAULT.tablepn; end
if ~isfield(pipeParams, 'neuTypeFilter'), pipeParams.neuTypeFilter = DEFAULT.neuTypeFilter; end
if ~isfield(pipeParams, 'predbsTime'), pipeParams.predbsTime       = DEFAULT.predbsTime; end
if ~isfield(pipeParams, 'dbsTime'), pipeParams.dbsTime             = DEFAULT.dbsTime; end
if ~isfield(pipeParams, 'hzThresh'), pipeParams.hzThresh           = DEFAULT.hzThresh; end
if ~isfield(pipeParams, 'ordH'), pipeParams.ordH                   = DEFAULT.ordH ; end
if ~isfield(pipeParams, 'binsPD'), pipeParams.binsPD               = DEFAULT.binsPD; end



%% LOAD NECESSARY METADATA TABLES
% Table names shortened for each of use in code

% load NEXprocfiles_XXX table, where metadata for individual NEX files is
% stored:
load([pipeParams.tablepn, '\', 'NEXprocfiles_', pipeParams.subjID, '.mat']);
NEX = NEXprocfiles; 

% load Sorted_Units_XXX table, where info related to Spike-Sorting is stored:
load([pipeParams.tablepn, '\', 'SortedUnits_', pipeParams.subjID, '.mat']);
Sort = SortedUnits;

% load SweepAnalysisTrials4Paper2018_XXX table, where info on DBS
% parameters is stored
load([pipeParams.tablepn, '\', 'SweepAnalysisTrials4Paper2018_', pipeParams.subjID, '.mat']);
TrialInfo = SweepAnalysisTrials4Paper2018;



%% CREATE/LOAD a table with all data to be analyzed for entropy

nexLabel = 'NEXprocfiles';


% Create a new name for the intermediate table based on the pipeline params
name = buildDataSelectTableName(pipeParams);
IntTableName = [nexLabel, 'Selected', name];

% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createSelectionTableIfNeeded(NEX, Sort, IntTableName, pipeParams);


% LOAD the specified dataSelect Table: NEX_2analyze
load([pipeParams.tablepn, '\', IntTableName, '.mat' ]) 



%% CREATE/LOAD a table with Direct-Entropy estimates for both pre- and DBS periods

% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyDirectISI';
name = buildEntropyDirectISITableName(pipeParams);
H_ISItableName = [baseLabel, name];


% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createEntropyDirectISITableIfNeeded(NEX_2analyze, H_ISItableName, pipeParams);


% LOAD the specified Table: H_DirectISIResults
load([pipeParams.tablepn, '\', H_ISItableName, '.mat']);



%% CREATE TABLE with PSTH-Entropy estimates for both pre- and DBS periods

% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyPsthLetter';
name = buildEntropyPsthLetterTableName(pipeParams);
H_PSTHtableName = [baseLabel, name];


% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createEntropyPsthLetterTableIfNeeded(NEX_2analyze, H_PSTHtableName, pipeParams);


% LOAD the specified table: H_letterResults
load([pipeParams.tablepn, '\', H_PSTHtableName, '.mat' ]);



%% DIVIDE both Entropy tables into Contact-Sweep and Frequency-Sweep
% Note: both tables have rows that correspond exactly to each other

R_ISI = H_DirectISIResults; % simplify variable name
R_PSTH = H_letterResults;


% Divide Analyze table into Contact-sweep and Frequency-sweep
isCswp = strcmp(R_ISI.sweepType(:), 'Contact');
R_ISI_Cswp = R_ISI(isCswp, :);
R_PSTH_Cswp = R_PSTH(isCswp, :);


% Check to make sure only 130Hz is in Contact-Sweep
freqs = R_ISI_Cswp.dbsFrequency;
isRemove = freqs ~= 130;
R_ISI_Cswp(isRemove, :) = [];

freqs = R_PSTH_Cswp.dbsFrequency;
isRemove = freqs ~= 130;
R_PSTH_Cswp(isRemove, :) = [];


% Now Get Frequency Sweeps
isFswp = strcmp(R_ISI.sweepType(:), 'Frequency');
R_ISI_Fswp = R_ISI(isFswp, :);
R_PSTH_Fswp = R_PSTH(isFswp, :);


% Check to make sure only C0 is in Frequency-Sweep
cElec = R_ISI_Fswp.dbsElectrode;
isRemove = ~strcmp(cElec, 'C0');
R_ISI_Fswp(isRemove, :) = [];

cElec = R_PSTH_Fswp.dbsElectrode;
isRemove = ~strcmp(cElec, 'C0');
R_PSTH_Fswp(isRemove, :) = [];



%% EXTRACT rows containing all unique single units that have at least ONE pattern-mod trial

% FOR ANTIDROMIC UPDATE CHANGE THIS PART

alpha = pipeParams.pValAlpha;

% for Csweep:
[R_ISI_PhsLckNeu_Cswp, R_ISI_PatModNeu_Cswp] = separateISItable_PhsLckVsPatMod(R_ISI_Cswp, R_PSTH_Cswp, alpha);


% for Fsweep:
[R_ISI_PhsLckNeu_Fswp, R_ISI_PatModNeu_Fswp] = separateISItable_PhsLckVsPatMod(R_ISI_Fswp, R_PSTH_Fswp, alpha);

   

%% CALCULATE Delta-Entropies for ISI-based Entropy of each trial

deltaH_PhsLckNeu_Cswp = calcDeltaEntropyISI_allTableRows(R_ISI_PhsLckNeu_Cswp);
deltaH_PatModNeu_Cswp = calcDeltaEntropyISI_allTableRows(R_ISI_PatModNeu_Cswp);

deltaH_PhsLckNeu_Fswp = calcDeltaEntropyISI_allTableRows(R_ISI_PhsLckNeu_Fswp);
deltaH_PatModNeu_Fswp = calcDeltaEntropyISI_allTableRows(R_ISI_PatModNeu_Fswp);



%% GET DBS-trial labels for each row for plotting


[Clabels_PhsLckNeu, isRowClabel_PhsLckNeu] = getAlldbsElectrodeLabels(R_ISI_PhsLckNeu_Cswp);
[Clabels_PatModkNeu, isRowClabel_PatModNeu] = getAlldbsElectrodeLabels(R_ISI_PatModNeu_Cswp);


[Flabels_PhsLckNeu, isRowFlabel_PhsLckNeu] = getAlldbsFrequencyLabels(R_ISI_PhsLckNeu_Fswp);
[Flabels_PatModkNeu, isRowFlabel_PatModNeu] = getAlldbsFrequencyLabels(R_ISI_PatModNeu_Fswp);

% [Flabels, isRowFlabel] = getAlldbsFrequencyLabels(R_Fswp);


%% PLOT boxplots of delta Entropy 

% For Contact-sweep:
rowClabel = R_ISI_PhsLckNeu_Cswp.dbsElectrode;
f1 = figure;
dispDeltaH_Cswp_Boxplot(deltaH_PhsLckNeu_Cswp, rowClabel, Clabels_PhsLckNeu);
title('ISI-Entropy change in Phase-locked Contact-Sweep trials');


rowClabel = R_ISI_PatModNeu_Cswp.dbsElectrode;
f2 = figure;
dispDeltaH_Cswp_Boxplot(deltaH_PatModNeu_Cswp, rowClabel, Clabels_PatModkNeu);
title('ISI-Entropy change in Pattern-Mod Contact-Sweep trials');


% For Frequency-sweep:
rowFrequency = R_ISI_PhsLckNeu_Fswp.dbsFrequency(:);
nRows = size(rowFrequency, 1);
rowFlabel = cell(nRows, 1);
for i = 1:nRows, rowFlabel{i,1} = num2str(rowFrequency(i,1)); end

f3 = figure;
dispDeltaH_Fswp_Boxplot(deltaH_PhsLckNeu_Fswp, rowFlabel, Flabels_PhsLckNeu);
title('ISI-Entropy change in Phase-locked Frequency-Sweep trials');



rowFrequency = R_ISI_PatModNeu_Fswp.dbsFrequency(:);
nRows = size(rowFrequency, 1);
rowFlabel = cell(nRows, 1);
for i = 1:nRows, rowFlabel{i,1} = num2str(rowFrequency(i,1)); end

f4 = figure;
dispDeltaH_Fswp_Boxplot(deltaH_PatModNeu_Fswp, rowFlabel, Flabels_PatModkNeu);
title('ISI-Entropy change in Pattern-Mod Frequency-Sweep trials');



%
%% Final optional saving of figures

if pipeParams.finalFigs.save == true
    savfPn = pipeParams.finalFigs.savepn;
    
    % save as .eps files
    saveas(f1, [savfPn, '\ISIH_deltaH_PhsLck_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f2, [savfPn, '\ISIH_deltaH_PatMod_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f3, [savfPn, '\ISIH_deltaH_PhsLck_Fswp_', pipeParams.subjID ],'epsc');
    saveas(f4, [savfPn, '\ISIH_deltaH_PatMod_Fswp_', pipeParams.subjID ],'epsc');
    
    % save as .jpg files
    saveas(f1, [savfPn, '\ISIH_deltaH_PhsLck_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f2, [savfPn, '\ISIH_deltaH_PatMod_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f3, [savfPn, '\ISIH_deltaH_PhsLck_Fswp_', pipeParams.subjID ],'jpg');
    saveas(f4, [savfPn, '\ISIH_deltaH_PatMod_Fswp_', pipeParams.subjID ],'jpg');

end





end % END function