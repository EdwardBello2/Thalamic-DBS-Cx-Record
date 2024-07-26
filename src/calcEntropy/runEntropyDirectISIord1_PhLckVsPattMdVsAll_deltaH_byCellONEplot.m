function runEntropyDirectISIord1_PhLckVsPattMdVsAll_deltaH_byCellONEplot(pipeParams)
% Compare Entropies for all data using ISI-estimation and PSTH Entropy 
% to look at deltaH with one plot to compare all vs phase-lock in each sweep Pipeline
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
% Note: the two functions below also take unique neuron identities into
% account, so that all of one cell's rows are included in one table but not
% the other. 

alpha = pipeParams.pValAlpha;

% for Csweep:
[R_ISI_PhsLckNeu_Cswp, R_ISI_NOTPhsLckNeu_Cswp] = separateISItable_PhsLckVsNOTPhsLck(R_ISI_Cswp, R_PSTH_Cswp, alpha);


% for Fsweep:
[R_ISI_PhsLckNeu_Fswp, R_ISI_NOTPhsLckNeu_Fswp] = separateISItable_PhsLckVsNOTPhsLck(R_ISI_Fswp, R_PSTH_Fswp, alpha);

   

%% CALCULATE Delta-Entropies for ISI-based Entropy of each trial

deltaH_PhsLckNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_PhsLckNeu_Cswp);
deltaH_NOTPhsLckNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_NOTPhsLckNeu_Cswp);
deltaH_AllNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_Cswp);

deltaH_PhsLckNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_PhsLckNeu_Fswp);
deltaH_NOTPhsLckNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_NOTPhsLckNeu_Fswp);
deltaH_AllNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_Fswp);



%% GET DBS-trial labels for each row for plotting


[Clabels_PhsLckNeu, isRowClabel_PhsLckNeu] = getAlldbsElectrodeLabels(R_ISI_PhsLckNeu_Cswp);
[Clabels_NOTPhsLckNeu, isRowClabel_NOTPhsLckNeu] = getAlldbsElectrodeLabels(R_ISI_NOTPhsLckNeu_Cswp);
[Clabels_AllNeu, isRowClabel_AllNeu] = getAlldbsElectrodeLabels(R_ISI_Cswp);


[Flabels_PhsLckNeu, isRowFlabel_PhsLckNeu] = getAlldbsFrequencyLabels(R_ISI_PhsLckNeu_Fswp);
[Flabels_NOTPhsLckNeu, isRowFlabel_NOTPhsLckNeu] = getAlldbsFrequencyLabels(R_ISI_NOTPhsLckNeu_Fswp);
[Flabels_AllNeu, isRowFlabel_AllNeu] = getAlldbsFrequencyLabels(R_ISI_Fswp);

% [Flabels, isRowFlabel] = getAlldbsFrequencyLabels(R_Fswp);


%% PLOT boxplots of delta Entropy 

% For Contact-sweep:

% Join the two groups of deltaH's to display
deltaHgroups = [deltaH_PhsLckNeu_Cswp; deltaH_AllNeu_Cswp];

ClabelsGroups = [R_ISI_PhsLckNeu_Cswp.dbsElectrode; R_ISI_Cswp.dbsElectrode];

phslck = cell(numel(deltaH_PhsLckNeu_Cswp), 1);
phslck(:,1) = {'P'};

all = cell(numel(deltaH_AllNeu_Cswp), 1);
all(:,1) = {'A'};

respTypeGroups = [phslck; all];

groups = {ClabelsGroups; respTypeGroups};


% Plot them together in a boxplot
f1 = figure;
dispDeltaH_Cswp_BoxplotGroups(deltaHgroups, groups, Clabels_PhsLckNeu);








% Display boxplots & ANOVA results of Phase-locked cells 
f1 = figure;
dispDeltaH_Cswp_BoxplotGroups(deltaH_PhsLckNeu_Cswp, R_ISI_PhsLckNeu_Cswp.dbsElectrode, Clabels_PhsLckNeu);
title('ISI-Entropy change in Phase-locked Contact-Sweep trials');
[pC, tabC, statsC] = anova1(deltaH_PhsLckNeu_Cswp, R_ISI_PhsLckNeu_Cswp.dbsElectrode);
set(f1, 'Position', [168 400 560 420])

% Display boxplots & ANOVA results of NON-Phase-locked cells 
f2 = figure;
dispDeltaH_Cswp_Boxplot(deltaH_NOTPhsLckNeu_Cswp, R_ISI_NOTPhsLckNeu_Cswp.dbsElectrode, Clabels_NOTPhsLckNeu);
title('ISI-Entropy change in NONPhase-locked Contact-Sweep trials');
[pC, tabC, statsC] = anova1(deltaH_NOTPhsLckNeu_Cswp, R_ISI_NOTPhsLckNeu_Cswp.dbsElectrode);
set(f2, 'Position', [753 398 560 420])


% Display boxplots & ANOVA results of ALL cells
f3 = figure;
dispDeltaH_Cswp_Boxplot(deltaH_AllNeu_Cswp, R_ISI_Cswp.dbsElectrode, Clabels_AllNeu);
title('ISI-Entropy change in ALL Contact-Sweep trials');
[pC, tabC, statsC] = anova1(deltaH_AllNeu_Cswp, R_ISI_Cswp.dbsElectrode);
set(f3, 'Position', [1336 397 560 420])



% For Frequency-sweep:

% Display boxplots & ANOVA results of Phase-locked cells 
f4 = figure;
dispDeltaH_Fswp_Boxplot(deltaH_PhsLckNeu_Fswp, convert2stringArray(R_ISI_PhsLckNeu_Fswp.dbsFrequency(:)), ... 
                        Flabels_PhsLckNeu);
title('ISI-Entropy change in Phase-locked Frequency-Sweep trials');
[pF, tabF, statsF] = anova1(deltaH_PhsLckNeu_Fswp, convert2stringArray(R_ISI_PhsLckNeu_Fswp.dbsFrequency(:)));
set(f4, 'Position', [1971 416 560 420])


% Display boxplots & ANOVA results of Phase-locked cells 
f5 = figure;
dispDeltaH_Fswp_Boxplot(deltaH_NOTPhsLckNeu_Fswp, convert2stringArray(R_ISI_NOTPhsLckNeu_Fswp.dbsFrequency(:)), ... 
                        Flabels_NOTPhsLckNeu);
title('ISI-Entropy change in NONPhase-locked Frequency-Sweep trials');
[pF, tabF, statsF] = anova1(deltaH_NOTPhsLckNeu_Fswp, convert2stringArray(R_ISI_NOTPhsLckNeu_Fswp.dbsFrequency(:)));
set(f5, 'Position', [2554 415 560 420])


% Display boxplots & ANOVA results of Phase-locked cells 
f6 = figure;
dispDeltaH_Fswp_Boxplot(deltaH_AllNeu_Fswp, convert2stringArray(R_ISI_Fswp.dbsFrequency(:)), ...
                        Flabels_AllNeu);
title('ISI-Entropy change in ALL Frequency-Sweep trials');
[pF, tabF, statsF] = anova1(deltaH_AllNeu_Fswp, convert2stringArray(R_ISI_Fswp.dbsFrequency(:)));
set(f6, 'Position', [3137 414 560 420])








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

function rowFlabel = convert2stringArray(rowFrequency)

% rowFrequency = R_ISI_PhsLckNeu_Fswp.dbsFrequency(:);
nRows = size(rowFrequency, 1);
rowFlabel = cell(nRows, 1);
for i = 1:nRows, rowFlabel{i,1} = num2str(rowFrequency(i,1)); end


end