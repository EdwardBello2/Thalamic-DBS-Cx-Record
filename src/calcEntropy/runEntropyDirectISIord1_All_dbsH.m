function runEntropyDirectISIord1_All_userVarH(pipeParams)
% Compare Entropies for all data looking at only Entropy during DBS 
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



%% REMOVE individual neurons from analysis...
% ... if they are not present for AT LEAST "nTrialMin" number of trials

nTrialMin = pipeParams.neuronMinimumTrials;

% 1) Create tables to assess individual neuron presence for each DBS trial
% 2) Remove those neurons that don't have at least nTrialMin number of
% trials


% for Csweep:
cT = isNeuPresentForTrialCswp(R_ISI_Cswp);
totTrials = sum(table2array(cT), 2);
hasTrialMin = totTrials >= nTrialMin;
neuIDs = cT.Properties.RowNames;
keepNeurons = neuIDs(hasTrialMin);

% keep only those neurons that have minimum trial-presence
nRows = size(R_ISI_Cswp, 1);
keepRows = false(nRows, 1);
for iRow = 1:nRows
    keepRows(iRow,1) = any(strcmp(R_ISI_Cswp.Unit_objectID{iRow}, keepNeurons));
    
end
R_ISI_Cswp = R_ISI_Cswp(keepRows,:);
R_PSTH_Cswp = R_PSTH_Cswp(keepRows,:);


% for Fsweep:
fT = isNeuPresentForTrialFswp(R_ISI_Fswp);
totTrials = sum(table2array(fT), 2);
hasTrialMin = totTrials >= nTrialMin;
neuIDs = fT.Properties.RowNames;
keepNeurons = neuIDs(hasTrialMin);

% keep only those neurons that have minimum trial-presence
nRows = size(R_ISI_Fswp, 1);
keepRows = false(nRows, 1);
for iRow = 1:nRows
    keepRows(iRow,1) = any(strcmp(R_ISI_Fswp.Unit_objectID{iRow}, keepNeurons));
    
end
R_ISI_Fswp = R_ISI_Fswp(keepRows,:);
R_PSTH_Fswp = R_PSTH_Fswp(keepRows,:);


% %% EXTRACT rows containing all unique single units that have at least ONE pattern-mod trial
% % Note: the two functions below also take unique neuron identities into
% % account, so that all of one cell's rows are included in one table but not
% % the other. 
% 
% alpha = pipeParams.pValAlpha;
% 
% % for Csweep:
% [R_ISI_PhsLckNeu_Cswp, ...
%  R_ISI_NOTPhsLckNeu_Cswp] = separateISItable_PhsLckVsNOTPhsLck(R_ISI_Cswp, R_PSTH_Cswp, alpha);
% 
% 
% % for Fsweep:
% [R_ISI_PhsLckNeu_Fswp, ...
%  R_ISI_NOTPhsLckNeu_Fswp] = separateISItable_PhsLckVsNOTPhsLck(R_ISI_Fswp, R_PSTH_Fswp, alpha);
% 


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

% Display boxplots & ANOVA results of Phase-locked cells 
f1 = figure;
ax1 = axes;
dispDeltaH_Cswp_Boxplot(deltaH_PhsLckNeu_Cswp, R_ISI_PhsLckNeu_Cswp.dbsElectrode, Clabels_PhsLckNeu);
title('ISI-Entropy change in Phase-locked Contact-Sweep trials');
[pC, tabC, statsC] = anova1(deltaH_PhsLckNeu_Cswp, R_ISI_PhsLckNeu_Cswp.dbsElectrode);
[pC, tabC, statsC] = kruskalwallis(deltaH_PhsLckNeu_Cswp, R_ISI_PhsLckNeu_Cswp.dbsElectrode);

% set(f1, 'Position', [168 400 560 420])
set(f1, 'Position', [592 576 560 420])
if isfield(pipeParams, 'CsweepBoxplotXLim')
    ax1.XLim = pipeParams.CsweepBoxplotXLim;
end


% Display boxplots & ANOVA results of NON-Phase-locked cells 
f2 = figure;
ax2 = axes;
dispDeltaH_Cswp_Boxplot(deltaH_NOTPhsLckNeu_Cswp, R_ISI_NOTPhsLckNeu_Cswp.dbsElectrode, Clabels_NOTPhsLckNeu);
title('ISI-Entropy change in NONPhase-locked Contact-Sweep trials');
[pC, tabC, statsC] = anova1(deltaH_NOTPhsLckNeu_Cswp, R_ISI_NOTPhsLckNeu_Cswp.dbsElectrode);
set(f2, 'Position', [753 398 560 420])
if isfield(pipeParams, 'CsweepBoxplotXLim')
    ax2.XLim = pipeParams.CsweepBoxplotXLim;
end


% Display boxplots & ANOVA results of ALL cells
f3 = figure;
ax3 = axes;
dispDeltaH_Cswp_Boxplot(deltaH_AllNeu_Cswp, R_ISI_Cswp.dbsElectrode, Clabels_AllNeu);
title('ISI-Entropy change in ALL Contact-Sweep trials');
[pC, tabC, statsC] = anova1(deltaH_AllNeu_Cswp, R_ISI_Cswp.dbsElectrode);
[pC, tabC, statsC] = kruskalwallis(deltaH_AllNeu_Cswp, R_ISI_Cswp.dbsElectrode);

% set(f3, 'Position', [1336 397 560 420])
set(f3, 'Position', [    15   575   560   420])
if isfield(pipeParams, 'CsweepBoxplotXLim')
    ax3.XLim = pipeParams.CsweepBoxplotXLim;
end

% % rescale phase-locked and All figurse to match
% XLimScaled(1) = min(ax1.XLim(1), ax3.XLim(1));
% XLimScaled(2) = max(ax1.XLim(2), ax3.XLim(2));
% ax1.XLim = XLimScaled;
% ax3.XLim = XLimScaled;

% For Frequency-sweep:

% Display boxplots & ANOVA results of Phase-locked cells 
f4 = figure;
ax4 = axes;
dispDeltaH_Fswp_Boxplot(deltaH_PhsLckNeu_Fswp, convert2stringArray(R_ISI_PhsLckNeu_Fswp.dbsFrequency(:)), ... 
                        Flabels_PhsLckNeu);
title('ISI-Entropy change in Phase-locked Frequency-Sweep trials');
[pF, tabF, statsF] = anova1(deltaH_PhsLckNeu_Fswp, convert2stringArray(R_ISI_PhsLckNeu_Fswp.dbsFrequency(:)));
[pF, tabF, statsF] = kruskalwallis(deltaH_PhsLckNeu_Fswp, convert2stringArray(R_ISI_PhsLckNeu_Fswp.dbsFrequency(:)));

% set(f4, 'Position', [1971 416 560 420])
set(f4, 'Position', [591 59 560 422])

if isfield(pipeParams, 'FsweepBoxplotYLim')
    ax4.YLim = pipeParams.FsweepBoxplotYLim;
end


% Display boxplots & ANOVA results of NON-Phase-locked cells 
f5 = figure;
ax5 = axes;
dispDeltaH_Fswp_Boxplot(deltaH_NOTPhsLckNeu_Fswp, convert2stringArray(R_ISI_NOTPhsLckNeu_Fswp.dbsFrequency(:)), ... 
                        Flabels_NOTPhsLckNeu);
title('ISI-Entropy change in NONPhase-locked Frequency-Sweep trials');
[pF, tabF, statsF] = anova1(deltaH_NOTPhsLckNeu_Fswp, convert2stringArray(R_ISI_NOTPhsLckNeu_Fswp.dbsFrequency(:)));
set(f5, 'Position', [2554 415 560 420])
if isfield(pipeParams, 'FsweepBoxplotYLim')
    ax5.YLim = pipeParams.FsweepBoxplotYLim;
end


% Display boxplots & ANOVA results of ALL cells 
f6 = figure;
ax6 = axes;
dispDeltaH_Fswp_Boxplot(deltaH_AllNeu_Fswp, convert2stringArray(R_ISI_Fswp.dbsFrequency(:)), ...
                        Flabels_AllNeu);
title('ISI-Entropy change in ALL Frequency-Sweep trials');
[pF, tabF, statsF] = anova1(deltaH_AllNeu_Fswp, convert2stringArray(R_ISI_Fswp.dbsFrequency(:)));
[pF, tabF, statsF] = kruskalwallis(deltaH_AllNeu_Fswp, convert2stringArray(R_ISI_Fswp.dbsFrequency(:)));

% set(f6, 'Position', [3137 414 560 420])
set(f6, 'Position', [14 59 560 420])
if isfield(pipeParams, 'FsweepBoxplotYLim')
    ax6.YLim = pipeParams.FsweepBoxplotYLim;
end

% % rescale phase-locked and All figurse to match
% YLimScaled(1) = min(ax4.YLim(1), ax6.YLim(1));
% YLimScaled(2) = max(ax4.YLim(2), ax6.YLim(2));
% ax4.YLim = YLimScaled;
% ax6.YLim = YLimScaled;



%% PRINT info about final dataset results

% Print number of Neurons and number of data points for both "all neurons"
% and just "phase-locked neurons", for Contact-sweep:
c = newline;

disp('--------------- Contact Sweep ---------------');


disp('All Neurons:');

nNeurons = numel(unique(R_ISI_Cswp.Unit_objectID));
nDatapoints = size(R_ISI_Cswp, 1);
disp(['    # Neurons: ', num2str(nNeurons)])
disp(['# data-points: ', num2str(nDatapoints)]);
disp(c)


disp('Phs-Lck Neurons:')

nNeurons = numel(unique(R_ISI_PhsLckNeu_Cswp.Unit_objectID));
nDatapoints = size(R_ISI_PhsLckNeu_Cswp, 1);
disp(['    # Neurons: ', num2str(nNeurons)])
disp(['# data-points: ', num2str(nDatapoints)]);
disp(c)


disp('--------------- Frequency Sweep ---------------');


disp('All Neurons:');

nNeurons = numel(unique(R_ISI_Fswp.Unit_objectID));
nDatapoints = size(R_ISI_Fswp, 1);
disp(['    # Neurons: ', num2str(nNeurons)])
disp(['# data-points: ', num2str(nDatapoints)]);
disp(c)


disp('Phs-Lck Neurons:')

nNeurons = numel(unique(R_ISI_PhsLckNeu_Fswp.Unit_objectID));
nDatapoints = size(R_ISI_PhsLckNeu_Fswp, 1);
disp(['    # Neurons: ', num2str(nNeurons)])
disp(['# data-points: ', num2str(nDatapoints)]);
disp(c)







%
%% Final optional saving of figures

if pipeParams.finalFigs.save == true
    savfPn = pipeParams.finalFigs.savepn;
    
    % save as .eps files
    saveas(f1, [savfPn, '\ISI_deltaH_PhsLck_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f2, [savfPn, '\ISI_deltaH_NotPhsLck_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f3, [savfPn, '\ISI_deltaH_All_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f4, [savfPn, '\ISI_deltaH_PhsLck_Fswp_', pipeParams.subjID ],'epsc');
    saveas(f5, [savfPn, '\ISI_deltaH_NotPhsLck_Fswp_', pipeParams.subjID ],'epsc');
    saveas(f6, [savfPn, '\ISI_deltaH_All_Fswp_', pipeParams.subjID ],'epsc');

    % save as .jpg files
    saveas(f1, [savfPn, '\ISI_deltaH_PhsLck_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f2, [savfPn, '\ISI_deltaH_NotPhsLck_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f3, [savfPn, '\ISI_deltaH_All_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f4, [savfPn, '\ISI_deltaH_PhsLck_Fswp_', pipeParams.subjID ],'jpg');
    saveas(f5, [savfPn, '\ISI_deltaH_NotPhsLck_Fswp_', pipeParams.subjID ],'jpg');
    saveas(f6, [savfPn, '\ISI_deltaH_All_Fswp_', pipeParams.subjID ],'jpg');

end





end % END function

function rowFlabel = convert2stringArray(rowFrequency)

% rowFrequency = R_ISI_PhsLckNeu_Fswp.dbsFrequency(:);
nRows = size(rowFrequency, 1);
rowFlabel = cell(nRows, 1);
for i = 1:nRows, rowFlabel{i,1} = num2str(rowFrequency(i,1)); end


end