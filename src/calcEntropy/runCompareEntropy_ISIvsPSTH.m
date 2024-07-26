function runCompareentropy_ISIvsPSTH(pipeParams)
% Compare Entropies for all data using ISI-estimation and PSTH Entropy Pipeline
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


% Divide Analyze table into Contact-sweep and Current-sweep
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


% Now Frequency Sweep:
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



%% FIND rows that are Modulated, each according to PSTH or ISI based Entropy
% Note: Anything that is phase-locked (PSTH) is pattern modulated, but not
% vice-versa

% ANTIDROMIC UPDATE WILL GO HERE

alpha = pipeParams.pValAlpha;

% For Csweep:
isPhaseLock_Cswp = R_PSTH_Cswp.pVal(:) < alpha;
isBoth_Cswp = (R_ISI_Cswp.pVal(:) < alpha) & (R_PSTH_Cswp.pVal(:) < alpha);
isPattMod_Cswp = (R_ISI_Cswp.pVal(:) < alpha) & (~isBoth_Cswp);


% For Fsweep:
isPhaseLock_Fswp = R_PSTH_Fswp.pVal(:) < alpha;
isBoth_Fswp = (R_ISI_Fswp.pVal(:) < alpha) & (R_PSTH_Fswp.pVal(:) < alpha);
isPattMod_Fswp = (R_ISI_Fswp.pVal(:) < alpha) & (~isBoth_Fswp);



%% GET DBS-trial labels for each row, and modulated rows, for plotting

% For Csweep:
[Clabels, isRowClabel] = getAlldbsElectrodeLabels(R_ISI_Cswp);
isRowClabel_PhaseLock = isRowClabel & isPhaseLock_Cswp;
isRowClabel_PattMod = isRowClabel & isPattMod_Cswp; 


% For Fsweep:
[Flabels, isRowFlabel] = getAlldbsFrequencyLabels(R_ISI_Fswp);
isRowFlabel_PhaseLock = isRowFlabel & isPhaseLock_Fswp;
isRowFlabel_PattMod = isRowFlabel & isPattMod_Fswp;



%% PLOT boxplots of entropy results for both sweeps-types

% Display both sweeps % modulated results:
f1 = figure;
dispModPercent_ISIvsPSTH_Cswp(Clabels, isRowClabel, isRowClabel_PhaseLock, ...
                                                    isRowClabel_PattMod);
title('Observations of Significant Modulation in Contact-Sweep trials');
legend('Phase-Locked', 'No Mod', 'Pattern-Mod');
xlabel('% of all Observations');
ylabel('DBS Contact');


f2 = figure;
dispModPercent_ISIvsPSTH_Fswp(Flabels, isRowFlabel, isRowFlabel_PhaseLock, ...
                                                    isRowFlabel_PattMod);
title('Observations of Significant Modulation in Frequency-Sweep trials');
legend('Phase-Locked', 'No Mod', 'Pattern-Mod');
ylabel('% of all Observations');
xlabel('DBS Frequency (Hz)');


%% Final optional saving of figures

if pipeParams.finalFigs.save == true
    savfPn = pipeParams.finalFigs.savepn;
    
    % save as .eps files
    saveas(f1, [savfPn, '\CompareEntropy_ISIvsPSTH_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f2, [savfPn, '\CompareEntropy_ISIvsPSTH_Fswp_', pipeParams.subjID ],'epsc');
    
    % save as .jpg files
    saveas(f1, [savfPn, '\CompareEntropy_ISIvsPSTH_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f2, [savfPn, '\CompareEntropy_ISIvsPSTH_Fswp_', pipeParams.subjID ],'jpg');

end





end % END function