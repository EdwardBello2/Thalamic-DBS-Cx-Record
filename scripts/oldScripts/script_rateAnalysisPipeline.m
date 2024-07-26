% script with inputs to Rate Analysis pipelines


%% EntropyDirectISIord1_PhLckVsNonPhLckVsAll_deltaH_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation for 1) Phase-locked cells and 2) general pattern-modulation
% cells

% *************************************************************************

%% Include toolbox: "Thalamic-DBS-Cx-Record" latest version from local git-repo

toolboxPath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record';
addpath(genpath(toolboxPath));



%% Code

clear; 

pipeParams.subjID = 'Uva';
pipeParams.trialType = 'FrequencySweep'; % 'ContactSweep' | 'FrequencySweep'


% If a given intermediate table already exists, overwrite it anyway
% WARNING, may take a long time to redo some tables, hours
pipeParams.overwriteTables = false; % true/false

pipeParams.tablepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 15; % num of bins per decade of log-spaced bins for ISI
pipeParams.neuronMinimumTrials = 4;
pipeParams.intDatapn = ['L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\intermediateData\spkRate\', pipeParams.subjID];
% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = false;


% Specify Range of Entropies that the boxplots will show
% pipeParams.CsweepBoxplotXLim = [-1.5 0.5];
% pipeParams.FsweepBoxplotYLim = [-4.00,-1.00];

sp = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190226_PatternModFigs\FsweepLowerOutliers';
pipeParams.finalFigs.savepn = sp;

rateData = buildAnalysisTable_Rates(pipeParams);

unitID = 'nr17080101-ch01a';
% unitID = 'nr17080101-ch05a';

isUnit = strcmp(unitID, rateData.Unit_objectID);
rateDataSU = rateData(isUnit,:);
exploreNeuronRatesFsweep(rateDataSU);

%% 
% runRateChangeAllTrials(pipeParams);

