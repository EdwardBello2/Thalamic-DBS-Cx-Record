%% EntropyDirectISIord1_All_userVarH pipeline
% Shows dbs-Entropies according to ISI-based (log-binned) Entropy
% estimation for three options:
% 1) dbsH (Calculated entropies for the during-DBS period of each neuron)
% 2) preH (Calculated entropies for the preDBS period of each neuron)
% 3) deltaH (dbsH - preH)


clear; 

pipeParams.subjID = 'Kramer';

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

% for each neuron observed, exclude it from analysis if it does not stick
% around for at least "neuronMininumTrials" of dbs-trials (i.e. 4 trials)
pipeParams.neuronMinimumTrials = 4; 

pipeParams.entropyType = '%deltaH'; % choose: 'dbsH', 'preH', 'deltaH', '%deltaH'
pipeParams.trialType = 'FrequencySweep'; % choose: 'ContactSweep', 'FrequencySweep'

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.boxplotPres = 'cleanWoutliers';

pipeParams.finalFigs.save = false;

% Specify Range of Entropies that the boxplots will show
% pipeParams.CsweepBoxplotXLim = [-1.5 0.5];
% pipeParams.FsweepBoxplotYLim = [-4.00,-1.00];

sp = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190226_PatternModFigs\FsweepLowerOutliers';
pipeParams.finalFigs.savepn = sp;
% [Entropy_AllNeu, Entropy_grpLabel] = runEntropyDirectISIord1_All_userVarH(pipeParams);
[Entropy_AllNeu, Entropy_grpLabel] = runEntropyDirectISIord1_All_userVarH_percPlot(pipeParams);

% runPipeline('EntropyDirectISIord1_All_userVarH', pipeParams);




