% try basic inputs
%
% ***NOTE***
% EntropyDirectISIord1_All_userVarH pipeline is my preferred final pipeline


%% EntropyDirect ISI pipeline

clear; 

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 20; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = true;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190122_Group figures';

runPipeline('EntropyDirectISI', pipeParams);



%% Entropy PSTH-letter pipeline

clear;

pipeParams.tablepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anywaypipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz

% PSTH-bin parameters
pipeParams.binEdgeLeft = 0; % seconds
pipeParams.binWidth = 0.5/1000; % seconds
pipeParams.binNum = 15;
pipeParams.trimPSTH = true; % if true, the first bin in the PSTH is removed; if false, no change

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% Final figures to be saved somwhere?
pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190122_Group figures';

runPipeline('EntropyPSTHletter', pipeParams);



%% Entropy PSTH-letter by-PSTHclass pipeline

clear;

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anywaypipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz

% PSTH-bin parameters
pipeParams.binEdgeLeft = 0; % seconds
pipeParams.binWidth = 0.5/1000; % seconds
pipeParams.binNum = 15;
pipeParams.trimPSTH = true; % if true, the first bin in the PSTH is removed; if false, no change


% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% Final figures to be saved somwhere as both eps and jpg
pipeParams.finalFigs.save = true;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190122_Group figures';


runPipeline('EntropyPSTHletter_byPSTHclass', pipeParams);



%% EntropyDirect ISI by PSTHclass pipeline

clear; 

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anywaypipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz


% Direct-Entropy ISI log-bin parameters
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 20; % num of bins per decade of log-spaced bins for ISI


% Phase-lock entropy PSTH-bin parameters
pipeParams.binEdgeLeft = 0; % seconds
pipeParams.binWidth = 0.5/1000; % seconds
pipeParams.binNum = 15;


% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;


% save figures
pipeParams.maxdH = -1.5;
pipeParams.dispPsthType = 'other';
pipeParams.finalFigs.save = true;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190122_Group figures';

runPipeline('EntropyDirectISI_byPSTHclass', pipeParams);



%% Stim-Evoked Spike Latency analysis pipeline

clear; 

pipeParams.subjID = 'Kramer';
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anywaypipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz


% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;


% Phase-lock entropy PSTH-bin parameters
pipeParams.binEdgeLeft = 0; % seconds
pipeParams.binWidth = 0.5/1000; % seconds
pipeParams.binNum = 15;
pipeParams.trimPSTH = true; % if true, the first bin in the PSTH is removed; if false, no change

% Latency Win defines the window of time after each pulse that we look for
% spike latencies
pipeParams.latencyWin = [0.5/1000, 2/1000]; %seconds 



% save figures
pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190122_Group figures';


runPipeline('EvokedSpikeLatency', pipeParams);




%% Compare ISI Entropies to PSTH Entropies pipeline

clear; 

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz

% ISI-Entropy parameters
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 20; % num of bins per decade of log-spaced bins for ISI

% PSTH-bin parameters
pipeParams.binEdgeLeft = 0; % seconds
pipeParams.binWidth = 0.5/1000; % seconds
pipeParams.binNum = 15;
pipeParams.trimPSTH = true; % if true, the first/last bins in the PSTH are removed; if false, no change

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% % save figures
% pipeParams.maxdH = -1.5;

% Final figures to be saved
pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

runPipeline('CompareEntropy_ISIvsPSTH', pipeParams);



%% Delta-Entropy by Cell PSTH-letter pipeline
% Shows delta-Entropies according to PSTH-letter based Entropy

clear;

pipeParams.subjID = 'Kramer';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anywaypipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz

% PSTH-bin parameters
pipeParams.binEdgeLeft = 0; % seconds
pipeParams.binWidth = 0.5/1000; % seconds
pipeParams.binNum = 15;
pipeParams.trimPSTH = true; % if true, the first bin in the PSTH is removed; if false, no change

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% Final figures to be saved somwhere?
pipeParams.finalFigs.save = true;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

runPipeline('EntropyPSTHletter_byCell', pipeParams);



%% EntropyDirectISI_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation.

clear; 

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 20; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = true;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

runPipeline('EntropyDirectISI_byCell', pipeParams);



%% EntropyDirectISI_PhLckVsPattMd_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation for 1) Phase-locked cells and 2) general pattern-modulation
% cells

clear; 

pipeParams.subjID = 'Kramer';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 20; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = true;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

runPipeline('EntropyDirectISI_PhLckVsPattMd_byCell', pipeParams);



%% EntropyDirectISIord1_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation.

clear; 

pipeParams.subjID = 'Kramer';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 15; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

runPipeline('EntropyDirectISIord1_byCell', pipeParams);



%% EntropyDirectISIord1_PhLckVsPattMd_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation for 1) Phase-locked cells and 2) general pattern-modulation
% cells

clear; 

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 20; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

runPipeline('EntropyDirectISIord1_PhLckVsPattMd_byCell', pipeParams);



%% EntropyDirectISIord1_PhLckVsNonPhLckVsAll_deltaH_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation for 1) Phase-locked cells and 2) general pattern-modulation
% cells

% *************************************************************************

clear; 

pipeParams.subjID = 'Uva';

% If a given intermediate table already exists, overwrite it anyway
% WARNING, may take a long time to redo some tables, hours
pipeParams.overwriteTables = false; % true/false

pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 15; % num of bins per decade of log-spaced bins for ISI
pipeParams.neuronMinimumTrials = 4;

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

sp = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190226_PatternModFigs\FsweepLowerOutliers';
pipeParams.finalFigs.savepn = sp;
runPipeline('EntropyDirectISIord1_PhLckVsNonPhLckVsAll_deltaH_byCell', pipeParams);



%% EntropyDirectISIord1_PhLckVsPattMdVsAll_deltaH_byCellONEplot pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation for 1) Phase-locked cells and 2) general pattern-modulation
% cells

clear; 

pipeParams.subjID = 'Kramer';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 15; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

runPipeline('EntropyDirectISIord1_PhLckVsPattMdVsAll_deltaH_byCellONEplot', pipeParams);



%% EntropyDirectISIord1_PhLckVsPattMdVsAll_dbsH_byCell pipeline
% Shows "raw" dbs-Entropies according to ISI-based (log-binned) Entropy
% estimation for 1) Phase-locked cells and 2) all cells

clear; 

pipeParams.subjID = 'Kramer';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 20; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

runPipeline('EntropyDirectISIord1_PhLckVsPattMdVsAll_dbsH_byCell', pipeParams);



%% Entropy PSTH-letter display of Phase-locked population pipeline

clear;

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anywaypipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz

% PSTH-bin parameters
pipeParams.binEdgeLeft = 0; % seconds
pipeParams.binWidth = 0.5/1000; % seconds
pipeParams.binNum = 15;
pipeParams.trimPSTH = true; % if true, the first bin in the PSTH is removed; if false, no change

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% Final figures to be saved somwhere
pipeParams.PhsLckPercentXLim = [0, 30];
pipeParams.PhsLckPercentYLim = [0, 50];

pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190122_Group figures';

runPipeline('EntropyPSTH_percentPopPhaseLocked', pipeParams);



%% EntropyDirectISIord1_All_userVarH pipeline
% Shows dbs-Entropies according to ISI-based (log-binned) Entropy
% estimation for three options:
% 1) dbsH (Calculated entropies for the during-DBS period of each neuron)
% 2) preH (Calculated entropies for the preDBS period of each neuron)
% 3) deltaH (dbsH - preH)


clear; 

pipeParams.subjID = 'Uva';

% If a given intermediate table already exists, overwrite it anyway
% WARNING, may take a long time to redo some tables, hours
pipeParams.overwriteTables = false; % true/false

pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 15; % num of bins per decade of log-spaced bins for ISI

% for each neuron observed, exclude it from analysis if it does not stick
% around for at least "neuronMininumTrials" of dbs-trials (i.e. 4 trials)
pipeParams.neuronMinimumTrials = 4; 

pipeParams.entropyType = '%deltaH'; % choose: 'dbsH', 'preH', 'deltaH', 'pctDeltaH'
pipeParams.trialType = 'FrequencySweep'; % choose: 'ContactSweep', 'FrequencySweep'

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

sp = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190226_PatternModFigs\FsweepLowerOutliers';
pipeParams.finalFigs.savepn = sp;
runPipeline('EntropyDirectISIord1_All_userVarH', pipeParams);

