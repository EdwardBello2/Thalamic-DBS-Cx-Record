% script for initializing parameters to affect the overall behavior of the
% pipeline


% full path on local PC where project folder is (don't include subfolders here,
% that's tracked within the appropriate tables)
ppar.projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record'; % string


% full path on local PC where tables to be loaded are kept
% ppar.tablePath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string
ppar.tablePath = 'C:\Users\bello043\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string


% Table selection related parameters
ppar.subjID = 'Kramer'; % string, 'Kramer' | 'Uva' | 'bothSubj'
ppar.trialType = 'ContactSweep'; % string, 'ContactSweep' | 'FrequencySweep'
ppar.neuronType = 'SU'; % string, 'SU' | 'MU' | 
ppar.hzThresh = 2; % Hz, numeric
ppar.neuronMinimumTrials = 3; % integer, minimum number of DBS trials that a unique neuron must be present for; otherwise all neuron data will be excluded from analysis.


% NEXfile related parameters
ppar.preDbsTime = 60; % seconds
ppar.dbsTime = 60; % seconds


% Row re-ordering related parameters
ppar.sortRows.byColumn = 'pVal_Hpsth'; % string, any variableName in table
ppar.sortRows.sortType = 'ascend'; % string, 'ascend' | 'descend'


% Data-grouping related parameters
% for scatter and pieChart(need to update)
ppar.groups.dbsFrequency.groupFreqs = {[10, 20], [100, 130]};
ppar.groups.dbsFrequency.groupLabels = {'LFS', 'HFS'};

% for boxplot & barplot
% ppar.groups.from.vars = {'dbsFrequency'};
% ppar.groups.from.select = {[10, 20], [100, 130]};
% ppar.groups.newLabel = {'LFS', 'HFS'};



% Rate-related parameters
ppar.rateType = 'meanDiff'; % string, 'meanDBS' | 'meanDiff'
ppar.pAlphaRateChange = 0.05; % double, alpha (p-value) for differentiating significant rate change categories
ppar.rateChangeType = 'all'; % string, 'excite' | 'inhib' | 'noChange' | 'all'


% PSTH-related parameters
ppar.trimPSTH = true; % TF indicating whether to remove the first and last bins
ppar.psthTimeBeg = 0; % seconds
ppar.psthBinWidth = 0.5 / 1000; %seconds
ppar.psthNumBins = 15;
ppar.pAlphaPSTH = 0.05;


% % Entropy-related parameters
ppar.entropyType = 'isiPercDiffH'; % string, 'isiDbsH' | 'isiPreH' | 'isiDiffH' | 'isiPercDiffH'


% log-transform related parameters (for multiple entries add another
% element to cell array {'one', 'two', ...}
% ppar.logTransform = {'H_DBSemp_Hisi', 'H_PREbootAv'}; % string, 'H_DBSemp_Hisi' | 'H_PREbootAv' | 'finalH'


% log-ISI related parameters
ppar.binsPerDecade = 15;
ppar.pAlphaISI = 0.05 / (2 * 6);


% Final display parameters
% ppar.boxplot.dataLim = [-0.10, 0.15];
ppar.dispPhaseLocked = 'both'; % 'phaseLck' | 'nonPhaseLck' | 'both'
ppar.boxplot.dataLim = [-0.10, 0.15]; % for %DiffH
% ppar.boxplot.dataLim = [-Inf, Inf]; % for Dbs/PreH
% ppar.boxplot.dataLim = [-15, 15]; % for diffRates



% Statistical test 
ppar.statTest = 'anova'; % string, 'anova' | 'kruskalwallis' | 'ttest' | 



