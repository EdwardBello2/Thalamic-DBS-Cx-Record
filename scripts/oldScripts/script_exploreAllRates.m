% script with inputs to Rate Analysis pipelines


%% EntropyDirectISIord1_PhLckVsNonPhLckVsAll_deltaH_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation for 1) Phase-locked cells and 2) general pattern-modulation
% cells

% *************************************************************************

%% Include project folder: "Thalamic-DBS-Cx-Record" latest version from local git-repo

% toolboxPath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record';
toolboxPath = 'C:\Users\bello043\Documents\GitHub\Thalamic-DBS-Cx-Record';
addpath(genpath(toolboxPath));



%% Code

clear; 

pipeParams.subjID = 'Kramer';
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
pipeParams.intDataPn = ['L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\intermediateData\spkRate\', pipeParams.subjID];
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



%% Load relevant tables and join them for analysis
% rateData = buildAnalysisTable_Rates(pipeParams);

subjID = pipeParams.subjID;


% load NEXProcfiles table
load([pipeParams.tablepn, '\', 'NEXprocfiles_', subjID], 'NEXprocfiles');
N = NEXprocfiles;

% load SortedUnits table
load([pipeParams.tablepn, '\', 'SortedUnits_', subjID], 'SortedUnits');
U = SortedUnits;

% load SweepAnalysisTrials4Paper2018 table
load([pipeParams.tablepn, '\', 'SweepAnalysisTrials4Paper2018_', subjID], ...
                                          'SweepAnalysisTrials4Paper2018');
S = SweepAnalysisTrials4Paper2018;

% load RateBins_60sec_1secBins table
load([pipeParams.tablepn, '\', 'RateBins_60sec_1secBins_', subjID], 'RateBins');
R = RateBins;


% join tables to NEXprocfiles
N2 = join(N, U, 'LeftKeys', 2, 'RightKeys', 1);

N3 = join(N2, S, 'LeftKeys', 3, 'RightKeys', 1);

N4 = join(N3, R, 'LeftKeys', 6, 'RightKeys', 1);



%% Specify a selection of the above joined table for final analysis

% basically this is like a database query:
isSU = strcmp(N4.NeuronType, 'SU');
isFreqSwp = N4.FrequencySweep == 1;
isNotSparse = (N4.meanRatePRE >= 2) | (N4.meanRateDBS >= 2);
Nselect = N4((isSU & isFreqSwp & isNotSparse), :);


% Remove any trials where frequency sweep was not done on electrode C0
isC0 = strcmp(Nselect.dbsContact, 'C0');
Nselect(~isC0,:) = [];


% Remove any neurons from analysis that do not have at least 4 trials
% present in the current selection table

% count the number of times each neuron shows up and store in new table
% called "NeuronCounts"
neurons = Nselect.Unit_objectID;
[uNeurons, ~, uNeuIdx] = unique(neurons);

nNeurons = length(uNeurons); % number of unique neurons
neuCount = zeros(nNeurons, 1);
for iNeu = 1:nNeurons
    neuCount(iNeu) = sum(uNeuIdx == iNeu);
    
end

NeuronCounts = [table(uNeurons), table(neuCount)];

% join this table to current selection table
Nselect = join(Nselect, NeuronCounts, 'LeftKeys', 2, 'RightKeys', 1);

% remove those rows with under 4 trials
isPresentforTrials = (Nselect.neuCount >= 4);
NselectF = Nselect;
NselectF(~isPresentforTrials,:) = [];



%% cycle thru a display of each individual neuron's rate chagnes
% in response to various frequencies of stimulation



% get all unique neuron unitIDs

uniqueIDs = unique(NselectF.Unit_objectID);
numNeurons = length(uniqueIDs);
for iNeu = 1:numNeurons
    disp(iNeu);
    
    
    % get subselection of rateData table for iNeuron
    unitID = uniqueIDs{iNeu};

    isUnit = strcmp(unitID, NselectF.Unit_objectID);
    rateDataSU = NselectF(isUnit,:);
    
    
    % display rate responses
    f1 = exploreNeuronRatesFsweep(rateDataSU, 'blockNum');
    f2 = exploreNeuronRatesFsweep(rateDataSU, 'dbsFrequency');
    
    f1.Position = [1947 272 808 682];
    f2.Position = [2900 253 808 682];

    
    % hold for user input
    pause()
    close(f1)
    close(f2)
    
end


%% 
% runRateChangeAllTrials(pipeParams);

