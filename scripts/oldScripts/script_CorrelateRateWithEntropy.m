% script to analyze rate-changes



%% Include toolbox: "Thalamic-DBS-Cx-Record" latest version from local git-repo

% toolboxPath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record';
toolboxPath = 'C:\Users\bello043\Documents\GitHub\Thalamic-DBS-Cx-Record';
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

% load EntropyDirectISI_<subjID>_SU_2HzThresh_60pre_60dbs_ordH2_15bins
load([pipeParams.tablepn, '\', 'EntropyDirectISI_', subjID, '_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots']);
H = H_DirectISIResults(:,{'objectID', 'H_PREemp', 'H_PREboot', 'H_DBS', 'pVal'});

% join tables to NEXprocfiles
N2 = join(N, U, 'LeftKeys', 2, 'RightKeys', 1);

N3 = join(N2, S, 'LeftKeys', 3, 'RightKeys', 1);

N4 = join(N3, R, 'LeftKeys', 6, 'RightKeys', 1);

Ncombo = join(H, N4, 'LeftKeys', 1, 'RightKeys', 1); % for now, just take advantage that both tables have identical Primary Keys...



%% Specify a selection of the above joined table for final analysis

% basically this is like a database query:
isSU = strcmp(Ncombo.NeuronType, 'SU');
isFreqSwp = Ncombo.FrequencySweep == 1;
isNotSparse = (Ncombo.meanRatePRE >= 2) | (Ncombo.meanRateDBS >= 2);
Nselect = Ncombo((isSU & isFreqSwp & isNotSparse), :);


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
Nselect = join(Nselect, NeuronCounts, 'LeftKeys', 6, 'RightKeys', 1);

% remove those rows with under 4 trials
isPresentforTrials = (Nselect.neuCount >= 4);
NselectF = Nselect;
NselectF(~isPresentforTrials,:) = [];



%% ScatterPlot Rate vs Entropy


for i = 1:size(NselectF, 1)
    H_DBS = NselectF.H_DBS{i,1};
    H1_DBS(i,1) = H_DBS(1,1);
    
end


figure; scatter(NselectF.meanRateDBS, H1_DBS)
xlabel('Average spikes / second');
ylabel('Entropy (bits/spk)');
title([subjID, ': DBS-ON spk-Rate vs. Entropy'])








% %% View Population average rates, grouped by DBS frequency
% 
% N = NselectF;
% 
% % Initialize matrix of dbsFrequency x (PRE or DBS)
% freqs = unique(N.dbsFrequency);
% nFreqs = 6;
% barMeans = zeros(nFreqs, 2);
% barStdvs = zeros(nFreqs, 2);
% 
% for iFreq = 1:nFreqs
%     isFreq = N.dbsFrequency == freqs(iFreq);
%     
%     ratesPre = N.meanRatePRE(isFreq);
%     barMeans(iFreq,1) = mean(ratesPre);
%     barStdvs(iFreq,1) = std(ratesPre);
%     
%     ratesDbs = N.meanRateDBS(isFreq);
%     barMeans(iFreq,2) = mean(ratesDbs);
%     barStdvs(iFreq,2) = std(ratesDbs);
%     
% end
%     
% 
% % Symmetric Example:
% % y = randn(3,4);         % random y values (3 groups of 4 parameters) 
% % errY = 0.1.*y;          % 10% error
% figure;
% h = barwitherr(barStdvs, barMeans);% Plot with errorbars
% 
% set(gca,'XTickLabel',freqs)
% legend('PRE-DBS','DBS-on')
% ylabel('Spikes / second')
% % set(h(1),'FaceColor','k');
% 
% 
% % figure; bar(N.meanRateDBS, N.dbsFrequency)
% % figure; boxplot(N.meanRatePRE, N.dbsFrequency)
% 
% 
% 
% %% Paired t-tests for each DBS frequency
% 
% result = zeros(nFreqs, 1);
%      p = zeros(nFreqs, 1);
%  ciLow = zeros(nFreqs, 1);
% ciHigh = zeros(nFreqs, 1);
%  tstat = zeros(nFreqs, 1);
%     df = zeros(nFreqs, 1);
%     sd = zeros(nFreqs, 1);
% sampNum = zeros(nFreqs, 1);
% 
% for iFreq = 1:nFreqs
%     isFreq = N.dbsFrequency == freqs(iFreq);
%     
%     ratesPre = N.meanRatePRE(isFreq);
% %     barMeans(iFreq,1) = mean(ratesPre);
% %     barStdvs(iFreq,1) = std(ratesPre);
%     
%     ratesDbs = N.meanRateDBS(isFreq);
% %     barMeans(iFreq,2) = mean(ratesDbs);
% %     barStdvs(iFreq,2) = std(ratesDbs);
% 
% if numel(ratesPre) ~= numel(ratesDbs)
%     disp('PRE and DBS samples not equal!')
%     
% end
% 
% disp(['for ', num2str(freqs(iFreq)), 'Hz, 2-sample ttest results:']);
% 
% [result(iFreq,1), p(iFreq,1), ci, stats] = ttest(ratesDbs, ratesPre);
% 
% ciLow(iFreq,1) = ci(1);
% ciHigh(iFreq,1) = ci(2);
% tstat(iFreq,1) = stats.tstat;
% df(iFreq,1) = stats.df;
% sd(iFreq,1) = stats.sd;
% sampNum(iFreq,1) = numel(ratesPre);
% 
% 
%     
% end
% 
% ResTtest = [table(result), table(p), table(ciLow), table(ciHigh), ...
%                    table(tstat), table(df), table(sd), table(sampNum)]
% 
% 
% 
% 

