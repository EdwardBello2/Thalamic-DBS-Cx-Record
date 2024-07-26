% script to analyze rate-changes in individual trial basis



%% Include toolbox: "Thalamic-DBS-Cx-Record" latest version from local git-repo

toolboxPath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record';
addpath(genpath(toolboxPath));



%% Code

clear; 

pValThresh = 0.01;
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



%% Calculate Rate changes for selected rows

N = NselectF;


binSeconds = 1;
nRows = size(N, 1);
dbsRateChange = cell(nRows, 1);
dbsExciteRate = zeros(nRows, 1);
dbsInhibRate = zeros(nRows, 1);
dbsNoChangeRate = zeros(nRows, 1);

for iNex = 1:nRows
    % load the row's matfile with intermediate data
    matfn = N.ratebinFn{iNex};
    matpn = N.ratebinPn{iNex};
    counts = load([matpn, '\', matfn]);
    
    ratesPRE = counts.bins.pre / binSeconds;
    ratesDBS = counts.bins.dbs / binSeconds;
    
    % rescale rates with log-transform for more normal distribution
    ratesPRElog = log(ratesPRE + 1);
    ratesDBSlog = log(ratesDBS + 1);
    
    
    % perform stat test to see if significant rate change
    [pVal, h, stats] = ranksum(ratesPRE, ratesDBS);
%     [pVal, h, stats] = ranksum(ratesPRElog, ratesDBSlog);

    
    % label row appropriately
    if pVal <= pValThresh
        if median(ratesDBS) > median(ratesPRE)
            dbsExciteRate(iNex,1) = 1;
            dbsRateChange{iNex,1} = 'excite';
            
        else
            dbsInhibRate(iNex,1) = 1;
            dbsRateChange{iNex,1} = 'inhib';
            
        end
               
    else
        dbsNoChangeRate(iNex,1) = 1;
        dbsRateChange{iNex,1} = 'noChange';
    end
    
end
    
    
N = [N, table(dbsExciteRate), table(dbsInhibRate), table(dbsNoChangeRate), ...
                                                   table(dbsRateChange)];





%% PLOT boxplots of delta Entropy 

 % Gather all unique DBS Frequencies
freqs = unique(N.dbsFrequency);
nFreqs = numel(freqs);


% Gather up all rateChange types coutns for each frequency
Fcount = zeros(nFreqs, 3);
for i = 1:nFreqs
    isHz = N.dbsFrequency == freqs(i);
    HzR = N(isHz,:);
    nExc = sum(HzR.dbsExciteRate);
    nInb = sum(HzR.dbsInhibRate);
    nNon = sum(HzR.dbsNoChangeRate);

    Fcount(i,1:3) = [nExc, nNon, nInb];

end

Fpercent = 100 * (Fcount ./ sum(Fcount, 2));
% Fpercent(:,2) = 100 * (1 - (Fmod ./ Ftot));


figure; 
b = bar(Fpercent, 'stacked');
b(1).FaceColor = [1.0, 1.0, 1.0];
b(2).FaceColor = [0.5, 0.5, 0.5];
b(3).FaceColor = [0.0, 0.0, 0.0];
set(gca, 'xticklabel', freqs)
legend('Excite', 'noChange', 'Inhib', 'Location', 'northeastoutside')

title(['Rate-Changes in ', pipeParams.subjID, ' for ', pipeParams.trialType]);
ylabel('% neurons recorded')
xlabel('DBS Frequency (Hz)')



%% Perform Fisher's exact test on previous tabulated data 
% to see if proportions among groups are significantly different

[tblAll, chi2All, pAll] = crosstab(N.dbsRateChange, N.dbsFrequency)


% excited vs not-excited, across DBS frequencies:
[tblExciteVsAll, chi2ExciteVsAll, pExciteVsAll] = crosstab(N.dbsExciteRate, N.dbsFrequency)


% inhib'd vs not-inhib'd, across DBS frequencies:
[tblInhibVsAll, chi2InhibVsAll, pInhibVsAll] = crosstab(N.dbsInhibRate, N.dbsFrequency)



% use Fisher's exact test to assess one frequency vs all the rest 
for i = 1:nFreqs
    isHz = N.dbsFrequency == freqs(i);
%     HzR = N(isHz,:);
    [tbl, chi2(i,1), pChi(i,1)] = crosstab(isHz, N.dbsExciteRate);
    [~, pFisher(i,1), statsFisher] = fishertest(tbl);
    ciLoFisher(i,1) = statsFisher.ConfidenceInterval(1);
    ciHiFisher(i,1) = statsFisher.ConfidenceInterval(2);
    
    

end
ExciteTab = [table(freqs), table(chi2), table(pChi), table(pFisher), ...
             table(ciLoFisher), table(ciHiFisher)]
         
         
% use Fisher's exact test to assess one frequency vs all the rest 
for i = 1:nFreqs
    isHz = N.dbsFrequency == freqs(i);
%     HzR = N(isHz,:);
    [tbl, chi2(i,1), pChi(i,1)] = crosstab(isHz, N.dbsInhibRate);
    [~, pFisher(i,1), statsFisher] = fishertest(tbl);
    ciLoFisher(i,1) = statsFisher.ConfidenceInterval(1);
    ciHiFisher(i,1) = statsFisher.ConfidenceInterval(2);
    
    

end
InhibTab = [table(freqs), table(chi2), table(pChi), table(pFisher), ...
             table(ciLoFisher), table(ciHiFisher)]         
        
         
% use Fisher's exact test to assess one frequency vs all the rest 
for i = 1:nFreqs
    isHz = N.dbsFrequency == freqs(i);
%     HzR = N(isHz,:);
    [tbl, chi2(i,1), pChi(i,1)] = crosstab(isHz, ~N.dbsNoChangeRate);
    [~, pFisher(i,1), statsFisher] = fishertest(tbl);
    ciLoFisher(i,1) = statsFisher.ConfidenceInterval(1);
    ciHiFisher(i,1) = statsFisher.ConfidenceInterval(2);
    
    

end
EitherChangeTab = [table(freqs), table(chi2), table(pChi), table(pFisher), ...
             table(ciLoFisher), table(ciHiFisher)]           
         



%% Old code

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
% %% Pairwise 2-sample t-tests for each DBS frequency
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
% [result(iFreq,1), p(iFreq,1), ci, stats] = ttest2(ratesPre, ratesDbs);
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
%                    table(tstat), table(df), table(sd), table(sampNum)];
% 
% 
% 
% 
% 
% %% Old code
% 
% 
% % %% cycle thru a display of each individual neuron's rate chagnes
% % % in response to various frequencies of stimulation
% % 
% % 
% % 
% % % get all unique neuron unitIDs
% % 
% % uniqueIDs = unique(NselectF.Unit_objectID);
% % numNeurons = length(uniqueIDs);
% % for iNeu = 1:numNeurons
% %     disp(iNeu);
% %     
% %     
% %     % get subselection of rateData table for iNeuron
% %     unitID = uniqueIDs{iNeu};
% % 
% %     isUnit = strcmp(unitID, NselectF.Unit_objectID);
% %     rateDataSU = NselectF(isUnit,:);
% %     
% %     
% %     % display rate responses
% %     f1 = exploreNeuronRatesFsweep(rateDataSU, 'blockNum');
% %     f2 = exploreNeuronRatesFsweep(rateDataSU, 'dbsFrequency');
% %     
% %     f1.Position = [1947 272 808 682];
% %     f2.Position = [2900 253 808 682];
% % 
% %     
% %     % hold for user input
% %     pause()
% %     close(f1)
% %     close(f2)
% %     
% % end
% % 
% % 
% % %% 
% % % runRateChangeAllTrials(pipeParams);
% % 
