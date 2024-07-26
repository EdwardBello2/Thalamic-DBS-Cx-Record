% Display a scatter of rate-change vs entropy-change 

% Author: Ed Bello
% Created: 2019/07/17

%% pipeline Parameters and script Constants

% EXAMPLE PIPELINE PARAMETERS:
% % full path on local PC where project folder is (don't include subfolders here,
% % that's tracked within the appropriate tables)
% ppar.projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record'; % string
% 
% % full path on local PC where tables to be loaded are kept
% ppar.tablePath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string
% 
% % Table selection related parameters
% ppar.subjID = 'Kramer'; % string, 'Kramer' | 'Uva'
% ppar.trialType = 'FrequencySweep'; % string, 'ContactSweep' | 'FrequencySweep'
% ppar.neuronType = 'SU'; % string, SU | MU | 
% ppar.hzThresh = 2; % Hz, numeric
% ppar.neuronMinimumTrials = 4; % integer
% 
% % NEXfile related parameters
% ppar.preDbsTime = 60; % seconds
% ppar.dbsTime = 60; % seconds
% 
% % PSTH-related parameters
% ppar.trimPSTH = true; % TF indicating whether to remove the first and last bins
% ppar.psthTimeBeg = 0; % seconds
% ppar.psthBinWidth = 0.5 / 1000; %seconds
% ppar.psthNumBins = 15;
% ppar.pAlphaPSTH = 0.05;
% 
% % log-ISI related parameters
% ppar.binsPerDecade = 15;

clear; 

% Initialize fields in ppar struct in an accompanying script:
script_pipelineParams

 % CONSTANTS
CURR_FUNC = 'script_displayPopRates_PREandDBS_barplot'; 



%% LOAD all relevant tables and MERGE them

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



%% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);



%% View Population average rates, grouped by DBS frequency

N = Tselect;

% Initialize matrix of dbsFrequency x (PRE or DBS)
freqs = unique(N.dbsFrequency);
nFreqs = 6;
barMeans = zeros(nFreqs, 2);
barStdvs = zeros(nFreqs, 2);

for iFreq = 1:nFreqs
    isFreq = N.dbsFrequency == freqs(iFreq);
    
    ratesPre = N.meanRatePRE(isFreq);
    barMeans(iFreq,1) = mean(ratesPre);
    barStdvs(iFreq,1) = std(ratesPre);
    
    ratesDbs = N.meanRateDBS(isFreq);
    barMeans(iFreq,2) = mean(ratesDbs);
    barStdvs(iFreq,2) = std(ratesDbs);
    
end
    

% Symmetric Example:
% y = randn(3,4);         % random y values (3 groups of 4 parameters) 
% errY = 0.1.*y;          % 10% error
figure;
h = barwitherr(barStdvs, barMeans);% Plot with errorbars

set(gca,'XTickLabel',freqs)
legend('PRE-DBS','DBS-on')
ylabel('Spikes / second')
xlabel('DBS frequency (all @ C0)')

title([ppar.subjID, ':  Rate Comparison'])
% set(h(1),'FaceColor','k');


% figure; bar(N.meanRateDBS, N.dbsFrequency)
% figure; boxplot(N.meanRatePRE, N.dbsFrequency)



%% Paired t-tests for each DBS frequency

result = zeros(nFreqs, 1);
     p = zeros(nFreqs, 1);
 ciLow = zeros(nFreqs, 1);
ciHigh = zeros(nFreqs, 1);
 tstat = zeros(nFreqs, 1);
    df = zeros(nFreqs, 1);
    sd = zeros(nFreqs, 1);
sampNum = zeros(nFreqs, 1);

for iFreq = 1:nFreqs
    isFreq = N.dbsFrequency == freqs(iFreq);
    
    ratesPre = N.meanRatePRE(isFreq);
%     barMeans(iFreq,1) = mean(ratesPre);
%     barStdvs(iFreq,1) = std(ratesPre);
    
    ratesDbs = N.meanRateDBS(isFreq);
%     barMeans(iFreq,2) = mean(ratesDbs);
%     barStdvs(iFreq,2) = std(ratesDbs);

if numel(ratesPre) ~= numel(ratesDbs)
    disp('PRE and DBS samples not equal!')
    
end

disp(['for ', num2str(freqs(iFreq)), 'Hz, 2-sample ttest results:']);

[result(iFreq,1), p(iFreq,1), ci, stats] = ttest(ratesDbs, ratesPre);

ciLow(iFreq,1) = ci(1);
ciHigh(iFreq,1) = ci(2);
tstat(iFreq,1) = stats.tstat;
df(iFreq,1) = stats.df;
sd(iFreq,1) = stats.sd;
sampNum(iFreq,1) = numel(ratesPre);


    
end

ResTtest = [table(result), table(p), table(ciLow), table(ciHigh), ...
                   table(tstat), table(df), table(sd), table(sampNum)]





%% Old code


% %% cycle thru a display of each individual neuron's rate chagnes
% % in response to various frequencies of stimulation
% 
% 
% 
% % get all unique neuron unitIDs
% 
% uniqueIDs = unique(NselectF.Unit_objectID);
% numNeurons = length(uniqueIDs);
% for iNeu = 1:numNeurons
%     disp(iNeu);
%     
%     
%     % get subselection of rateData table for iNeuron
%     unitID = uniqueIDs{iNeu};
% 
%     isUnit = strcmp(unitID, NselectF.Unit_objectID);
%     rateDataSU = NselectF(isUnit,:);
%     
%     
%     % display rate responses
%     f1 = exploreNeuronRatesFsweep(rateDataSU, 'blockNum');
%     f2 = exploreNeuronRatesFsweep(rateDataSU, 'dbsFrequency');
%     
%     f1.Position = [1947 272 808 682];
%     f2.Position = [2900 253 808 682];
% 
%     
%     % hold for user input
%     pause()
%     close(f1)
%     close(f2)
%     
% end
% 
% 
% %% 
% % runRateChangeAllTrials(pipeParams);
% 
