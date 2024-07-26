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
CURR_FUNC = 'script_displayRateChangeCategory_barplotStack'; 



%% LOAD all relevant tables and MERGE them

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



%% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);



%% Calculate Rate changes for selected rows

N = Tselect;


binSeconds = 1;
nRows = size(N, 1);
dbsRateChange = cell(nRows, 1);
dbsExciteRate = zeros(nRows, 1);
dbsInhibRate = zeros(nRows, 1);
dbsNoChangeRate = zeros(nRows, 1);

for iNex = 1:nRows
    % load the row's matfile with intermediate data
    matfn = N.matFile_RateBins{iNex};
    matpn = N.matFileFolder_RateBins{iNex};
    counts = load([ppar.projRootPath, '\', matpn, '\', matfn]);
    
    ratesPRE = counts.bins.pre / binSeconds;
    ratesDBS = counts.bins.dbs / binSeconds;
    
    % rescale rates with log-transform for more normal distribution
    ratesPRElog = log(ratesPRE + 1);
    ratesDBSlog = log(ratesDBS + 1);
    
    
    % perform stat test to see if significant rate change
    [pVal, h, stats] = ranksum(ratesPRE, ratesDBS);
%     [pVal, h, stats] = ranksum(ratesPRElog, ratesDBSlog);

    
    % label row appropriately
    if pVal <= ppar.pAlphaRateChange
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





%% PLOT barplots of % population rate changes (stacked)

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

Fpercent = 100 * (Fcount ./ sum(Fcount, 2))
% Fpercent(:,2) = 100 * (1 - (Fmod ./ Ftot));


figure; ax1 = axes;
b = bar(Fpercent, 'stacked');
b(1).FaceColor = [0.0, 0.0, 0.0];
b(2).FaceColor = [1.0, 1.0, 1.0];
b(3).FaceColor = [0.5, 0.5, 0.5];
set(gca, 'xticklabel', freqs)
legend('Excite', 'noChange', 'Inhib', 'Location', 'northeastoutside')

title([ppar.subjID, ':  Rate-Change Categories']);

ylabel('% neurons recorded')
xlabel('DBS Frequency (Hz)')
ax1.YLim = [0, 100];



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
