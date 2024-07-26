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

binSeconds = 1;
nRows = size(Tselect, 1);
dbsRateChange = cell(nRows, 1);
dbsExciteRate = zeros(nRows, 1);
dbsInhibRate = zeros(nRows, 1);
dbsNoChangeRate = zeros(nRows, 1);

for iNex = 1:nRows
    % load the row's matfile with intermediate data
    matfn = Tselect.matFile_RateBins{iNex};
    matpn = Tselect.matFileFolder_RateBins{iNex};
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
    
    
Tselect = [Tselect, table(dbsExciteRate), table(dbsInhibRate), table(dbsNoChangeRate), ...
                                                   table(dbsRateChange)];





% %% PLOT boxplots of delta Entropy 
% 
%  % Gather all unique DBS Frequencies
% freqs = unique(Tselect.dbsFrequency);
% nFreqs = numel(freqs);
% 
% 
% % Gather up all rateChange types coutns for each frequency
% Fcount = zeros(nFreqs, 3);
% for i = 1:nFreqs
%     isHz = Tselect.dbsFrequency == freqs(i);
%     HzR = Tselect(isHz,:);
%     nExc = sum(HzR.dbsExciteRate);
%     nInb = sum(HzR.dbsInhibRate);
%     nNon = sum(HzR.dbsNoChangeRate);
% 
%     Fcount(i,1:3) = [nExc, nNon, nInb];
% 
% end
% 
% Fpercent = 100 * (Fcount ./ sum(Fcount, 2));
% % Fpercent(:,2) = 100 * (1 - (Fmod ./ Ftot));
% 
% 
% figure; ax1 = axes;
% b = bar(Fpercent, 'stacked');
% b(1).FaceColor = [1.0, 1.0, 1.0];
% b(2).FaceColor = [0.5, 0.5, 0.5];
% b(3).FaceColor = [0.0, 0.0, 0.0];
% set(gca, 'xticklabel', freqs)
% legend('Excite', 'noChange', 'Inhib', 'Location', 'northeastoutside')
% 
% title([ppar.subjID, ':  Rate-Change Categories']);
% 
% ylabel('% neurons recorded')
% xlabel('DBS Frequency (Hz)')
% ax1.YLim = [0, 100];
% 
% 
% 
% %% Perform Fisher's exact test on previous tabulated data 
% % to see if proportions among groups are significantly different
% 
% [tblAll, chi2All, pAll] = crosstab(Tselect.dbsRateChange, Tselect.dbsFrequency)
% 
% 
% % excited vs not-excited, across DBS frequencies:
% [tblExciteVsAll, chi2ExciteVsAll, pExciteVsAll] = crosstab(Tselect.dbsExciteRate, Tselect.dbsFrequency)
% 
% 
% % inhib'd vs not-inhib'd, across DBS frequencies:
% [tblInhibVsAll, chi2InhibVsAll, pInhibVsAll] = crosstab(Tselect.dbsInhibRate, Tselect.dbsFrequency)
% 
% 
% 
% % use Fisher's exact test to assess one frequency vs all the rest 
% for i = 1:nFreqs
%     isHz = Tselect.dbsFrequency == freqs(i);
% %     HzR = Tselect(isHz,:);
%     [tbl, chi2(i,1), pChi(i,1)] = crosstab(isHz, Tselect.dbsExciteRate);
%     [~, pFisher(i,1), statsFisher] = fishertest(tbl);
%     ciLoFisher(i,1) = statsFisher.ConfidenceInterval(1);
%     ciHiFisher(i,1) = statsFisher.ConfidenceInterval(2);
%     
%     
% 
% end
% ExciteTab = [table(freqs), table(chi2), table(pChi), table(pFisher), ...
%              table(ciLoFisher), table(ciHiFisher)]
%          
%          
% % use Fisher's exact test to assess one frequency vs all the rest 
% for i = 1:nFreqs
%     isHz = Tselect.dbsFrequency == freqs(i);
% %     HzR = Tselect(isHz,:);
%     [tbl, chi2(i,1), pChi(i,1)] = crosstab(isHz, Tselect.dbsInhibRate);
%     [~, pFisher(i,1), statsFisher] = fishertest(tbl);
%     ciLoFisher(i,1) = statsFisher.ConfidenceInterval(1);
%     ciHiFisher(i,1) = statsFisher.ConfidenceInterval(2);
%     
%     
% 
% end
% InhibTab = [table(freqs), table(chi2), table(pChi), table(pFisher), ...
%              table(ciLoFisher), table(ciHiFisher)]         
%         
%          
% % use Fisher's exact test to assess one frequency vs all the rest 
% for i = 1:nFreqs
%     isHz = Tselect.dbsFrequency == freqs(i);
% %     HzR = Tselect(isHz,:);
%     [tbl, chi2(i,1), pChi(i,1)] = crosstab(isHz, ~Tselect.dbsNoChangeRate);
%     [~, pFisher(i,1), statsFisher] = fishertest(tbl);
%     ciLoFisher(i,1) = statsFisher.ConfidenceInterval(1);
%     ciHiFisher(i,1) = statsFisher.ConfidenceInterval(2);
%     
%     
% 
% end
% EitherChangeTab = [table(freqs), table(chi2), table(pChi), table(pFisher), ...
%              table(ciLoFisher), table(ciHiFisher)]           
%          
% 

%% Calculate final Entropy analysis based on ppar.entropyType
% ppar.rateType = 'meanDiff'; % string, 'meanDBS' | 'meanDiff'

% Collect Entropy values for comparison
switch ppar.rateType
    case 'meanDBS'
%             EntropyChange_AllNeu = TbyGroup.H_DBSemp_Hisi;
        RateChange_AllNeu = Tselect.meanRateDBS;
        rateTitle = 'Mean Firing Rate';
        rateAxisLabel = 'spikes / s';

    case 'meanDiff'
%             EntropyChange_AllNeu = TbyGroup.H_DBSemp_Hisi - TbyGroup.H_PREbootAv;
        RateChange_AllNeu = Tselect.meanRateDBS - Tselect.meanRatePRE;
        rateTitle = '\Delta Mean Firing Rate';
        rateAxisLabel = '\Delta spikes / s';

%     case 'isiPercDiffH'
% %             EntropyChange_AllNeu = (TbyGroup.H_DBSemp_Hisi - TbyGroup.H_PREbootAv) ./ ...
% %                         TbyGroup.H_PREbootAv;
%         EntropyChange_AllNeu = (Tselect.H_DBSemp_Hisi - Tselect.H_PREbootAv) ./ ...
%                     Tselect.H_PREbootAv;
%         rateTitle = '%\DeltaEntropy';
%         rateAxisLabel = '%\Delta Entropy';

    otherwise
        error('wrong input for ppar.entropyType');

end



%%

% Display boxplots & ANOVA results of ALL cells 
f1 = figure;
ax1 = axes;
switch ppar.trialType
    
    case 'FrequencySweep'
        % Get Frequency labels of all neurons as array of strings
        freqs = Tselect.dbsFrequency;
        Rate_grpLabel = cellstr(num2str(freqs));
        Flabels_AllNeu = cellstr(num2str(unique(freqs)));
        
        dispDeltaH_Fswp_Boxplot(RateChange_AllNeu, Rate_grpLabel, ...
                                Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
       ylabel(rateAxisLabel)

    case 'ContactSweep'
        % Get Electrode Contact labels for all neurons as array of strings
        contacts = Tselect.dbsContact;
        Rate_grpLabel = contacts;
        Clabels_AllNeu = unique(contacts);
        
        dispDeltaH_Cswp_Boxplot(RateChange_AllNeu, Rate_grpLabel, ...
                                Clabels_AllNeu, 'Pres', 'cleanWoutliers');
        xlabel(rateAxisLabel)
        
    otherwise
        error('Oh snap! wrong trialType')
        
end

% f1.Position = FIG_POSITION;
title([ppar.subjID, ':  ', rateTitle]); 





% [pF, tabF, statsF] = anova1(Entropy_AllNeu_Fswp, convert2stringArray(R_ISI_Fswp.dbsFrequency(:)));
[pV, tab, stats] = kruskalwallis(RateChange_AllNeu, Rate_grpLabel);


% Final table summarizing the number of recordings for each condition:


% % Final "deltaH_cell" detailing data points according to group order in "label_cell" 
% labels_cell= unique(Rate_grpLabel);
%  
% nLabs = numel(labels_cell);
% for iLab = 1:nLabs
%     isLabel = strcmp(labels_cell{iLab}, Rate_grpLabel);
%     deltaH_cell{iLab} = RateChange_AllNeu(isLabel);
%     
% end



%% Now group LFS and HFS

switch ppar.trialType
    
    case 'FrequencySweep'
        % Choose the two DBS frequency labels from the data for comparison
        % strings: 10 | 20 | 30 | 50 | 100 | 130
        LFSfreq1 = ' 10'; 
        LFSfreq2 = ' 20';

        HFSfreq1 = '100';
        HFSfreq2 = '130';


        % Sift out all frequencies but the two choices
        isLFS1 = strcmp(LFSfreq1, Rate_grpLabel);
        isLFS2 = strcmp(LFSfreq2, Rate_grpLabel);

        isHFS1 = strcmp(HFSfreq1, Rate_grpLabel);
        isHFS2 = strcmp(HFSfreq2, Rate_grpLabel);

        % isInclude = isLFS1 | isLFS2 | isHFS1 | isHFS2;
        % 
        % H_include = Entropy_AllNeu(isInclude);
        % Freq_include = Entropy_grpLabel(isInclude);


        % combine lows and combine highs
        isLFS = isLFS1 | isLFS2;
        isHFS = isHFS1 | isHFS2;


        % Divide Entropy data between the two groups
        Hlfs = RateChange_AllNeu(isLFS);
        Hhfs = RateChange_AllNeu(isHFS);

        [hTtest, pTtest, ci, statsTtest] = ttest2(Hlfs, Hhfs)

        [pRanksum, hRanksum, statsRanksum] = ranksum(Hlfs, Hhfs)
        
    otherwise % do nothing, I don't have this figured out for ContactSweep case...
        
end