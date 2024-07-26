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

% Tselect = Tselect;


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


%% If rateChangeType group is specified, whittle down table to have only mean rate values of that category
% 'excite' | 'inhib' | 'noChange'

rateChangeDataType = 'all'; % default value

if isfield(ppar, 'rateChangeType')
    str = ppar.rateChangeType;
    
    if strcmp(str, 'excite') | strcmp(str, 'inhib') | strcmp(str, 'noChange')
        isRateChangeType = strcmp(str, Tselect.dbsRateChange);
        Tselect = Tselect(isRateChangeType,:);
        rateChangeDataType = str;
        
    elseif strcmp(str, 'all')
        % do nothing
        
    else
        error('Oh snap! wrong input string for ppar.rateChangeType')
            
    end
     
end



%% Create User-defined Group column, with TF indices for which rows have data that belongs to said group

% ppar.groups.from.vars = {'dbsFrequency'};
% ppar.groups.from.select = {10};
% ppar.groups.newLabel = {'Hz10'};

if isfield(ppar.groups, 'from') & isfield(ppar.groups, 'newLabel')
    % Get data from user-specified existing column:
    varsData = Tselect{:, ppar.groups.from.vars{1}};

    % find user-defined selection from vars in table
    
    isUserGrp = false(size(varsData, 1), 1);
    for i = 1:numel(ppar.groups.from.select)
        isUserGrp = isUserGrp | (varsData == ppar.groups.from.select{i});
        
    end

    Grp = table(isUserGrp);
    Grp.Properties.VariableNames{1} = ppar.groups.newLabel{1};


    Tselect = [Tselect, Grp];

end




%% View Population average rates, grouped by DBS frequency
% ppar.groups.dbsFrequency.groupFreqs = {[10, 20], [100, 130]};
% ppar.groups.dbsFrequency.groupLabels = {'LFS', 'HFS'};

ratesDBSdata = [];
ratesDBSgrps = {};
if isfield(ppar.groups, 'from') % case where I specified subgroups of frequency (LFS/HFS)
    grps = ppar.groups.dbsFrequency.groupLabels;
    nFreqs = numel(grps);
    % barMeans = zeros(nFreqs, 2);
    % barStdvs = zeros(nFreqs, 2);
    barMeans = zeros(nFreqs, 1);
    barStdvs = zeros(nFreqs, 1);


    
    for iFreq = 1:nFreqs
        isGrp = strcmp(Tselect{:, 'group'}, grps{iFreq});

    %     ratesPre = Tselect.meanRatePRE(isFreq);
    %     barMeans(iFreq,1) = mean(ratesPre);
    %     barStdvs(iFreq,1) = std(ratesPre);

        ratesDbs = Tselect.meanRateDBS(isGrp);
        barMeans(iFreq,1) = mean(ratesDbs);
        barStdvs(iFreq,1) = std(ratesDbs);
        
        tempData = Tselect.meanRateDBS(isGrp);
        ratesDBSdata = [ratesDBSdata; tempData];
        tempGrps = cell(size(tempData, 1), 1); 
        tempGrps(:) = {grps{iFreq}};
        ratesDBSgrps = [ratesDBSgrps; tempGrps];

    end
    
    figure;
    h = barwitherr(barStdvs, barMeans);% Plot with errorbars

    set(gca,'XTickLabel',grps)
    
%     % data for the statistical test:
%     freqs = Tselect{:, 'dbsFrequency'};
%     grpLabel = cellstr(num2str(freqs));
    
    
else % case where the default 6 groups are the groups
    freqs = unique(Tselect.dbsFrequency);
    nFreqs = numel(freqs);
    % barMeans = zeros(nFreqs, 2);
    % barStdvs = zeros(nFreqs, 2);
    barMeans = zeros(nFreqs, 1);
    barStdvs = zeros(nFreqs, 1);

    for iFreq = 1:nFreqs
        isFreq = Tselect.dbsFrequency == freqs(iFreq);

    %     ratesPre = Tselect.meanRatePRE(isFreq);
    %     barMeans(iFreq,1) = mean(ratesPre);
    %     barStdvs(iFreq,1) = std(ratesPre);

        ratesDbs = Tselect.meanRateDBS(isFreq);
        barMeans(iFreq,1) = mean(ratesDbs);
        barStdvs(iFreq,1) = std(ratesDbs);

    end
    
    figure;
    h = barwitherr(barStdvs, barMeans);% Plot with errorbars

    set(gca,'XTickLabel',freqs)
    
    % data for the statistical test:
    freqs = Tselect{:, 'dbsFrequency'};
    ratesDBSgrps = cellstr(num2str(freqs));
    ratesDBSdata = Tselect.meanRateDBS;
    
end

ylabel('Spikes / second')
xlabel('DBS frequency (all @ C0)')
title([ppar.subjID, ':  Rate Comparison of DBS-on spkRates for ', ...
       rateChangeDataType, ' neural responses'])
% set(h(1),'FaceColor','k');

%             % get sub-selection of group labels and entropy values 
%             % according to user-grouping
%             idxGrp = Tfinal{:, ppar.groups.newLabel{1}};
%             freqs = Tfinal{idxGrp, 'dbsFrequency'};
%             
%             % Whittle down the Entropy results to user-defined subselection
%             EntropyChange_AllNeu = EntropyChange_AllNeu(idxGrp);
%             Entropy_grpLabel = cell(size(EntropyChange_AllNeu, 1), 1);
%             Entropy_grpLabel(:) = {ppar.groups.newLabel{1}};
%             
%                        
%         else
%             freqs = Tfinal{:, 'dbsFrequency'};
%             Entropy_grpLabel = cellstr(num2str(freqs));
% 
%              
%         end
%         
%         Flabels_AllNeu = unique(Entropy_grpLabel);
%         
%         dispDeltaH_Fswp_Boxplot(EntropyChange_AllNeu, Entropy_grpLabel, ...
%                                 Flabels_AllNeu, 'Pres', 'cleanWoutliers');                        
%         ylabel(entropyAxisLabel)

%% Perform ANOVA to see if group means significantly differ



[pF, tabF, statsF] = anova1(ratesDBSdata, ratesDBSgrps)



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
% 

