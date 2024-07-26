function runRateChangeAllTrials(pipeParams)
% Loads specifiied data from experiment based on table metadata, and plots
% results of rate changes as "excited", "inhibited", or "nochange"
%
% This pipeline assumes that the following tables exist and are correct:
%
% NEXprocfiles_subjID.mat
% SortedUnits_sibjID.mat
% SweepAnalysisTrials4Paper2018_subjID.mat
%
% where "subjID" is the nhp name (i.e. 'Uva')

% Author: Ed Bello
% Created: 2019/04/25
%
% TO-DO
% - Control for the case that there may only be 30 sec of DBS...
% - go back thru development and see if RateBins table creates rates or
%   counts...



%% DEFAULT PARAMETERS

DEFAULT.subjID        = 'XXX';
DEFAULT.tablepn       = '\';
DEFAULT.neuTypeFilter = 'SU';
DEFAULT.predbsTime    = 60; % seconds
DEFAULT.dbsTime       = 60; % seconds
DEFAULT.hzThresh      = 2; % Hz
% DEFAULT.ordH          = 2;
% DEFAULT.binsPD        = 20;
% DEFAULT.entropyType   = 'dbsH';
DEFAULT.trialType     = 'ContactSweep';


%% SET DEFAULT PARAMETERS IF USER HAS NOT SET THEM

% Check for specific User-defined pipeline parameter inputs:
if ~isfield(pipeParams, 'subjID'), pipeParams.subjID               = DEFAULT.subjID; end
if ~isfield(pipeParams, 'tablepn'), pipeParams.tablepn             = DEFAULT.tablepn; end
if ~isfield(pipeParams, 'neuTypeFilter'), pipeParams.neuTypeFilter = DEFAULT.neuTypeFilter; end
if ~isfield(pipeParams, 'predbsTime'), pipeParams.predbsTime       = DEFAULT.predbsTime; end
if ~isfield(pipeParams, 'dbsTime'), pipeParams.dbsTime             = DEFAULT.dbsTime; end
if ~isfield(pipeParams, 'hzThresh'), pipeParams.hzThresh           = DEFAULT.hzThresh; end
% if ~isfield(pipeParams, 'ordH'), pipeParams.ordH                   = DEFAULT.ordH ; end
% if ~isfield(pipeParams, 'binsPD'), pipeParams.binsPD               = DEFAULT.binsPD; end
% if ~isfield(pipeParams, 'entropyType'), pipeParams.entropyType     = DEFAULT.entropyType; end
if ~isfield(pipeParams, 'trialType'), pipeParams.trialType     = DEFAULT.trialType; end


%% LOAD NECESSARY METADATA TABLES

% Table names shortened for each of use in code

% load NEXprocfiles_XXX table, where metadata for individual NEX files is
% stored:
load([pipeParams.tablepn, '\', 'NEXprocfiles_', pipeParams.subjID, '.mat']);
NEX = NEXprocfiles; 

% load Sorted_Units_XXX table, where info related to Spike-Sorting is stored:
load([pipeParams.tablepn, '\', 'SortedUnits_', pipeParams.subjID, '.mat']);
Sort = SortedUnits;

% load SweepAnalysisTrials4Paper2018_XXX table, where info on DBS
% parameters is stored
load([pipeParams.tablepn, '\', 'SweepAnalysisTrials4Paper2018_', pipeParams.subjID, '.mat']);
TrialInfo = SweepAnalysisTrials4Paper2018;


%% As an intermediate step, read in all nexfiles and find their rate bins
% then save each as an intermediate processed data file

nNexfiles = size(NEX, 1);
intDataFn = cell(nNexfiles, 1);
intDataPn = cell(nNexfiles, 1);

% for iNex = 1:nNexfiles
%     % load one nexfile at a time
%     nexfn = NEX.Filename{iNex};
%     nexpn = NEX.Pathname{iNex};
%     nexFile = readNexFile([nexpn, '\', nexfn]);
% 
% 
%     % Get the spike times and DBS times
%     spkTimes = nexFile.neurons{1,1}.timestamps;
%     
%     evDBS = nexGetEvent(nexFile, 'DBS_stims');
%     dbsTimes = evDBS.timestamps;
%     
%     
%     % Create all rate bins for 1) dbs epoc, 2) predbs epoc
%     binWidth = 1; % second 
%     totTime = 60; % seconds
%     
%     dbsOnset = dbsTimes(1);
%     binEdgesDBS = dbsOnset:binWidth:(dbsOnset + totTime);
%     binCtsDBS = binTimeEvents(spkTimes, binEdgesDBS);
%     binRatesDBS = binCtsDBS / binWidth;
%     
% 
%     tStartPRE = dbsOnset - totTime;
%     tStartPRE = max(tStartPRE, nexFile.tbeg); % in case totTime exceeds amount of time in PRE condition
%     binEdgesPRE = dbsOnset:-binWidth:tStartPRE;
%     binCtsPRE = binTimeEvents(spkTimes, binEdgesPRE);
%     binRatesPRE = binCtsPRE / binWidth;
%     
% 
%     % Save intermediate mat file to intermediate processing location, and track
%     % this saved location in an "Intermediate" table
%     bins.dbs = binRatesDBS;
%     bins.pre = binRatesPRE;
%     
%     [~, nexID, ~] = fileparts(nexfn);
%     saveIntFn = [nexID, '_60sec_1secBins'];
%     saveIntPn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\intermediateData\spkRate';
%     save([saveIntPn, '\', saveIntFn],'bins');
%     
%     intDataFn{iNex,1} = saveIntFn;
%     intDataPn{iNex,1} = saveIntPn;
%     
% end
%     
% % Create Intermediate Data Table, with pathways for easily loading data
% 
% IFN = table(intDataFn);
% IPN = table(intDataPn);
% 
% RateBins = [NEX(:,1:3), IFN, IPN];
% 
% save([saveIntPn, '\', 'RateBins_60sec_1secBins'], 'RateBins');
    






%% CREATE/LOAD a table with SELECTED data to be analyzed for rate changes

nexLabel = 'NEXprocfiles';

Sort.Properties.VariableNames{1,1} = 'Unit_objectID';
NEXselect = join(NEX, Sort);
TrialInfo.Properties.VariableNames{1,1} = 'Trial_objectID';
NEXselect = join(NEXselect, TrialInfo);

% % SELECT only SU's 
% isSU = strcmp(Sort.NeuronType, 'SU');
% 
% SortSU = Sort(isSU,:);
% 
% % MERGE the tables where their primary keys are common
% isNexSU = ismember(NEX.Unit_objectID, SortSU.objectID);
% NexSU = NEX(isNexSU, :);

load([pipeParams.intDatapn, '\RateBins_60sec_1secBins']);

% SELECT only trials that had over 2Hz rate in at least one epoc
nNexfiles = size(RateBins, 1);
isAboveThresh = false(nNexfiles, 1);

tic
for iNex = 1:nNexfiles
%     clear counts
    % load intermediate data with rate bins: counts
    matfn = RateBins.intDataFn{iNex};
    matpn = RateBins.intDataPn{iNex};
    counts = load([matpn, '\', matfn]);
    
      
    % Get average spike rate for preDBS and DBSon:
    countsPRE = counts.bins.pre;
    avRatePRE = sum(countsPRE) / numel(countsPRE);
    
    countsDBS = counts.bins.dbs;
    avRateDBS = sum(countsDBS) / numel(countsDBS); % total spks / total seconds
    
    rateThresh = 2; %Hz
    isAboveThresh(iNex,1) = (avRatePRE > rateThresh) || (avRateDBS > rateThresh);
    
end
toc
% add boolean threshold check to RateBins table
RateBinsThresh = [RateBins, table(isAboveThresh)];


% MERGE Rate table with NEX-singleUnit table
rateData = join(NEXselect, RateBinsThresh);


% SELECT the subset of rows that pertain to SU's and above 2Hz rates
isSU = strcmp('SU', rateData.NeuronType);

isSelect = isSU & rateData.isAboveThresh; 

rateDataSelect = rateData(isSelect,:);



% 
% 
%     % for one nexFile, determine if BOTH pre and dbs are below 2Hz:
% 
%     % Get time intervals pertaining to dbs-on and pre-dbs
%     
%     
% PREDBS_TIME = pipeParams.predbsTime; % seconds
% DBS_TIME = pipeParams.dbsTime; % seconds
% HZ_THRESH = pipeParams.hzThresh; % Hz
% 
% % dbsOnset = nexFile
% 
% dbsEvent = nexGetEvent(nexFile, 'DBS_stims');
% dbsT = dbsEvent.timestamps; 
% dbsOnset = dbsT(1);
% 
% % specify pre and dbs intervals
% dbsInterv = [dbsOnset, dbsOnset + DBS_TIME];
% preInterv = [dbsOnset - PREDBS_TIME, dbsOnset];
% 
% % correct for the unlikely case that the specified PREDBS_TIME extends
% % further back in time than the file has (i.e. -20 seconds):
% preInterv(1) = max(preInterv(1), nexFile.tbeg);
% 
% 
% % Test whether both of these intervals BOTH have below 2Hz average spkRate:
% spkTimes = nexFile.neurons{1,1}.timestamps;
% intervBelowThrsh = detectIntervalSpkRateThresh(spkTimes, 2, [preInterv; dbsInterv], ...
%                             'ThreshCross', 'below');
%                         
% includeNexFile(iNex) = ~all(intervBelowThrsh);
%                         
% 
% 
% % MERGE the SU table with the above-2Hz-table
% 
% NEXsu = NEX(isNexSU,:);


% EXCLUDE instances where BOTH pre and dbs sections are <2Hz in spk-rate



%% CREATE/LOAD a table with Direct-Entropy estimates for both pre- and DBS periods
% 
% % Make name for analysis table depending on pipeParams
% baseLabel = 'EntropyDirectISI';
% name = buildEntropyDirectISITableName(pipeParams);
% H_ISItableName = [baseLabel, name];
% 
% 
% % Check if the intermediate table to be created already exists. If it does
% % not exist, runPipeline will proceed to create the table
% createEntropyDirectISITableIfNeeded(NEX_2analyze, H_ISItableName, pipeParams);
% 
% 
% % LOAD the specified Table: H_DirectISIResults
% load([pipeParams.tablepn, '\', H_ISItableName, '.mat']);
% 
% 
% 
% %% CREATE TABLE with PSTH-Entropy estimates for both pre- and DBS periods
% 
% % Make name for analysis table depending on pipeParams
% baseLabel = 'EntropyPsthLetter';
% name = buildEntropyPsthLetterTableName(pipeParams);
% H_PSTHtableName = [baseLabel, name];
% 
% 
% % Check if the intermediate table to be created already exists. If it does
% % not exist, runPipeline will proceed to create the table
% createEntropyPsthLetterTableIfNeeded(NEX_2analyze, H_PSTHtableName, pipeParams);
% 
% 
% % LOAD the specified table: H_letterResults
% load([pipeParams.tablepn, '\', H_PSTHtableName, '.mat' ]);
% 


%% SELECT sweep type: 'ContactSweep' | 'FrequencySweep'

switch pipeParams.trialType
    case 'ContactSweep'
        isCsweep = logical(rateDataSelect.ContactSweep);
        rateDataSweep = rateDataSelect(isCsweep,:);
        
        % Check to make sure only 130Hz is in Contact-Sweep
        freqs = rateDataSweep.dbsFrequency;
        isRemove = freqs ~= 130;
        rateDataSweep(isRemove, :) = [];
        
    case 'FrequencySweep'
        isFsweep = logical(rateDataSelect.FrequencySweep);
        rateDataSweep = rateDataSelect(isFsweep,:);
        
        % Check to make sure only C0 is in Frequency-Sweep
        cElec = rateDataSweep.dbsContact;
        isRemove = ~strcmp(cElec, 'C0');
        rateDataSweep(isRemove, :) = [];
        
    otherwise
        error('invalid input for pipeParams.trialType')
        
end



%% REMOVE individual neurons from analysis...
% ... if they are not present for AT LEAST "nTrialMin" number of trials

nTrialMin = pipeParams.neuronMinimumTrials;

% 1) Create tables to assess individual neuron presence for each DBS trial
% 2) Remove those neurons that don't have at least nTrialMin number of
% trials


% 1)
uniqueUnits = unique(rateDataSweep.Unit_objectID);
nUnits = numel(uniqueUnits);
disp(['unique SUs: ', num2str(nUnits)]);
unitCount = zeros(nUnits, 1);
for k = 1:nUnits
    unitCount(k) = sum(strcmp(rateDataSweep.Unit_objectID, uniqueUnits(k)));
    
end

UnitCounts = table(uniqueUnits, unitCount);
UnitCounts.Properties.VariableNames{1,1} = 'Unit_objectID';

% Add unitCount info to table; 
rateDataSweep = join(rateDataSweep, UnitCounts);


% 2) keep only those neurons that have minimum trial-presence
isUnitCountAbove = rateDataSweep.unitCount >= nTrialMin;
rateDataSweep = rateDataSweep(isUnitCountAbove,:);
nRemainingUniqueUnits = numel(unique(rateDataSweep.Unit_objectID));
disp(['unique SUs with at leat ', num2str(nTrialMin), ...
     ' trials: ', num2str(nRemainingUniqueUnits)]);
 
 
 
%% Calculate Rate changes for selected rows

binSeconds = 1;
nRows = size(rateDataSweep, 1);
dbsRateChange = cell(nRows, 1);
for iNex = 1:nRows
    % load the row's matfile with intermediate data
    matfn = RateBins.intDataFn{iNex};
    matpn = RateBins.intDataPn{iNex};
    counts = load([matpn, '\', matfn]);
    
    ratesPRE = counts.bins.pre / binSeconds;
    ratesDBS = counts.bins.dbs / binSeconds;
    
    % rescale rates with log-transform for more normal distribution
    ratesPRElog = log(ratesPRE + 1);
    ratesDBSlog = log(ratesDBS + 1);
    
    
    % perform stat test to see if significant rate change
    [pVal, h, stats] = ranksum(ratesPRElog, ratesDBSlog);

    
    % label row appropriately
    if pVal <= 0.05
        if median(ratesDBS) > median(ratesPRE)
            dbsRateChange{iNex,1} = 'excite';
            
        else
            dbsRateChange{iNex,1} = 'inhib';
            
        end
               
    else
        dbsRateChange{iNex,1} = 'noChange';
        
    end
    
end
    
    
rateDataSweep = [rateDataSweep, table(dbsRateChange)];



% 
% 
% % for Fsweep:
% fT = isNeuPresentForTrialFswp(R_ISI_Fswp);
% totTrials = sum(table2array(fT), 2);
% hasTrialMin = totTrials >= nTrialMin;
% neuIDs = fT.Properties.RowNames;
% keepNeurons = neuIDs(hasTrialMin);
% 
% % keep only those neurons that have minimum trial-presence
% nRows = size(R_ISI_Fswp, 1);
% keepRows = false(nRows, 1);
% for iRow = 1:nRows
%     keepRows(iRow,1) = any(strcmp(R_ISI_Fswp.Unit_objectID{iRow}, keepNeurons));
%     
% end
% R_ISI_Fswp = R_ISI_Fswp(keepRows,:);
% R_PSTH_Fswp = R_PSTH_Fswp(keepRows,:);
% 


%% EXTRACT rows containing all unique single units that have at least ONE pattern-mod trial
% Note: the two functions below also take unique neuron identities into
% account, so that all of one cell's rows are included in one table but not
% the other. 

% alpha = pipeParams.pValAlpha;
% 
% % for Csweep:
% [R_ISI_PhsLckNeu_Cswp, ...
%  R_ISI_NOTPhsLckNeu_Cswp] = separateISItable_PhsLckVsNOTPhsLck(R_ISI_Cswp, R_PSTH_Cswp, alpha);
% 
% 
% % for Fsweep:
% [R_ISI_PhsLckNeu_Fswp, ...
%  R_ISI_NOTPhsLckNeu_Fswp] = separateISItable_PhsLckVsNOTPhsLck(R_ISI_Fswp, R_PSTH_Fswp, alpha);
% 


%% CALCULATE Delta-Entropies for ISI-based Entropy of each trial

% deltaH_PhsLckNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_PhsLckNeu_Cswp);
% deltaH_NOTPhsLckNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_NOTPhsLckNeu_Cswp);
% Entropy_AllNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_Cswp);
% 
% deltaH_PhsLckNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_PhsLckNeu_Fswp);
% deltaH_NOTPhsLckNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_NOTPhsLckNeu_Fswp);
% Entropy_AllNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_Fswp);


%% GATHER Entropy values Entropy of each trial, according to:
% 1)   'dbsH': just gather Entropy values during DBS
% 2)   'preH': just gather Entropy values before DBS
% 3) 'deltaH': dbsH - preH

% switch pipeParams.entropyType
%     case 'dbsH'
%         % get during-DBS H1-estimate for all rows:
%         Entropy_AllNeu_Cswp = extractH1(R_ISI_Cswp.H_DBS);
%         Entropy_AllNeu_Fswp = extractH1(R_ISI_Fswp.H_DBS);
%         
%         
%     case 'preH'
%         % get pre-DBS H1-estimate for all rows:
%         Entropy_AllNeu_Cswp = extractH1pre(R_ISI_Cswp);
%         Entropy_AllNeu_Fswp = extractH1pre(R_ISI_Fswp);
%            
%         
%     case 'deltaH'
%         % get deltaH (dbsH - preH):
% %         deltaH_PhsLckNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_PhsLckNeu_Cswp);
% %         deltaH_NOTPhsLckNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_NOTPhsLckNeu_Cswp);
%         Entropy_AllNeu_Cswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_Cswp);
% 
% %         deltaH_PhsLckNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_PhsLckNeu_Fswp);
% %         deltaH_NOTPhsLckNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_NOTPhsLckNeu_Fswp);
%         Entropy_AllNeu_Fswp = calcDeltaEntropyISIord1_allTableRows(R_ISI_Fswp);
%         
%         
%     case '%deltaH'
%         deltaHC = calcDeltaEntropyISIord1_allTableRows(R_ISI_Cswp);
%         preHC = extractH1pre(R_ISI_Cswp);
%         Entropy_AllNeu_Cswp = 100*(deltaHC./preHC);
%         
%         deltaHF = calcDeltaEntropyISIord1_allTableRows(R_ISI_Fswp);
%         preHF = extractH1pre(R_ISI_Fswp);
%         Entropy_AllNeu_Fswp = 100*(deltaHF./preHF);        
%         
%         
%     otherwise % default case: 'dbsH'
%         error('Cannnot gather entropy values: invalid value for pipeParams.entropyType')
% 
%         
% end


%% GET DBS-trial labels for each row for plotting

% 
% [Clabels_PhsLckNeu, isRowClabel_PhsLckNeu] = getAlldbsElectrodeLabels(R_ISI_PhsLckNeu_Cswp);
% [Clabels_NOTPhsLckNeu, isRowClabel_NOTPhsLckNeu] = getAlldbsElectrodeLabels(R_ISI_NOTPhsLckNeu_Cswp);
% [Clabels_AllNeu, isRowClabel_AllNeu] = getAlldbsElectrodeLabels(R_ISI_Cswp);
% 
% 
% [Flabels_PhsLckNeu, isRowFlabel_PhsLckNeu] = getAlldbsFrequencyLabels(R_ISI_PhsLckNeu_Fswp);
% [Flabels_NOTPhsLckNeu, isRowFlabel_NOTPhsLckNeu] = getAlldbsFrequencyLabels(R_ISI_NOTPhsLckNeu_Fswp);
% [Flabels_AllNeu, isRowFlabel_AllNeu] = getAlldbsFrequencyLabels(R_ISI_Fswp);

% [Flabels, isRowFlabel] = getAlldbsFrequencyLabels(R_Fswp);


%% PLOT boxplots of delta Entropy 

R = rateDataSweep;

switch pipeParams.trialType
    case 'ContactSweep'
        % Gather all unique DBS Contacts
        contacts = unique(R.dbsContact);
        nConts = numel(contacts);
        
        % Gather up all rateChange types coutns for each frequency
        Ccount = zeros(nConts, 3);
        for i = 1:nConts
            isC = strcmp(contacts(i), R.dbsContact);
            cR = R(isC,:);
            nExc = sum(strcmp('excite', cR.dbsRateChange));
            nInb = sum(strcmp('inhib', cR.dbsRateChange));
            nNon = sum(strcmp('noChange', cR.dbsRateChange));
        
            Ccount(i,1:3) = [nExc, nNon, nInb];
        
        end
        
        Cpercent = 100 * (Ccount ./ sum(Ccount, 2));


        figure; 
        b = barh(Cpercent, 'stacked');
        b(1).FaceColor = [1.0, 1.0, 1.0];
        b(2).FaceColor = [0.5, 0.5, 0.5];
        b(3).FaceColor = [0.0, 0.0, 0.0];
        set(gca, 'yticklabel', contacts)
        set(gca, 'XLim', [0, 100])
        legend('Excite', 'nCh', 'Inhib', 'Location', 'northeastoutside')
        
        title(['Rate-Changes in ', pipeParams.subjID, ' for ', pipeParams.trialType]);
        xlabel('% neurons recorded')
        ylabel('DBS Contact')
        
        
        
        
        
        


        
    case 'FrequencySweep'
        % Gather all unique DBS Frequencies
        freqs = unique(R.dbsFrequency);
        nFreqs = numel(freqs);
        
        % Gather up all rateChange types coutns for each frequency
        Fcount = zeros(nFreqs, 3);
        for i = 1:nFreqs
            isHz = R.dbsFrequency == freqs(i);
            HzR = R(isHz,:);
            nExc = sum(strcmp('excite', HzR.dbsRateChange));
            nInb = sum(strcmp('inhib', HzR.dbsRateChange));
            nNon = sum(strcmp('noChange', HzR.dbsRateChange));
        
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
        
        
    otherwise
        error('Invalid value for pipeParams.trialType')
        
end
        



%% PRINT info about final dataset results

% % Print number of Neurons and number of data points for both "all neurons"
% % and just "phase-locked neurons", for Contact-sweep:
% c = newline;
% 
% disp('--------------- Contact Sweep ---------------');
% 
% 
% disp('All Neurons:');
% 
% nNeurons = numel(unique(R_ISI_Cswp.Unit_objectID));
% nDatapoints = size(R_ISI_Cswp, 1);
% disp(['    # Neurons: ', num2str(nNeurons)])
% disp(['# data-points: ', num2str(nDatapoints)]);
% disp(c)
% 
% 
% % disp('Phs-Lck Neurons:')
% % 
% % nNeurons = numel(unique(R_ISI_PhsLckNeu_Cswp.Unit_objectID));
% % nDatapoints = size(R_ISI_PhsLckNeu_Cswp, 1);
% % disp(['    # Neurons: ', num2str(nNeurons)])
% % disp(['# data-points: ', num2str(nDatapoints)]);
% % disp(c)
% 
% 
% disp('--------------- Frequency Sweep ---------------');
% 
% 
% disp('All Neurons:');
% 
% nNeurons = numel(unique(R_ISI_Fswp.Unit_objectID));
% nDatapoints = size(R_ISI_Fswp, 1);
% disp(['    # Neurons: ', num2str(nNeurons)])
% disp(['# data-points: ', num2str(nDatapoints)]);
% disp(c)
% 
% 
% % disp('Phs-Lck Neurons:')
% % 
% % nNeurons = numel(unique(R_ISI_PhsLckNeu_Fswp.Unit_objectID));
% % nDatapoints = size(R_ISI_PhsLckNeu_Fswp, 1);
% % disp(['    # Neurons: ', num2str(nNeurons)])
% % disp(['# data-points: ', num2str(nDatapoints)]);
% % disp(c)
% 
% 
% 
% 



%
%% Final optional saving of figures

if pipeParams.finalFigs.save == true
    savfPn = pipeParams.finalFigs.savepn;
    
    % save as .eps files
    saveas(f1, [savfPn, '\ISI_deltaH_PhsLck_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f2, [savfPn, '\ISI_deltaH_NotPhsLck_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f3, [savfPn, '\ISI_deltaH_All_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f4, [savfPn, '\ISI_deltaH_PhsLck_Fswp_', pipeParams.subjID ],'epsc');
    saveas(f5, [savfPn, '\ISI_deltaH_NotPhsLck_Fswp_', pipeParams.subjID ],'epsc');
    saveas(f6, [savfPn, '\ISI_deltaH_All_Fswp_', pipeParams.subjID ],'epsc');

    % save as .jpg files
    saveas(f1, [savfPn, '\ISI_deltaH_PhsLck_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f2, [savfPn, '\ISI_deltaH_NotPhsLck_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f3, [savfPn, '\ISI_deltaH_All_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f4, [savfPn, '\ISI_deltaH_PhsLck_Fswp_', pipeParams.subjID ],'jpg');
    saveas(f5, [savfPn, '\ISI_deltaH_NotPhsLck_Fswp_', pipeParams.subjID ],'jpg');
    saveas(f6, [savfPn, '\ISI_deltaH_All_Fswp_', pipeParams.subjID ],'jpg');

end





end % END function

function rowFlabel = convert2stringArray(rowFrequency)

% rowFrequency = R_ISI_PhsLckNeu_Fswp.dbsFrequency(:);
nRows = size(rowFrequency, 1);
rowFlabel = cell(nRows, 1);
for i = 1:nRows, rowFlabel{i,1} = num2str(rowFrequency(i,1)); end


end