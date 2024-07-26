function [rateDataSweep] = genTable_intDataRateBins_XXsec_XsecBins_subjID(parentTable, pipeParams)
% Generates table of intermediate data where spike times from individual
% nex file records are rate-binned in time. 
%
% Arranges all relevant data to be analyzed into one large table, with
% accompanying metadata.
% 
% This pipeline assumes that the following tables exist and are correct:
%
% NEXprocfiles_subjID.mat
% SortedUnits_sibjID.mat
% SweepAnalysisTrials4Paper2018_subjID.mat
%
% where "subjID" is the nhp name (i.e. 'Uva')

% Author: Ed Bello
% Created: 2019/06/28
%
% TO-DO
% 



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

% % Table names shortened for each of use in code
% 
% % load NEXprocfiles_XXX table, where metadata for individual NEX files is
% % stored:
% load([pipeParams.tablepn, '\', 'NEXprocfiles_', pipeParams.subjID, '.mat']);
% NEX = NEXprocfiles; 
% 
% % load Sorted_Units_XXX table, where info related to Spike-Sorting is stored:
% load([pipeParams.tablepn, '\', 'SortedUnits_', pipeParams.subjID, '.mat']);
% Sort = SortedUnits;
% 
% % load SweepAnalysisTrials4Paper2018_XXX table, where info on DBS
% % parameters is stored
% load([pipeParams.tablepn, '\', 'SweepAnalysisTrials4Paper2018_', ...
%       pipeParams.subjID, '.mat']);
% TrialInfo = SweepAnalysisTrials4Paper2018;

% pipeParams.ratebinWidth = 1; % seconds
% pipeParams.totTime = 60; % seconds

%% As an intermediate step, read in all nexfiles and find their rate bins
% then save each as an intermediate processed data file

nNexfiles = size(parentTable, 1);
intDataFn = cell(nNexfiles, 1);
intDataPn = cell(nNexfiles, 1);


% Create all rate bins for 1) dbs epoc, 2) predbs epoc
binWidth = pipeParams.ratebinWidth; % second 
totTime = pipeParams.totTime; % seconds


IntTableName = ['intDataRateBins_', num2str(totTime), 'sec_', num2str(binWidth), ...
                'secBins_', pipeParams.subjID];
            
% primary key of the new table should be identical to the parent table here
primKey = parentTable(:,'objectID');
% if ~exist([pipeParams.intDataPn, '\', IntTableName, '.mat']) 

for iNex = 1:nNexfiles
    % load one nexfile at a time
    nexfn = parentTable.Filename{iNex};
    nexpn = parentTable.Pathname{iNex};
    nexFile = readNexFile([nexpn, '\', nexfn]);
    pipeParams.tbeg = nexFile.tbeg;
    pipeParams.tend = nexFile.tend;

    % Get the spike times and DBS times
    spkTimes = nexFile.neurons{1,1}.timestamps;

    evDBS = nexGetEvent(nexFile, 'DBS_stims');
    dbsTimes = evDBS.timestamps;


    % Gather up spk times into bins using the function below:
    bins = genIntData_spkRate_XXsec_XsecBins(spkTimes, dbsTimes, pipeParams);

    
    % Save the data into an individual intermediate data file:
    [~, nexID, ~] = fileparts(nexfn);
    saveIntFn = [nexID, '_60sec_1secBins'];
    save([pipeParams.intDataPn, '\', saveIntFn],'bins');

    intDataFn{iNex,1} = saveIntFn;
    intDataPn{iNex,1} = pipeParams.intDataPn;

end

% Create Intermediate Data Table, with pathways for easily loading data

IFN = table(intDataFn);
IPN = table(intDataPn);


RateBins = [primKey, IFN, IPN];

save([pipeParams.tablepn, '\', IntTableName], 'RateBins');




% %% CREATE/LOAD a table with SELECTED data to be analyzed for rate changes
% 
% nexLabel = 'NEXprocfiles';
% 
% Sort.Properties.VariableNames{1,1} = 'Unit_objectID';
% NEXselect = join(NEX, Sort);
% TrialInfo.Properties.VariableNames{1,1} = 'Trial_objectID';
% NEXselect = join(NEXselect, TrialInfo);
% 
% % % SELECT only SU's 
% % isSU = strcmp(Sort.NeuronType, 'SU');
% % 
% % SortSU = Sort(isSU,:);
% % 
% % % MERGE the tables where their primary keys are common
% % isNexSU = ismember(NEX.Unit_objectID, SortSU.objectID);
% % NexSU = NEX(isNexSU, :);
% 
% % load([pipeParams.intDataPn, '\RateBins_60sec_1secBins']);
% 
% % SELECT only trials that had over 2Hz rate in at least one epoc
% nNexfiles = size(RateBins, 1);
% isAboveThresh = false(nNexfiles, 1);
% 
% tic
% for iNex = 1:nNexfiles
% %     clear counts
%     % load intermediate data with rate bins: counts
%     matfn = RateBins.intDataFn{iNex};
%     matpn = RateBins.intDataPn{iNex};
%     counts = load([matpn, '\', matfn]);
%     
%       
%     % Get average spike rate for preDBS and DBSon:
%     countsPRE = counts.bins.pre;
%     avRatePRE = sum(countsPRE) / numel(countsPRE);
%     
%     countsDBS = counts.bins.dbs;
%     avRateDBS = sum(countsDBS) / numel(countsDBS); % total spks / total seconds
%     
%     rateThresh = 2; %Hz
%     isAboveThresh(iNex,1) = (avRatePRE > rateThresh) || (avRateDBS > rateThresh);
%     
% end
% toc
% % add boolean threshold check to RateBins table
% RateBinsThresh = [RateBins, table(isAboveThresh)];
% 
% 
% % MERGE Rate table with NEX-singleUnit table
% rateData = join(NEXselect, RateBinsThresh);
% 
% 
% % SELECT the subset of rows that pertain to SU's and above 2Hz rates
% isSU = strcmp('SU', rateData.NeuronType);
% 
% isSelect = isSU & rateData.isAboveThresh; 
% 
% rateDataSelect = rateData(isSelect,:);
% 
% 
% 
% %% SELECT sweep type: 'ContactSweep' | 'FrequencySweep'
% 
% switch pipeParams.trialType
%     case 'ContactSweep'
%         isCsweep = logical(rateDataSelect.ContactSweep);
%         rateDataSweep = rateDataSelect(isCsweep,:);
%         
%         % Check to make sure only 130Hz is in Contact-Sweep
%         freqs = rateDataSweep.dbsFrequency;
%         isRemove = freqs ~= 130;
%         rateDataSweep(isRemove, :) = [];
%         
%     case 'FrequencySweep'
%         isFsweep = logical(rateDataSelect.FrequencySweep);
%         rateDataSweep = rateDataSelect(isFsweep,:);
%         
%         % Check to make sure only C0 is in Frequency-Sweep
%         cElec = rateDataSweep.dbsContact;
%         isRemove = ~strcmp(cElec, 'C0');
%         rateDataSweep(isRemove, :) = [];
%         
%     otherwise
%         error('invalid input for pipeParams.trialType')
%         
% end
% 
% 
% 
% %% REMOVE individual neurons from analysis...
% % ... if they are not present for AT LEAST "nTrialMin" number of trials
% 
% nTrialMin = pipeParams.neuronMinimumTrials;
% 
% % 1) Create tables to assess individual neuron presence for each DBS trial
% % 2) Remove those neurons that don't have at least nTrialMin number of
% % trials
% 
% 
% % 1)
% uniqueUnits = unique(rateDataSweep.Unit_objectID);
% nUnits = numel(uniqueUnits);
% disp(['unique SUs: ', num2str(nUnits)]);
% unitCount = zeros(nUnits, 1);
% for k = 1:nUnits
%     unitCount(k) = sum(strcmp(rateDataSweep.Unit_objectID, uniqueUnits(k)));
%     
% end
% 
% UnitCounts = table(uniqueUnits, unitCount);
% UnitCounts.Properties.VariableNames{1,1} = 'Unit_objectID';
% 
% % Add unitCount info to table; 
% rateDataSweep = join(rateDataSweep, UnitCounts);
% 
% 
% % 2) keep only those neurons that have minimum trial-presence
% isUnitCountAbove = rateDataSweep.unitCount >= nTrialMin;
% rateDataSweep = rateDataSweep(isUnitCountAbove,:);
% nRemainingUniqueUnits = numel(unique(rateDataSweep.Unit_objectID));
% disp(['unique SUs with at leat ', num2str(nTrialMin), ...
%      ' trials: ', num2str(nRemainingUniqueUnits)]);
%  
%  
%  
% end