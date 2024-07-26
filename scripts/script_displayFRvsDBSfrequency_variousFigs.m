% script for creating bootstrapped entropy comparisons for 30-sec stretches
% of DBS trial
%NOTE: still need to put in method for Contact-sweep too

% Cited:
% Moran, A., Stein, E., Tischler, H., Belelovsky, K. & Bar-Gad, 
% I. Dynamic Stereotypic Responses of Basal Ganglia Neurons to Subthalamic
% Nucleus High-Frequency Stimulation in the Parkinsonian Primate. 
% Frontiers in Systems Neuroscience 5, 21–21 (2011).



%% build pipeline-parameter struct "ppar"

clear; 

% PIPELINE PARAMETERS
% full path on local PC where project folder is (don't include subfolders here,
% that's tracked within the appropriate tables)
ppar.projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record'; % string

% full path on local PC where tables are to be loaded from (or saved to)
% ppar.tablePath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string
ppar.tablePath = 'C:\Users\bello043\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string

ppar.intDataRootFolder = 'DataProcessing\intermediateData';

% For getting sub-selection of data from tables
ppar.neuronType = 'SU';
ppar.subjID = 'bothSubj'; % 'Kramer' | 'Uva' | 'bothSubj'
ppar.trialType = 'FrequencySweep'; % 'FrequencySweep' | 'ContactSweep'
% ppar.groups.dbsFrequency.groupFreqs = {10, 20, 30, 50, 100, 130};
% ppar.groups.dbsFrequency.groupLabels = {'hz10', 'hz20', 'hz30', 'hz50', 'hz100', 'hz130'};



% CONSTANTS

% Spike-sort quality:
REFRAC = 1 / 1000; % seconds, neuron refractory period in which spkes shouldn't occur
PERC_REFRACTHRESH = 1.0; % threshold for rejecting a unit with x% spikes within the refractory period


% Filter the data:
NEURON_MIN_TRIALS = 3;
SPARSE_HZ = 2;


% Define time ranges of various peri-DBS events: (in seconds, make sure 2nd
% element > 1st element
  SPKINTERV_PRE = [-25, 0]; % [-30, 0] | [-Inf, 0]
SPKINTERV_DBS_1 = [0, 30];   % [0, 30]
SPKINTERV_DBS_2 = [30, 60];  % [30, 60]
  SPKINTERV_POS = [60, 85]; % [60, 90] | [60, Inf]
  % NOTE: the -31 and 91 is to make it possible to collect bins for -30 and
  % 90 seconds...

BINWIDTH = 1; % seconds

% Entropy calculation:
BINS_PER_DECADE = 15;
ORD_H = 1;
NBOOTS = 10000;

% change significance for categorizing changes in decrH / incrH labels
SIG_PVAL = 0.0001 / 2; % to account for two-tailed distributions


% Final figure:
FIG_POSITION = [14 59 560 420];

t = now;
d = datetime(t, 'ConvertFrom', 'datenum');
YYYYMMDD = num2str(yyyymmdd(d));
YYMMDD = YYYYMMDD(3:end);


% Gaussian filter for smoothing time series
w = gausswin(5);
GAUSSW = w / sum(w);


%% OPTIONS

% gaussFilt the Average FR over time
gaussFiltAvFR = false;

% gaussFilt all the FR data over time
gaussFiltFR = true;



%%

% Get the name of the currently running script:
[scriptDirectoryFullPath, scriptName] = fileparts(mfilename('fullpath'));

% First make sure that this script has a folder within the project folder's
% intermediate data section
scriptIntermediateDataFolder = [ppar.projRootPath, '\', ppar.intDataRootFolder, '\', scriptName];
if ~exist(scriptIntermediateDataFolder, 'dir')
    mkdir(scriptIntermediateDataFolder)
    
end

% Perform the pipeline for every Row in tableRoot
% load root table "tableRoot"
load([ppar.tablePath, '\', 'tableRoot'], 'tableRoot');

% Get subselection of tableRoot for analysis
Tselect = filterTable_Master(tableRoot, scriptName, ppar);

% Check if a given unit has good-enough sort quality to include in
% analysis, remove those that are low quality:
[unitIDs, unitISIs, percInRefrac] = func_assessUnitSpkInRefrac(Tselect, REFRAC, ppar);
isQualityUnit = percInRefrac < PERC_REFRACTHRESH;
tabUnitQ = [table(unitIDs), table(isQualityUnit)];
T = join(Tselect, tabUnitQ, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'unitIDs'); 
T(~T.isQualityUnit,:) = []; % REMOVE THOSE ROWS PERTAINING TO POOR QUALITY UNIT
Tselect = T;
disp(['Removed ', num2str(sum(~isQualityUnit)), ' of ', ...
       num2str(numel(unitIDs)),' units due to low sort quality']);
   

nNex = height(Tselect);

% Initialize columns to add to Tanalyze for final firingRate analysis:
    FRpreCorr = zeros(nNex, 1);
 binFRpreCorr = cell(nNex, 1);
   FRdbsCorr1 = zeros(nNex, 1);
binFRdbsCorr1 = cell(nNex, 1);
   FRdbsCorr2 = zeros(nNex, 1);
binFRdbsCorr2 = cell(nNex, 1);
    FRposCorr = zeros(nNex, 1);
 binFRposCorr = cell(nNex, 1);
 

% % Initialize columns to add to Tanalyze for final Entropy analysis in this script:
% Hpre = zeros(nNex, 1);
% 
%           Hdbs1 = zeros(nNex, 1);
% Hdbs1_preBootAv = zeros(nNex, 1);
%   Hdbs1_EmpPval = zeros(nNex, 1);
% 
%           Hdbs2 = zeros(nNex, 1);
% Hdbs2_preBootAv = zeros(nNex, 1);
%   Hdbs2_EmpPval = zeros(nNex, 1);
% 
% Hpos = zeros(nNex, 1);

% For each nexfile contained within each record-row in Tselect, calculate
% the corrected Firing Rate and entropies for pre, dbs1, dbs2, and post
tic
for iNex = 1:nNex

    %% load nexfile of interest
%     iNex
    iNexObjectID = Tselect.objectID{iNex,1};
           nexfn = Tselect.nexFile{iNex,1};
           nexpn = Tselect.nexFileFolder{iNex,1};
         nexFile = readNexFile([ppar.projRootPath, '\', nexpn, '\', nexfn]);

    % extract spike times and DBS stim times from nexFile struct
    [spkTimes, StimTs] = parseNexFile(nexFile);

    % make sure spkTimes are all of time values that are monotonically
    % increasing:
    spkTimes = sort(spkTimes);



    %% Get PRE-DBS period spike rate, both observed and corrected for artifact
    % blanking

    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'FRpre';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'FRpre_%s_%ss-%ss_fromDbsOnset_%ssBins';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_PRE(1)), ...
                              num2str(SPKINTERV_PRE(2)), ...
                              num2str(BINWIDTH));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
        preDbsInterval = SPKINTERV_PRE; %s, seconds to gather spikes and virtual stims before DBS onset
        % check and make sure that recording has enough pre-DBS time as requested
        % in "preDbsInterval", crop preDBS interval if necessary
        if abs(preDbsInterval(1)) > StimTs.DBS(1)
            preDbsInterval(1) = -StimTs.DBS(1);

        end
        
        FRpre = calcFiringRates(spkTimes, StimTs, preDbsInterval);
        FRpre.intervalLimits = preDbsInterval;
        
        
        % Generate bin-counts of spks/sec for this interval
        binFRpre = getBinRates(spkTimes, StimTs.DBS(1), preDbsInterval, BINWIDTH);
        
       
        save(fullPathFn, 'FRpre', 'binFRpre');
    
    end
    
    FRpreCorr(iNex,1) = FRpre.rateCorrected;
    % correct for blanking-effect
    binFRpreCorr{iNex,1} = binFRpre * (FRpre.totIntervTime / (FRpre.totIntervTime - FRpre.totBlankTime));

            

    %% Get DBS-ON period spike rate, both observed and corrected for artifact
    % blanking for first 30 sec

    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'FRdbs';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'FRdbs_%s_%ss-%ss_fromDbsOnset_%ssBins';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_DBS_1(1)), ...
                              num2str(SPKINTERV_DBS_1(2)), ...
                              num2str(BINWIDTH));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
        FRdbs = calcFiringRates(spkTimes, StimTs, SPKINTERV_DBS_1);
        FRdbs.intervalLimits = SPKINTERV_DBS_1;
        
        % Generate bin-counts of spks/sec for this interval
        binFRdbs = getBinRates(spkTimes, StimTs.DBS(1), SPKINTERV_DBS_1, BINWIDTH);
        
        save(fullPathFn, 'FRdbs', 'binFRdbs');
    
    end
    
    FRdbsCorr1(iNex,1) = FRdbs.rateCorrected;
    % correct for blanking-effect
    binFRdbsCorr1{iNex,1} = binFRdbs * (FRdbs.totIntervTime / (FRdbs.totIntervTime - FRdbs.totBlankTime));

    
    
    %% Get DBS-ON period spike rate, both observed and corrected for artifact
    % blanking for SECOND 30 sec

    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'FRdbs';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'FRdbs_%s_%ss-%ss_fromDbsOnset_%ssBins';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_DBS_2(1)), ...
                              num2str(SPKINTERV_DBS_2(2)), ...
                              num2str(BINWIDTH));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
        FRdbs = calcFiringRates(spkTimes, StimTs, SPKINTERV_DBS_2);
        FRdbs.intervalLimits = SPKINTERV_DBS_2;
        
        % Generate bin-counts of spks/sec for this interval
        binFRdbs = getBinRates(spkTimes, StimTs.DBS(1), SPKINTERV_DBS_2, BINWIDTH);
        
        save(fullPathFn, 'FRdbs', 'binFRdbs');
            
    end
    
    FRdbsCorr2(iNex,1) = FRdbs.rateCorrected;
    % correct for blanking-effect
    binFRdbsCorr2{iNex,1} = binFRdbs * (FRdbs.totIntervTime / (FRdbs.totIntervTime - FRdbs.totBlankTime));

    
    
    %% Get POST-DBS period spike rate, both observed and corrected for artifact
    % blanking

    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'FRpos';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'FRpos_%s_%ss-%ss_fromDbsOnset_%ssBins';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_POS(1)), ...
                              num2str(SPKINTERV_POS(2)), ...
                              num2str(BINWIDTH));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
        posDbsInterval = SPKINTERV_POS; %s, seconds to gather spikes and virtual stims before DBS onset
        % check and make sure that recording has enough pre-DBS time as requested
        % in "preDbsInterval", crop preDBS interval if necessary
        if posDbsInterval(2) > nexFile.tend
            posDbsInterval(2) = nexFile.tend;

        end
        
        FRpos = calcFiringRates(spkTimes, StimTs, posDbsInterval);
        FRpos.intervalLimits = posDbsInterval;
        
        % Generate bin-counts of spks/sec for this interval
        binFRpos = getBinRates(spkTimes, StimTs.DBS(1), posDbsInterval, BINWIDTH);
        
        save(fullPathFn, 'FRpos', 'binFRpos');
    
    end
    
    FRposCorr(iNex,1) = FRpos.rateCorrected;
    % correct for blanking-effect
    binFRposCorr{iNex,1} = binFRpos * (FRpos.totIntervTime / (FRpos.totIntervTime - FRpos.totBlankTime));

    
    
end
toc



%% Collate FR bins over time for each row (put this before final table reduction)
% add to table for analysis

   timeFRall = cell(nNex, 1);
binFRallCorr = cell(nNex, 1);

for iNex = 1:nNex
    % First create time vectors for each section of the trial
     nBinsPre = length(binFRpreCorr{iNex});
    nBinsDbs1 = length(binFRdbsCorr1{iNex});
    nBinsDbs2 = length(binFRdbsCorr2{iNex});
     nBinsPos = length(binFRposCorr{iNex});

    timeFRdbs1 = (1:nBinsDbs1) * BINWIDTH;
    timeFRdbs2 = ((nBinsDbs1 + 1):(nBinsDbs1 + nBinsDbs2)) * BINWIDTH;
     timeFRpos = ((nBinsDbs1 + nBinsDbs2 + 1):(nBinsDbs1 + nBinsDbs2 + nBinsPos)) * BINWIDTH;
     timeFRpre = (-(nBinsPre - 1):1:0) * BINWIDTH;
    
    % Store result in cell array
    timeFRall{iNex} = [timeFRpre, timeFRdbs1, timeFRdbs2, timeFRpos];

    % Second group all bins for each trial in the correct order
    % Store this result cell array as well
    binFRallCorr{iNex} = [binFRpreCorr{iNex}, binFRdbsCorr1{iNex}, ...
                          binFRdbsCorr2{iNex}, binFRposCorr{iNex}];

end
    


%% Prepare data for final visualization and analysis

% Append data of interest to Tselect
Tfinal = [Tselect, ...
          table(FRpreCorr), table(FRdbsCorr1), table(FRdbsCorr2), table(FRposCorr), ...
          table(timeFRall), table(binFRallCorr)];

numel(unique(Tfinal.Unit_objectID)) % display current number units       
Tfinal2 = Tfinal;

% Tanalyze = Tfinal;

% Remove any rows that have a firing rate of less than min firing rate:
% filter out any rows that have too-sparse cells, based on 30-second
% portions
isPreRateBelow = Tfinal2{:, 'FRpreCorr'} < SPARSE_HZ;
% isDbs1RateBelow = Tfinal2{:, 'FRdbsCorr1'} < SPARSE_HZ;
isDbs2RateBelow = Tfinal2{:, 'FRdbsCorr2'} < SPARSE_HZ;
% isRemove = isPreRateBelow | isDbs1RateBelow | isDbs2RateBelow;
isRemove = isPreRateBelow & isDbs2RateBelow;

Tfinal2(isRemove, :) = [];
numel(unique(Tfinal2.Unit_objectID))


% Remove any rows that have neurons that were not recorded with enough time
% as specified by tStard and tEnd (where t = 0 is onset of DBS). Note this
% is not the same as excluding neurons that were recorded to have a FR of
% 0.
tStart = SPKINTERV_PRE(1) + 1;
tEnd = SPKINTERV_POS(2) - 1;
timeFRstandard = tStart:BINWIDTH:tEnd;
nTimeBins = numel(timeFRstandard);
nRows = height(Tfinal2);
hasEnoughTimeBins = false(nRows, 1);
T = Tfinal2;
for iRow = 1:nRows
    % check if time bins for this row contains the requisite start and end
    % times
    timeBins = T.timeFRall{iRow,1};
    if (timeBins(1) <= tStart) && (timeBins(end) >= tEnd)
        hasEnoughTimeBins(iRow) = true;
        
    end
        
end
T(~hasEnoughTimeBins,:) = [];
disp(['Reomved ', num2str(sum(~hasEnoughTimeBins)), ' of ', ...
       num2str(nRows),' table rows due to not enough time recorded']);
Tfinal2 = T;


% Remove any neurons that have less than minimum trials 
Tfinal3 = filter_neuronMinimumTrials(Tfinal2, NEURON_MIN_TRIALS);
numel(unique(Tfinal3.Unit_objectID))

Tanalyze = Tfinal3;


% Add in average FR line from -30 to 60 seconds, for each dbsFrequency


% Reassign tStart to the shortest pre-DBS start time in the data, so that
% all data will have equal length




%% Display FR bins over time, group axes by DBS frequency

Tdisp = Tanalyze;
freqs = [10, 20, 30, 50, 100, 130];
f1 = figure;

f1.Position = [1939 80 1833 903];

% if option is checked, gauss-filter all time-series data
if gaussFiltFR % smooth all the data itself
    Tdisp = tab_gaussFiltDataCell(Tdisp, 'binFRallCorr', GAUSSW);

end

% plot bins over time, add in vertical lines showing beginning and end of DBS

% choose subselection of data by 10 Hz
for dbsCond = 1:6
    ax(dbsCond) = plotFRoverTime(dbsCond, Tdisp, freqs);
    
end

set(ax, 'XLim', [-60, 120])
set(ax, 'YLim', [0, 40])
% set(ax, 'YLim', [0, 180])

for dbsCond = 1:6, plotDBSlines(ax(dbsCond), 0, 60); end
ylabel(ax(1), 'FR (spk/sec)')
ylabel(ax(4), 'FR (spk/sec)')
xlabel(ax(4), 'Time (seconds), 0 is DBS onset')




for dbsCond = 1:6
    % gather all table rows for current dbsCond:
    isFr = Tdisp.dbsFrequency == freqs(dbsCond);
    T = Tdisp(isFr,:);
    nTrials = height(T);
    
    % Fill matrix with values for each FR time series
    FRoverTimeAll = zeros(nTrials, nTimeBins);
    for iTr = 1:nTrials
        % for this trial, get data indices that correspond to -30 to 60 sec
        iTr_time = T.timeFRall{iTr};
        timeIdx = (iTr_time >= tStart) & (iTr_time <= tEnd);
        iTr_FRoverTime = T.binFRallCorr{iTr};
        FRoverTimeAll(iTr,:) = iTr_FRoverTime(timeIdx);
        
    end
    

    
    % get average FR from -30 to 60 seconds
    FRoverTimeAll_av = mean(FRoverTimeAll, 1);
    FRoverTimeAll_stdv = std(FRoverTimeAll, [], 1);
    FRoverTimeAll_sem = std(FRoverTimeAll, [], 1) ./ sqrt(size(FRoverTimeAll, 1));
    if gaussFiltAvFR % smooth the average of the data
          FRoverTimeAll_av = filtfilt(GAUSSW, 1, FRoverTimeAll_av); % smooth it
        FRoverTimeAll_stdv = filtfilt(GAUSSW, 1, FRoverTimeAll_stdv); % smooth it
         FRoverTimeAll_sem = filtfilt(GAUSSW, 1, FRoverTimeAll_sem); % smooth it
        
    end
    
    
    % plot average and added SEM lines in proper axes
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av, ...
         'LineWidth', 2, 'Color', 'k');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av + FRoverTimeAll_sem, ...
         'LineWidth', 1, 'Color', 'k');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av - FRoverTimeAll_sem, ...
         'LineWidth', 1, 'Color', 'k');    
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av + FRoverTimeAll_stdv, ...
         'LineWidth', 1, 'Color', 'r');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av - FRoverTimeAll_stdv, ...
         'LineWidth', 1, 'Color', 'r');     
   
    
end



%% Display FR bins over time normalized to PreDBS % change, group axes by DBS frequency

% Initial step normalizing bins to average 
T = Tanalyze;
nTrials = height(T);
binFRallCorr_norm = cell(nTrials, 1);

for iTr = 1:nTrials
    bins = T.binFRallCorr{iTr};
    preAv = T.FRpreCorr(iTr);
    binFRallCorr_norm{iTr} = 100 * (bins - preAv) ./ preAv; 
    
end
T = [T, table(binFRallCorr_norm)];
Tanalyze = T;


% display individual lines of DBS FR %change normalized by pre-DBS average
Tdisp = Tanalyze;
freqs = [10, 20, 30, 50, 100, 130];
f2 = figure;

f2.Position = [1939 80 1833 903];

% if option is checked, gauss-filter all time-series data
if gaussFiltFR % smooth all the data itself
    Tdisp = tab_gaussFiltDataCell(Tdisp, 'binFRallCorr_norm', GAUSSW);

end

% plot bins over time, add in vertical lines showing beginning and end of DBS

% choose subselection of data by 10 Hz
for dbsCond = 1:6
    ax(dbsCond) = plotFRoverTime_norm(dbsCond, Tdisp, freqs);
    
end

% Also draw an average line of all displayed lines, with SEM


% set all axes limits to be equal
set(ax, 'XLim', [-60, 120])
set(ax, 'YLim', [-100, 250])
% set(ax, 'YLim', [-100, 1500])

% draw red lines to show dbs onset and offset, and a horizontal line for 0%
% FR change
for dbsCond = 1:6
    plotDBSlines(ax(dbsCond), 0, 60);
    plotHorizLines(ax(dbsCond), 0);
    
end
ylabel(ax(1), '%\DeltaFR from pre-DBS av.')
ylabel(ax(4), '%\DeltaFR from pre-DBS av.')
xlabel(ax(4), 'Time (seconds), 0 is DBS onset')


% Add in average FR normalized-change line from -30 to 60 seconds, for each
% dbsFrequency
% tStart = -30;
% tEnd = 90;
% timeFRstandard = tStart:BINWIDTH:tEnd;
% nTimeBins = numel(timeFRstandard);
for dbsCond = 1:6
    % gather all table rows for current dbsCond:
    isFr = Tdisp.dbsFrequency == freqs(dbsCond);
    T = Tdisp(isFr,:);
    nTrials = height(T);
    
    % Fill matrix with values for each FR time series
    FRoverTimeAll = zeros(nTrials, nTimeBins);
    for iTr = 1:nTrials
        % for this trial, get data indices that correspond to -30 to 60 sec
        iTr_time = T.timeFRall{iTr};
        timeIdx = (iTr_time >= tStart) & (iTr_time <= tEnd);
        iTr_FRoverTime = T.binFRallCorr_norm{iTr};
        FRoverTimeAll(iTr,:) = iTr_FRoverTime(timeIdx);
        
    end
       
    % get average FR from -30 to 60 seconds
    FRoverTimeAll_av = mean(FRoverTimeAll, 1);
    FRoverTimeAll_stdv = std(FRoverTimeAll, [], 1);
    FRoverTimeAll_sem = std(FRoverTimeAll, [], 1) ./ sqrt(size(FRoverTimeAll, 1));
    if gaussFiltAvFR
          FRoverTimeAll_av = filtfilt(GAUSSW, 1, FRoverTimeAll_av); % smooth it
        FRoverTimeAll_stdv = filtfilt(GAUSSW, 1, FRoverTimeAll_stdv); % smooth it
         FRoverTimeAll_sem = filtfilt(GAUSSW, 1, FRoverTimeAll_sem); % smooth it
        
    end
    
    % plot average and added SEM lines in proper axes
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av, ...
         'LineWidth', 2, 'Color', 'k');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av + FRoverTimeAll_sem, ...
         'LineWidth', 1, 'Color', 'k');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av - FRoverTimeAll_sem, ...
         'LineWidth', 1, 'Color', 'k');    
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av + FRoverTimeAll_stdv, ...
         'LineWidth', 1, 'Color', 'r');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av - FRoverTimeAll_stdv, ...
         'LineWidth', 1, 'Color', 'r');     
    
end



%% Display FR bins over time normalized to PreDBS FR index, group axes by DBS frequency

% Initial step normalizing bins to average 
T = Tanalyze;
nTrials = height(T);
binFRallCorr_fridx = cell(nTrials, 1);

for iTr = 1:nTrials
    bins = T.binFRallCorr{iTr};
    preAv = T.FRpreCorr(iTr);
    binFRallCorr_fridx{iTr} = (bins - preAv) ./ (bins + preAv); 
    
end
T = [T, table(binFRallCorr_fridx)];
Tanalyze = T;


% display
Tdisp = Tanalyze;
freqs = [10, 20, 30, 50, 100, 130];
f2 = figure;

f2.Position = [1939 80 1833 903];

% if option is checked, gauss-filter all time-series data
if gaussFiltFR % smooth all the data itself
    Tdisp = tab_gaussFiltDataCell(Tdisp, 'binFRallCorr_norm', GAUSSW);

end

% plot bins over time, add in vertical lines showing beginning and end of DBS

% choose subselection of data by 10 Hz
for dbsCond = 1:6
    ax(dbsCond) = plotFRoverTime_fridx(dbsCond, Tdisp, freqs);
    
end

% set all axes limits to be equal
set(ax, 'XLim', [-60, 120])
set(ax, 'YLim', [-1, 1])

% draw red lines to show dbs onset and offset, and a horizontal line for 0%
% FR change
for dbsCond = 1:6
    plotDBSlines(ax(dbsCond), 0, 60);
    plotHorizLines(ax(dbsCond), 0);
    
end
ylabel(ax(1), '\DeltaFR index')
ylabel(ax(4), '\DeltaFR index')
xlabel(ax(4), 'Time (seconds), 0 is DBS onset')


% Add in average deltaFR index line from -30 to 60 seconds, for each
% dbsFrequency
% tStart = -30;
% tEnd = 90;
% timeFRstandard = tStart:BINWIDTH:tEnd;
% nTimeBins = numel(timeFRstandard);
for dbsCond = 1:6
    % gather all table rows for current dbsCond:
    isFr = Tdisp.dbsFrequency == freqs(dbsCond);
    T = Tdisp(isFr,:);
    nTrials = height(T);
    
    % Fill matrix with values for each FR time series
    FRoverTimeAll = zeros(nTrials, nTimeBins);
    for iTr = 1:nTrials
        % for this trial, get data indices that correspond to -30 to 60 sec
        iTr_time = T.timeFRall{iTr};
        timeIdx = (iTr_time >= tStart) & (iTr_time <= tEnd);
        iTr_FRoverTime = T.binFRallCorr_fridx{iTr};
        FRoverTimeAll(iTr,:) = iTr_FRoverTime(timeIdx);
        
    end
    
    
    % get average FR from -30 to 60 seconds
    FRoverTimeAll_av = mean(FRoverTimeAll, 1);
    FRoverTimeAll_stdv = std(FRoverTimeAll, [], 1);
    FRoverTimeAll_sem = std(FRoverTimeAll, [], 1) ./ sqrt(size(FRoverTimeAll, 1));
    if gaussFiltAvFR
          FRoverTimeAll_av = filtfilt(GAUSSW, 1, FRoverTimeAll_av); % smooth it
        FRoverTimeAll_stdv = filtfilt(GAUSSW, 1, FRoverTimeAll_stdv); % smooth it
         FRoverTimeAll_sem = filtfilt(GAUSSW, 1, FRoverTimeAll_sem); % smooth it
        
    end
    
    % plot average and added SEM lines in proper axes
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av, ...
         'LineWidth', 2, 'Color', 'k');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av + FRoverTimeAll_sem, ...
         'LineWidth', 1, 'Color', 'k');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av - FRoverTimeAll_sem, ...
         'LineWidth', 1, 'Color', 'k');    
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av + FRoverTimeAll_stdv, ...
         'LineWidth', 1, 'Color', 'r');
    plot(ax(dbsCond), timeFRstandard, FRoverTimeAll_av - FRoverTimeAll_stdv, ...
         'LineWidth', 1, 'Color', 'r');     
    
end



%% Show Average FR over time during DBS normalized to preDBS, with SEM to show variability
% 
% % Initial step normalizing bins to average 
% T = Tanalyze;
% nTrials = height(T);
% binFRallCorr_norm = cell(nTrials, 1);
% 
% for iTr = 1:nTrials
%     bins = T.binFRallCorr{iTr};
%     preAv = T.FRpreCorr(iTr);
%     binFRallCorr_norm{iTr} = 100 * (bins - preAv) ./ preAv; 
%     
% end
% T = [T, table(binFRallCorr_norm)];
% Tanalyze = T;
% 

% 
% 
% %% Show barplots of raw av firing rates during first 30 seconds and last 30 seconds of DBS
% T = Tanalyze; 
% % gather average FR values 
% 
% FR1 = T.FRdbsCorr1;
% FR2 = T.FRdbsCorr2;
% grps = T.dbsFrequency;
% 
% figure; 
% 
% ax1 = subplot(1,2,1);
% bar(FR1, grps);
% xlabel('DBS frequency');
% tit1 = title([ppar.subjID, ': FR1 vs. Frequency']);
% ax1.YGrid = 'on';
% ax1.YMinorGrid = 'on';
% 
% ax2 = subplot(1,2,2);
% boxplot(FR2, grps);
% xlabel('DBS frequency');
% tit2 = title([ppar.subjID, ': FR2 change by Frequency']);
% ax2.YGrid = 'on';
% ax2.YMinorGrid = 'on';
% 
% YLimMax = max(ax1.YLim(2), ax2.YLim(2));
% YLimMin = min(ax1.YLim(1), ax2.YLim(1));
% ax1.YLim = [YLimMin, YLimMax];
% ax2.YLim = [YLimMin, YLimMax];
% 
% 
% 
% 
% 
% 
% 
% %% For first 30 sec and last 30 sec, see firing rate changes by DBS freq
% 
% T = Tanalyze; 
% % gather deltaFR values 
% 
% deltaFR1 = T.FRdbsCorr1 - T.FRpreCorr;
% deltaFR2 = T.FRdbsCorr2 - T.FRpreCorr;
% grps = T.dbsFrequency;
% 
% figure; 
% 
% ax1 = subplot(1,2,1);
% boxplot(deltaFR1, grps);
% xlabel('DBS frequency');
% tit1 = title([ppar.subjID, ': FR1 change by Frequency']);
% ax1.YGrid = 'on';
% ax1.YMinorGrid = 'on';
% 
% ax2 = subplot(1,2,2);
% boxplot(deltaFR2, grps);
% xlabel('DBS frequency');
% tit2 = title([ppar.subjID, ': FR2 change by Frequency']);
% ax2.YGrid = 'on';
% ax2.YMinorGrid = 'on';
% 
% YLimMax = max(ax1.YLim(2), ax2.YLim(2));
% YLimMin = min(ax1.YLim(1), ax2.YLim(1));
% ax1.YLim = [YLimMin, YLimMax];
% ax2.YLim = [YLimMin, YLimMax];
% 
% 
% 
% %% Display progression of firing rate over time, from pre thru to post
% 
% 
% 
% 
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1, 3, 1);
% boxplot(diffFR); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(diffFR(:,1:3)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:3;
% ax2.XTickLabel = {'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot(diffFR(:,4:6)'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:3;
% ax3.XTickLabel = {'HFS1', 'HFS2', 'postHFS'};
% % 
% ax1.YLim = [-40, 100];
% ax2.YLim = [-40, 100];
% ax3.YLim = [-40, 100];

% 
% 
% 
% %% Hdbs1: For each trial, categorize as Entorpy incr, decr, and significance
% % based on bootstrapped distributions of pre-DBS entropy
% 
% dbsPortionStr = 'Hdbs1';
% 
% % % change significance
% % SIG_PVAL = 0.05 / 2; % to account for two-tailed distributions
% 
% Tdisp = Tanalyze; % just to make typing easier...
% nRows = height(Tdisp);
% 
%      wasDecrH = false(nRows, 1);
%      wasIncrH = false(nRows, 1);
% signifChangeH = false(nRows, 1);
% 
% for iRow = 1:nRows
%     if Tdisp.Hdbs1(iRow) <= Tdisp.Hdbs1_preBootAv(iRow) % entorpy decrease
%         wasDecrH(iRow) = true; 
%         
%     else % entropy increase
%         wasIncrH(iRow) = true;
%         
%     end
%     
%     if abs(Tdisp.Hdbs1_EmpPval(iRow)) < SIG_PVAL
%         signifChangeH(iRow) = true;
%         
%     end
%     
% end
% 
% Tdisp = [Tdisp, table(wasDecrH), table(wasIncrH), table(signifChangeH)];
% 
% 
% % Hdbs1: Display the group proportions according to frequency
% 
% % pre-allocate matrix for counts with each frequency being a row; columns:
% % wasDecrH, wasIncrH, noCh
% freqs = [10, 20, 30, 50, 100, 130];
% nFreqs = numel(freqs);
% catCounts = zeros(nFreqs, 3);
% 
% % Gather all counts
% for iFr = 1:nFreqs
%     % get subselection of table rows pertaining to frequency "iFr"
%     T_byFreq = Tdisp(Tdisp.dbsFrequency == freqs(iFr),:);
%     
%     % categorize each trial into "incrH", "decrH", "noCh" (in that order), 
%     % keep count
%     catCounts(iFr,:) = getCategoryCounts_table(T_byFreq);
%   
% end
%  
% % transfer to percentages
% tot = sum(catCounts, 2);
% catPercs = catCounts ./ tot;
% 
% % display results
% figure; ax = axes;
% bar(catPercs * 100)
% ax.XTickLabel = freqs;
% grid on
% xlabel('DBS frequency (Hz)')
% ylabel('% recorded cells')
% legend('H increase', 'H decrease', 'noChange', ...
%       'Location', 'northeastoutside') ;
% tit = title([ppar.subjID, ':  population percentage by ', dbsPortionStr,' Entropy change']);
% ax.YLim = [0, 100];
% 
% 
% 
% %% Hdbs2: For each trial, categorize as Entorpy incr, decr, and significance
% % based on bootstrapped distributions of pre-DBS entropy
% 
% dbsPortionStr = 'Hdbs2';
% 
% % % change significance
% % SIG_PVAL = 0.05 / 2; % to account for two-tailed distributions
% 
% Tdisp = Tanalyze; % just to make typing easier...
% nRows = height(Tdisp);
% 
%      wasDecrH = false(nRows, 1);
%      wasIncrH = false(nRows, 1);
% signifChangeH = false(nRows, 1);
% 
% for iRow = 1:nRows
%     if Tdisp.Hdbs2(iRow) <= Tdisp.Hdbs2_preBootAv(iRow) % entorpy decrease
%         wasDecrH(iRow) = true; 
%         
%     else % entropy increase
%         wasIncrH(iRow) = true;
%         
%     end
%     
%     if abs(Tdisp.Hdbs2_EmpPval(iRow)) < SIG_PVAL
%         signifChangeH(iRow) = true;
%         
%     end
%     
% end
% 
% Tdisp = [Tdisp, table(wasDecrH), table(wasIncrH), table(signifChangeH)];
% 
% % Hdbs2: Display the group proportions according to frequency
% % pre-allocate matrix for counts with each frequency being a row; columns:
% % wasDecrH, wasIncrH, noCh
% freqs = [10, 20, 30, 50, 100, 130];
% nFreqs = numel(freqs);
% catCounts = zeros(nFreqs, 3);
% 
% % Gather all counts
% for iFr = 1:nFreqs
%     % get subselection of table rows pertaining to frequency "iFr"
%     T_byFreq = Tdisp(Tdisp.dbsFrequency == freqs(iFr),:);
%     
%     % categorize each trial into "incrH", "decrH", "noCh" (in that order), 
%     % keep count
%     catCounts(iFr,:) = getCategoryCounts_table(T_byFreq);
%   
% end
%  
% % transfer to percentages
% tot = sum(catCounts, 2);
% catPercs = catCounts ./ tot;
% 
% % display results
% figure; ax = axes;
% bar(catPercs * 100)
% ax.XTickLabel = freqs;
% grid on
% xlabel('DBS frequency (Hz)')
% ylabel('% recorded cells')
% lgd = legend('H increase', 'H decrease', 'noChange', ...
%       'Location', 'northeastoutside') ;
% tit = title([ppar.subjID, ':  population percentage by ', dbsPortionStr,' Entropy change']);
% ax.YLim = [0, 100];



%%

% %% Display boxplot and individual neuron progression for LFS adn HFS groups
% % Bits/Spike
% 
% % 'FRpreCorr'    | 'FRdbsCorr1'    | 'FRdbsCorr2'    | 'FRposCorr'
% % 'Hpre'         | 'Hdbs1'         | 'Hdbs2'         | 'Hpos'
% % 'Hpre_BitpSec' | 'Hdbs1_BitpSec' | 'Hdbs2_BitpSec' | 'Hpos_BitpSec'
% % 'FRdbs1_idx'       | 'FRdbs2_idx'        | 'FRpos_idx'
% % 'Hdbs1_BitpSPKperc | 'Hdbs2_BitpSPKperc' | 'Hpos_BitpSPKperc'
% % 'Hdbs1_BitpSECperc | 'Hdbs2_BitpSECperc' | 'Hpos_BitpSECperc'
% 
% uniqueNeus = unique(Tfinal2.Unit_objectID);
% nNeus = numel(uniqueNeus);
% % pre-fill an array of data values
% baseline = zeros(nNeus, 1);
% LFS1 = zeros(nNeus, 1);
% LFS2 = zeros(nNeus, 1);
% postLFS = zeros(nNeus, 1);
% 
% HFS1 = zeros(nNeus, 1);
% HFS2 = zeros(nNeus, 1);
% postHFS = zeros(nNeus, 1);
% 
% freqs = zeros(nNeus, 2);
% grpLFS = [10, 20, 30];
% grpHFS = [50, 100, 130];
% for iNeu = 1:nNeus
%     clear T_iNeu freqLowest freqHighest
%     iNeuStr = uniqueNeus{iNeu};
%     % get sub-selection of table for current neuron
%     idx_iNeu = strcmp(iNeuStr, Tfinal2.Unit_objectID);
%     T_iNeu = Tfinal2(idx_iNeu,:);
%     
%     % First retreive the earliest pre-DBS period to use as Baseline
%     T_iNeu = sortrows(T_iNeu, 'blockNum');
%     baseline(iNeu) = T_iNeu{1, 'Hpre'};
%      
%     T_iNeu = sortrows(T_iNeu, 'dbsFrequency');
%     
%     % determine lowest frequency present within LFS group & Retrieve the 
%     % dataVariable pertaining to lowest LFS trial
%     freqLowest = T_iNeu.dbsFrequency(1);
%     if any(grpLFS == freqLowest)
%         freqs(iNeu,1) = freqLowest;
%         LFS1(iNeu) = T_iNeu{1,'Hdbs1'}; 
%         LFS2(iNeu) = T_iNeu{1,'Hdbs2'};
%         postLFS(iNeu) = T_iNeu{1,'Hpos'};
%         
%     else % if no LFS group frequency trial exists for this neuron
%         freqs(iNeu,1) = NaN;
%         LFS1(iNeu) = NaN;
%         LFS1(iNeu) = NaN;
%         postLFS(iNeu) = NaN;
%         
%     end
%     % determine highest frequency present within HFS group & Retrieve the 
%     % dataVariable pertaining to highest HFS trial
%     freqHighest = T_iNeu.dbsFrequency(end);
%     if any(grpHFS == freqHighest)
%         freqs(iNeu,2) = freqHighest;
%         HFS1(iNeu) = T_iNeu{end,'Hdbs1'};
%         HFS2(iNeu) = T_iNeu{end,'Hdbs2'};
%         postHFS(iNeu) = T_iNeu{end,'Hpos'};
%         
%     else % if no LFS group frequency trial exists for this neuron
%         freqs(iNeu,2) = NaN;
%         HFS1(iNeu) = NaN;     
%         HFS2(iNeu) = NaN;     
%         postHFS(iNeu) = NaN;
%         
%     end
% end
% 
% dataVar = [baseline, LFS1, LFS2, postLFS, HFS1, HFS2, postHFS];
% % Remove any pairs that don't have at least one LFS and one HFS
% remLFS = isnan(LFS1) | isnan(LFS2) | isnan(postLFS);
% remHFS = isnan(HFS1) | isnan(HFS2) | isnan(postHFS);
% freqs(remLFS | remHFS,:) = [];
% dataVar(remLFS | remHFS,:) = [];
% 
% nNeusUpdated = size(dataVar, 1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1657 376];
% ax1 = subplot(1, 3, 1);
% boxplot(dataVar); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'Baseline', 'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(dataVar(:,1:4)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:4;
% ax2.XTickLabel = {'Baseline', 'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot([dataVar(:,1), dataVar(:,5:7)]'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:4;
% ax3.XTickLabel = {'Baseline', 'HFS1', 'HFS2', 'postHFS'};
% 
% ax1.YLim = [0.5, 5.5];
% ax2.YLim = [0.5, 5.5];
% ax3.YLim = [0.5, 5.5];
% 
% 
% % Plot Change in Entropy
% % Add in extra step to get differences in Entropy referenced to Baseline
% diffEntropy = zeros(size(dataVar,1), size(dataVar,2) - 1);
% diffEntropy(:,1) = dataVar(:,2) - dataVar(:,1);
% diffEntropy(:,2) = dataVar(:,3) - dataVar(:,1);
% diffEntropy(:,3) = dataVar(:,4) - dataVar(:,1);
% diffEntropy(:,4) = dataVar(:,5) - dataVar(:,1);
% diffEntropy(:,5) = dataVar(:,6) - dataVar(:,1);
% diffEntropy(:,6) = dataVar(:,7) - dataVar(:,1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1, 3, 1);
% boxplot(diffEntropy); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(diffEntropy(:,1:3)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:3;
% ax2.XTickLabel = {'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot(diffEntropy(:,4:6)'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:3;
% ax3.XTickLabel = {'HFS1', 'HFS2', 'postHFS'};
% 
% ax1.YLim = [-4, 1.5];
% ax2.YLim = [-4, 1.5];
% ax3.YLim = [-4, 1.5];
% 
% 
% % Plot %change in Entropy
% % Add in extra step to get differences in Entropy referenced to Baseline
% pdiffEntropy = zeros(size(dataVar,1), size(dataVar,2) - 1);
% pdiffEntropy(:,1) = 100 * (dataVar(:,2) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,2) = 100 * (dataVar(:,3) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,3) = 100 * (dataVar(:,4) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,4) = 100 * (dataVar(:,5) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,5) = 100 * (dataVar(:,6) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,6) = 100 * (dataVar(:,7) - dataVar(:,1)) ./ dataVar(:,1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1, 3, 1);
% boxplot(pdiffEntropy); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(pdiffEntropy(:,1:3)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:3;
% ax2.XTickLabel = {'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot(pdiffEntropy(:,4:6)'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:3;
% ax3.XTickLabel = {'HFS1', 'HFS2', 'postHFS'};
% 
% ax1.YLim = [-80, 50];
% ax2.YLim = [-80, 50];
% ax3.YLim = [-80, 50];
% 
% 
% 
% %% Display boxplot and individual neuron progression for LFS adn HFS groups
% % Bits/Spike
% 
% uniqueNeus = unique(Tfinal2.Unit_objectID);
% nNeus = numel(uniqueNeus);
% % pre-fill an array of data values
% baseline = zeros(nNeus, 1);
% LFS1 = zeros(nNeus, 1);
% LFS2 = zeros(nNeus, 1);
% postLFS = zeros(nNeus, 1);
% 
% HFS1 = zeros(nNeus, 1);
% HFS2 = zeros(nNeus, 1);
% postHFS = zeros(nNeus, 1);
% 
% freqs = zeros(nNeus, 2);
% grpLFS = [10, 20, 30];
% grpHFS = [50, 100, 130];
% for iNeu = 1:nNeus
%     clear T_iNeu freqLowest freqHighest
%     iNeuStr = uniqueNeus{iNeu};
%     % get sub-selection of table for current neuron
%     idx_iNeu = strcmp(iNeuStr, Tfinal2.Unit_objectID);
%     T_iNeu = Tfinal2(idx_iNeu,:);
%     
%     % First retreive the earliest pre-DBS period to use as Baseline
%     T_iNeu = sortrows(T_iNeu, 'blockNum');
%     baseline(iNeu) = T_iNeu{1, 'FRpreCorr'};
%      
%     T_iNeu = sortrows(T_iNeu, 'dbsFrequency');
%     
%     % determine lowest frequency present within LFS group & Retrieve the 
%     % dataVariable pertaining to lowest LFS trial
%     freqLowest = T_iNeu.dbsFrequency(1);
%     if any(grpLFS == freqLowest)
%         freqs(iNeu,1) = freqLowest;
%         LFS1(iNeu) = T_iNeu{1,'FRdbsCorr1'}; 
%         LFS2(iNeu) = T_iNeu{1,'FRdbsCorr2'};
%         postLFS(iNeu) = T_iNeu{1,'FRposCorr'};
%         
%     else % if no LFS group frequency trial exists for this neuron
%         freqs(iNeu,1) = NaN;
%         LFS1(iNeu) = NaN;
%         LFS1(iNeu) = NaN;
%         postLFS(iNeu) = NaN;
%         
%     end
%     % determine highest frequency present within HFS group & Retrieve the 
%     % dataVariable pertaining to highest HFS trial
%     freqHighest = T_iNeu.dbsFrequency(end);
%     if any(grpHFS == freqHighest)
%         freqs(iNeu,2) = freqHighest;
%         HFS1(iNeu) = T_iNeu{end,'FRdbsCorr1'};
%         HFS2(iNeu) = T_iNeu{end,'FRdbsCorr2'};
%         postHFS(iNeu) = T_iNeu{end,'FRposCorr'};
%         
%     else % if no LFS group frequency trial exists for this neuron
%         freqs(iNeu,2) = NaN;
%         HFS1(iNeu) = NaN;     
%         HFS2(iNeu) = NaN;     
%         postHFS(iNeu) = NaN;
%         
%     end
% end
% 
% dataVar = [baseline, LFS1, LFS2, postLFS, HFS1, HFS2, postHFS];
% % Remove any pairs that don't have at least one LFS and one HFS
% remLFS = isnan(LFS1) | isnan(LFS2) | isnan(postLFS);
% remHFS = isnan(HFS1) | isnan(HFS2) | isnan(postHFS);
% freqs(remLFS | remHFS,:) = [];
% dataVar(remLFS | remHFS,:) = [];
% 
% nNeusUpdated = size(dataVar, 1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1657 376];
% ax1 = subplot(1, 3, 1);
% boxplot(dataVar); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'Baseline', 'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(dataVar(:,1:4)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:4;
% ax2.XTickLabel = {'Baseline', 'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot([dataVar(:,1), dataVar(:,5:7)]'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:4;
% ax3.XTickLabel = {'Baseline', 'HFS1', 'HFS2', 'postHFS'};
% 
% ax1.YLim = [0, 110];
% ax2.YLim = [0, 110];
% ax3.YLim = [0, 110];
% 
% 
% % Plot Change in FR
% % Add in extra step to get differences in Entropy referenced to Baseline
% diffFR = zeros(size(dataVar,1), size(dataVar,2) - 1);
% diffFR(:,1) = dataVar(:,2) - dataVar(:,1);
% diffFR(:,2) = dataVar(:,3) - dataVar(:,1);
% diffFR(:,3) = dataVar(:,4) - dataVar(:,1);
% diffFR(:,4) = dataVar(:,5) - dataVar(:,1);
% diffFR(:,5) = dataVar(:,6) - dataVar(:,1);
% diffFR(:,6) = dataVar(:,7) - dataVar(:,1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1, 3, 1);
% boxplot(diffFR); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(diffFR(:,1:3)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:3;
% ax2.XTickLabel = {'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot(diffFR(:,4:6)'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:3;
% ax3.XTickLabel = {'HFS1', 'HFS2', 'postHFS'};
% % 
% ax1.YLim = [-40, 100];
% ax2.YLim = [-40, 100];
% ax3.YLim = [-40, 100];
% 
% 
% % Plot %change in FR
% % Add in extra step to get differences in Entropy referenced to Baseline
% pdiffFR = zeros(size(dataVar,1), size(dataVar,2) - 1);
% pdiffFR(:,1) = 100 * (dataVar(:,2) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffFR(:,2) = 100 * (dataVar(:,3) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffFR(:,3) = 100 * (dataVar(:,4) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffFR(:,4) = 100 * (dataVar(:,5) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffFR(:,5) = 100 * (dataVar(:,6) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffFR(:,6) = 100 * (dataVar(:,7) - dataVar(:,1)) ./ dataVar(:,1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1, 3, 1);
% boxplot(pdiffFR); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(pdiffFR(:,1:3)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:3;
% ax2.XTickLabel = {'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot(pdiffFR(:,4:6)'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:3;
% ax3.XTickLabel = {'HFS1', 'HFS2', 'postHFS'};
% % 
% ax1.YLim = [-100, 800];
% ax2.YLim = [-100, 800];
% ax3.YLim = [-100, 800];
% 
% 
% % Plot FR change-index
% % Add in extra step to get differences in Entropy referenced to Baseline
% idxFR = zeros(size(dataVar,1), size(dataVar,2) - 1);
% idxFR(:,1) = calcFRindex(dataVar(:,1), dataVar(:,2));
% idxFR(:,2) = calcFRindex(dataVar(:,1), dataVar(:,3));
% idxFR(:,3) = calcFRindex(dataVar(:,1), dataVar(:,4));
% idxFR(:,4) = calcFRindex(dataVar(:,1), dataVar(:,5));
% idxFR(:,5) = calcFRindex(dataVar(:,1), dataVar(:,6));
% idxFR(:,6) = calcFRindex(dataVar(:,1), dataVar(:,7));
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1, 3, 1);
% boxplot(idxFR); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(idxFR(:,1:3)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:3;
% ax2.XTickLabel = {'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot(idxFR(:,4:6)'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:3;
% ax3.XTickLabel = {'HFS1', 'HFS2', 'postHFS'};
% % 
% ax1.YLim = [-1, 1];
% ax2.YLim = [-1, 1];
% ax3.YLim = [-1, 1];
% 
% 
% 
% %% Display boxplot and individual neuron progression for LFS adn HFS groups
% % Bits/Spike
% 
% % 'FRpreCorr'    | 'FRdbsCorr1'    | 'FRdbsCorr2'    | 'FRposCorr'
% % 'Hpre'         | 'Hdbs1'         | 'Hdbs2'         | 'Hpos'
% % 'Hpre_BitpSec' | 'Hdbs1_BitpSec' | 'Hdbs2_BitpSec' | 'Hpos_BitpSec'
% % 'FRdbs1_idx'       | 'FRdbs2_idx'        | 'FRpos_idx'
% % 'Hdbs1_BitpSPKperc | 'Hdbs2_BitpSPKperc' | 'Hpos_BitpSPKperc'
% % 'Hdbs1_BitpSECperc | 'Hdbs2_BitpSECperc' | 'Hpos_BitpSECperc'
% 
% uniqueNeus = unique(Tfinal2.Unit_objectID);
% nNeus = numel(uniqueNeus);
% % pre-fill an array of data values
% baseline = zeros(nNeus, 1);
% LFS1 = zeros(nNeus, 1);
% LFS2 = zeros(nNeus, 1);
% postLFS = zeros(nNeus, 1);
% 
% HFS1 = zeros(nNeus, 1);
% HFS2 = zeros(nNeus, 1);
% postHFS = zeros(nNeus, 1);
% 
% freqs = zeros(nNeus, 2);
% grpLFS = [10, 20, 30];
% grpHFS = [50, 100, 130];
% for iNeu = 1:nNeus
%     clear T_iNeu freqLowest freqHighest
%     iNeuStr = uniqueNeus{iNeu};
%     % get sub-selection of table for current neuron
%     idx_iNeu = strcmp(iNeuStr, Tfinal2.Unit_objectID);
%     T_iNeu = Tfinal2(idx_iNeu,:);
%     
%     % First retreive the earliest pre-DBS period to use as Baseline
%     T_iNeu = sortrows(T_iNeu, 'blockNum');
%     baseline(iNeu) = T_iNeu{1, 'Hpre_BitpSec'};
%      
%     T_iNeu = sortrows(T_iNeu, 'dbsFrequency');
%     
%     % determine lowest frequency present within LFS group & Retrieve the 
%     % dataVariable pertaining to lowest LFS trial
%     freqLowest = T_iNeu.dbsFrequency(1);
%     if any(grpLFS == freqLowest)
%         freqs(iNeu,1) = freqLowest;
%         LFS1(iNeu) = T_iNeu{1,'Hdbs1_BitpSec'}; 
%         LFS2(iNeu) = T_iNeu{1,'Hdbs2_BitpSec'};
%         postLFS(iNeu) = T_iNeu{1,'Hpos_BitpSec'};
%         
%     else % if no LFS group frequency trial exists for this neuron
%         freqs(iNeu,1) = NaN;
%         LFS1(iNeu) = NaN;
%         LFS1(iNeu) = NaN;
%         postLFS(iNeu) = NaN;
%         
%     end
%     % determine highest frequency present within HFS group & Retrieve the 
%     % dataVariable pertaining to highest HFS trial
%     freqHighest = T_iNeu.dbsFrequency(end);
%     if any(grpHFS == freqHighest)
%         freqs(iNeu,2) = freqHighest;
%         HFS1(iNeu) = T_iNeu{end,'Hdbs1_BitpSec'};
%         HFS2(iNeu) = T_iNeu{end,'Hdbs2_BitpSec'};
%         postHFS(iNeu) = T_iNeu{end,'Hpos_BitpSec'};
%         
%     else % if no LFS group frequency trial exists for this neuron
%         freqs(iNeu,2) = NaN;
%         HFS1(iNeu) = NaN;     
%         HFS2(iNeu) = NaN;     
%         postHFS(iNeu) = NaN;
%         
%     end
% end
% 
% dataVar = [baseline, LFS1, LFS2, postLFS, HFS1, HFS2, postHFS];
% % Remove any pairs that don't have at least one LFS and one HFS
% remLFS = isnan(LFS1) | isnan(LFS2) | isnan(postLFS);
% remHFS = isnan(HFS1) | isnan(HFS2) | isnan(postHFS);
% freqs(remLFS | remHFS,:) = [];
% dataVar(remLFS | remHFS,:) = [];
% 
% nNeusUpdated = size(dataVar, 1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1657 376];
% ax1 = subplot(1, 3, 1);
% boxplot(dataVar); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'Baseline', 'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(dataVar(:,1:4)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:4;
% ax2.XTickLabel = {'Baseline', 'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot([dataVar(:,1), dataVar(:,5:7)]'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:4;
% ax3.XTickLabel = {'Baseline', 'HFS1', 'HFS2', 'postHFS'};
% 
% ax1.YLim = [0, 500];
% ax2.YLim = [0, 500];
% ax3.YLim = [0, 500];
% 
% 
% % Plot Change in Entropy
% % Add in extra step to get differences in Entropy referenced to Baseline
% diffEntropy = zeros(size(dataVar,1), size(dataVar,2) - 1);
% diffEntropy(:,1) = dataVar(:,2) - dataVar(:,1);
% diffEntropy(:,2) = dataVar(:,3) - dataVar(:,1);
% diffEntropy(:,3) = dataVar(:,4) - dataVar(:,1);
% diffEntropy(:,4) = dataVar(:,5) - dataVar(:,1);
% diffEntropy(:,5) = dataVar(:,6) - dataVar(:,1);
% diffEntropy(:,6) = dataVar(:,7) - dataVar(:,1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1, 3, 1);
% boxplot(diffEntropy); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(diffEntropy(:,1:3)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:3;
% ax2.XTickLabel = {'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot(diffEntropy(:,4:6)'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:3;
% ax3.XTickLabel = {'HFS1', 'HFS2', 'postHFS'};
% 
% ax1.YLim = [-200, 300];
% ax2.YLim = [-200, 300];
% ax3.YLim = [-200, 300];
% 
% 
% % Plot %change in Entropy
% % Add in extra step to get differences in Entropy referenced to Baseline
% pdiffEntropy = zeros(size(dataVar,1), size(dataVar,2) - 1);
% pdiffEntropy(:,1) = 100 * (dataVar(:,2) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,2) = 100 * (dataVar(:,3) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,3) = 100 * (dataVar(:,4) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,4) = 100 * (dataVar(:,5) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,5) = 100 * (dataVar(:,6) - dataVar(:,1)) ./ dataVar(:,1);
% pdiffEntropy(:,6) = 100 * (dataVar(:,7) - dataVar(:,1)) ./ dataVar(:,1);
% 
% 
% 
% f = figure; hold on
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1, 3, 1);
% boxplot(pdiffEntropy); hold on;
% plot([ax1.XLim], [0 0], '--');
% ax1.XTickLabel = {'LFS1', 'LFS2', 'postLFS', 'HFS1', 'HFS2', 'postHFS'};
% title([ppar.subjID, ' ', ppar.trialType, ':  ', num2str(nNeusUpdated), ' neurons']);
% 
% ax2 = subplot(1, 3, 2);
% % plot(dataVar', '-o')
% plot(pdiffEntropy(:,1:3)'); hold on;
% plot([ax2.XLim], [0 0], '--');
% title('LFS')
% ax2.XTick = 1:3;
% ax2.XTickLabel = {'LFS1', 'LFS2', 'postLFS'};
% 
% ax3 = subplot(1, 3, 3);
% % plot(dataVar', '-o')
% plot(pdiffEntropy(:,4:6)'); hold on;
% plot([ax3.XLim], [0 0], '--');
% title('HFS')
% ax3.XTick = 1:3;
% ax3.XTickLabel = {'HFS1', 'HFS2', 'postHFS'};
% 
% ax1.YLim = [-100, 750];
% ax2.YLim = [-100, 750];
% ax3.YLim = [-100, 750];
% 
% 
% 
% 
% %% Display scatter of Depth vs. delta Firing Rate (corrected)
% 
% f = figure;
% f.Position = [1977 595 1405 376];
% ax1 = subplot(1,3,1);
% ax2 = subplot(1,3,2);
% ax3 = subplot(1,3,3);
% 
% % deltaFRdbsCorr1 = Tfinal2{:, 'FRdbsCorr1'} - Tfinal2{:, 'FRpreCorr'};
% % deltaFRdbsCorr2 = Tfinal2{:, 'FRdbsCorr2'} - Tfinal2{:, 'FRpreCorr'};
% % deltaFRposCorr = Tfinal2{:, 'FRposCorr'} - Tfinal2{:, 'FRpreCorr'};
% 
% scatter(ax1, Tfinal2.FRdbs1_idx, Tfinal2.FRpreCorr); grid(ax1, 'on');
% scatter(ax2, Tfinal2.FRdbs2_idx, Tfinal2.FRpreCorr); grid(ax2, 'on');
% scatter(ax3, Tfinal2.FRpos_idx, Tfinal2.FRpreCorr); grid(ax3, 'on');
% 
% ylabel(ax1, 'Pre-DBS firing rate (baseline) spk/sec')
% xlabel(ax2, 'FR change index')
% title(ax1, 'FRdbs1');
% title(ax2, 'FRdbs2');
% title(ax3, 'FRpost');
% 
% % freq sweep Uva
% ax1.XLim = [-1, 1];
% ax2.XLim = [-1, 1];
% ax3.XLim = [-1, 1];
% 
% % % contact sweep Uva
% % ax1.XLim = [-40, 140];
% % ax2.XLim = [-40, 140];
% % ax3.XLim = [-40, 140];
% 
% % all depths
% ax1.YLim = [0, 60];
% ax2.YLim = [0, 60];
% ax3.YLim = [0, 60];
% 
% 
%  
% %% Display scatterHIST of Depth vs. delta Firing Rate (corrected)
% 
% f = figure;
% f.Position = [21 563 1857 420];
% ax1 = subplot(1,3,1);
% ax2 = subplot(1,3,2);
% ax3 = subplot(1,3,3);
% 
% % deltaFRdbsCorr1 = Tfinal2{:, 'FRdbsCorr1'} - Tfinal2{:, 'FRpreCorr'};
% % deltaFRdbsCorr2 = Tfinal2{:, 'FRdbsCorr2'} - Tfinal2{:, 'FRpreCorr'};
% % deltaFRposCorr = Tfinal2{:, 'FRposCorr'} - Tfinal2{:, 'FRpreCorr'};
% 
% scatter(ax1, Tfinal2.FRdbs1_idx, Tfinal2.FRpreCorr); grid(ax1, 'on');
% scatter(ax2, Tfinal2.FRdbs2_idx, Tfinal2.FRpreCorr); grid(ax2, 'on');
% scatter(ax3, Tfinal2.FRpos_idx, Tfinal2.FRpreCorr); grid(ax3, 'on');
% 
% ylabel(ax1, 'Pre-DBS firing rate (baseline) spk/sec')
% xlabel(ax2, 'FR change index')
% title(ax1, 'FRdbs1');
% title(ax2, 'FRdbs2');
% title(ax3, 'FRpost');
% 
% % freq sweep Uva
% ax1.XLim = [-1, 1];
% ax2.XLim = [-1, 1];
% ax3.XLim = [-1, 1];
% 
% % % contact sweep Uva
% % ax1.XLim = [-40, 140];
% % ax2.XLim = [-40, 140];
% % ax3.XLim = [-40, 140];
% 
% % all depths
% ax1.YLim = [0, 60];
% ax2.YLim = [0, 60];
% ax3.YLim = [0, 60];



%% SUB-FUNCTIONS

% function [intervEvents] = getIntervalEvents(evTimes, refT, refInterval)
% % returns the observed firing rate within the desired time interval.
% % "refInterval" indicates the window within which to gather spike times 
% % works with the values in "spkTimes" and "StimTs" as inputs
% 
% rerefTimes = evTimes - refT;
% 
% % get observed spike rate "FRobs" based on count and time duration
% idxInInterv = (rerefTimes >= refInterval(1)) & (rerefTimes < refInterval(2));
% rerefSubselect = rerefTimes(idxInInterv);
% intervEvents = rerefSubselect + refT;
% 
% 
% 
% end
% 
% function [n] = calcFiringRates(spkTimes, StimTs, spkTimeInterval)
% % take the "spkTimes" and "StimTs" data from the nexFile and calculate the
% % observed spike firing rate and other information occurring within the 
% % time interval "spkTimeInterval" (1x2 array, seconds, refenced to DBS 
% % onset time).
% 
% % CONSTANTS
% % assumed to be true for both Uva and Kramer
% periStimBlank = 1 / 1000; %s, blanking time around each stim pulse
% 
% 
% % CODE
% refT = StimTs.DBS(1);
% allStimTs = [StimTs.VirtPre; StimTs.DBS; StimTs.VirtPost];
% 
% % Calculate observed spike rate for interval
% [spkTimesDbs] = getIntervalEvents(spkTimes, refT, spkTimeInterval);
% totIntervTime = spkTimeInterval(2) - spkTimeInterval(1);
% frDbsObs = numel(spkTimesDbs) / totIntervTime;
% 
% 
% % Estimate the "true" firing rate by correcting for DBS artifact blank
% % time, as done in Moran et al 2011
% stmTimesDbs = getIntervalEvents(allStimTs, refT, spkTimeInterval);
% totBlankTime = periStimBlank * numel(stmTimesDbs);
% frDbsCorr = frDbsObs * (totIntervTime / (totIntervTime - totBlankTime));
% 
% 
% n.refT = refT;
% n.spkTimesDbs = spkTimesDbs;
% n.stmTimesDbs = stmTimesDbs;
% n.rateObserved = frDbsObs; % spikes/second
% n.rateCorrected = frDbsCorr; % spikes/second
% n.totIntervTime = totIntervTime;
% n.totBlankTime = totBlankTime;
% n.intervalLimits = spkTimeInterval;
% 
% 
% end


function tabFilt = filter_neuronMinimumTrials(tab, neuronMinimumTrials)
% Remove any neurons from analysis that do not have at least 4 trials
% present in the current selection table

% count the number of times each neuron shows up and store in new table
% called "NeuronCounts"
neurons = tab.Unit_objectID;
[uNeurons, ~, uNeuIdx] = unique(neurons);

nNeurons = length(uNeurons); % number of unique neurons
neuCount = zeros(nNeurons, 1);
for iNeu = 1:nNeurons
    neuCount(iNeu) = sum(uNeuIdx == iNeu);

end

NeuronCounts = [table(uNeurons), table(neuCount)];

% join this table to current selection table
% leftKey = find(strcmp(Nselect.Properties.VariableNames, 'Unit_objectID'));
tab = join(tab, NeuronCounts, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'uNeurons');

% remove those rows with under "n" trials
isPresentforTrials = (tab.neuCount >= neuronMinimumTrials);
tabFilt = tab;
tabFilt(~isPresentforTrials,:) = [];

end

function FRidx = calcFRindex(FRbase, FRtest)
% firing rate index. 0 means no change; 1 means cell went from zero
% baseline to "something"; -1 means cell went from some activity to
% complete inhibition.

FRidx = (FRtest - FRbase) ./ (FRbase + FRtest);

end

function disp_NpointsVsFreq(labels)
% show the numbers of neuron observations for each grouping-label
uniqueLabels = sort(unique(labels));
numLabels = numel(uniqueLabels);

% prefill 

for iLab = 1:numel(uniqueLabels)
    iLabStr = uniqueLabels{iLab};
    num_iLab = sum(strcmp(iLabStr, labels));
    
    disp([iLabStr, ': ', num2str(num_iLab)])
    
    
end


end

function [EmpPval] = getEmpiricalPval(normDistrib, testValue)
% outputs a - or + number that tells how far "value" is from the mean of
% "normDistrib" much like how p-value works. Assumes that the randomly
% distributed data in "normDistrib" is normally distributed, but this is
% not strictly necessary. Interpret with caution. 

isLess = normDistrib < testValue;

Pval = sum(isLess) / numel(normDistrib);

if Pval > 0.5
    EmpPval = 1 - Pval;
    
else
    EmpPval = -Pval;
    
end


end

function [rowCounts] = getCategoryCounts_table(T)

% Gather all counts
TnoCh = T(~T.signifChangeH,:);
numNoCh = height(TnoCh);

Tchange = T(T.signifChangeH,:);

  TdecrH = Tchange(Tchange.wasDecrH,:);
numDecrH = height(TdecrH);

  TincrH = Tchange(Tchange.wasIncrH,:);
numIncrH = height(TincrH);

rowCounts = [numIncrH, numDecrH, numNoCh];

end

function [binCounts] = binSpkTimes(times, binEdges)
% "times" is nx1 vector , and binEdges is 1xm vector. Outputs a 1x(m-1)
% vector for counts of times found within each bin. Assumes that the values
% of "times" and "binEdges" are monotonically increasing. 

nBins = length(binEdges) - 1;

timesInBin = false(length(times), nBins);

for iBin = 1:nBins
    isInBin = (times >= binEdges(iBin)) & (times < binEdges(iBin+1));
    timesInBin(:,iBin) = isInBin;
    
end

binCounts = sum(timesInBin, 1);

end

function [binRates] = getBinRates(spkTimes, refT, refInterval, binWidth)
% Generate bin-counts of spks/sec for this interval
% spkTimes is assumed monotonically increasing, refT is intended to be the
% first DBS pulse time, refInterval is the time over which to look at spike
% times, and binWidth helps build the bin Edges. all units MUST BE in seconds.

% refT = StimTs.DBS(1);
% refInterval = preDbsInterval;
rerefTimes = spkTimes - refT;

% idxInInterv = (rerefTimes >= refInterval(1)) & (rerefTimes < refInterval(2));
% rerefSubselect = rerefTimes(idxInInterv);

% BINWIDTH = 1; % sec
binEdges = refInterval(1):binWidth:refInterval(end);

binCounts = binSpkTimes(rerefTimes, binEdges);

binRates = binCounts / binWidth; 

end

function [ax] = plotDBSlines(ax, xPoint1, xPoint2)
plot(ax, [xPoint1, xPoint1], ax.YLim, 'Color', 'b', 'LineWidth', 2.0);
plot(ax, [xPoint2, xPoint2], ax.YLim, 'Color', 'b', 'LineWidth', 2.0);

end

function [ax] =  plotFRoverTime(dbsCond, Tdisp, freqs)

% dbsCond = 1;
ax = subplot(2,3,dbsCond);
isFr = Tdisp.dbsFrequency == freqs(dbsCond);
T = Tdisp(isFr,:);
nTrials = height(T);
for iTr = 1:nTrials
    plot(T.timeFRall{iTr}, T.binFRallCorr{iTr}); hold on
    
end
% plotDBSlines(ax, 0, 60);
title([num2str(freqs(dbsCond)), 'Hz trials: ', num2str(nTrials)])

end

function [ax] =  plotFRoverTime_norm(dbsCond, Tdisp, freqs)

% dbsCond = 1;
ax = subplot(2,3,dbsCond);
isFr = Tdisp.dbsFrequency == freqs(dbsCond);
T = Tdisp(isFr,:);
nTrials = height(T);
for iTr = 1:nTrials
    plot(T.timeFRall{iTr}, T.binFRallCorr_norm{iTr}); hold on
    
end
% plotDBSlines(ax, 0, 60);
title([num2str(freqs(dbsCond)), 'Hz trials: ', num2str(nTrials)])

end

% function [ax] =  plotAvFRoverTime_norm(dbsCond, Tdisp, freqs)
% 
% % dbsCond = 1;
% ax = subplot(2,3,dbsCond);
% isFr = Tdisp.dbsFrequency == freqs(dbsCond);
% T = Tdisp(isFr,:);
% 
% % get average FR for each bin
% % binAvFRallCorr_norm = 
% 
% nTrials = height(T);
% for iTr = 1:nTrials
%     plot(T.timeFRall{iTr}, T.binFRallCorr_norm{iTr}); hold on
%     
% end
% % plotDBSlines(ax, 0, 60);
% title([num2str(freqs(dbsCond)), 'Hz trials: ', num2str(nTrials)])
% 
% end

function [ax] =  plotFRoverTime_fridx(dbsCond, Tdisp, freqs)

% dbsCond = 1;
ax = subplot(2,3,dbsCond);
isFr = Tdisp.dbsFrequency == freqs(dbsCond);
T = Tdisp(isFr,:);
nTrials = height(T);
for iTr = 1:nTrials
    plot(T.timeFRall{iTr}, T.binFRallCorr_fridx{iTr}); hold on
    
end
% plotDBSlines(ax, 0, 60);
title([num2str(freqs(dbsCond)), 'Hz trials: ', num2str(nTrials)])

end

function [ax] = plotHorizLines(ax, yPoint)
plot(ax, ax.XLim, [yPoint, yPoint], 'Color', 'k', 'LineStyle', '--');

end

function T = tab_gaussFiltDataCell(T, VariableName, GAUSSW)

    nRows = height(T);
    for iRow = 1:nRows
        dataCell = T{iRow, VariableName};
        dataDouble = dataCell{:};
        filtDouble = {filtfilt(GAUSSW, 1, dataDouble)};
        T{iRow, VariableName} = filtDouble;
        
    end
    
end
