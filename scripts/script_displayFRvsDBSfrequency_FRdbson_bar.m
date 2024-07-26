% script for looking at progression of DBS-on neuron Firing Rates as 
% Z-scores relative to baseline distribution of firing rateslooking at 
% creating bootstrapped entropy comparisons for 30-sec stretches of DBS 
% trial
% NOTE: still need to put in method for Contact-sweep too

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
ppar.subjID = 'Uva'; % 'Kramer' | 'Uva' | 'bothSubj'
ppar.trialType = 'FrequencySweep'; % 'FrequencySweep' | 'ContactSweep'
% ppar.groups.dbsFrequency.groupFreqs = {10, 20, 30, 50, 100, 130};
% ppar.groups.dbsFrequency.groupLabels = {'hz10', 'hz20', 'hz30', 'hz50', 'hz100', 'hz130'};



%%  CONSTANTS

% Perform analysis on either first of last 30 seconds of the 60-seconds of
% DBS
FRvar = 'FRdbsCorr2'; % str, 'FRdbsCorr1' | 'FRdbsCorr2'


% Spike-sort quality:
REFRAC = 1 / 1000; % seconds, neuron refractory period in which spkes shouldn't occur
PERC_REFRACTHRESH = 1.0; % threshold for rejecting a unit with x% spikes within the refractory period


% Filter the data:
NEURON_MIN_TRIALS = 3;
SPARSE_HZ = 2;


% Define time ranges of various peri-DBS events: (in seconds, make sure 2nd
% element > 1st element
  SPKINTERV_PRE = [-30, 0]; % [-30, 0] | [-Inf, 0]
SPKINTERV_DBS_1 = [0, 30];   % [0, 30]
SPKINTERV_DBS_2 = [30, 60];  % [30, 60]
  SPKINTERV_POS = [60, 90]; % [60, 90] | [60, Inf]
  % NOTE: the -31 and 91 is to make it possible to collect bins for -30 and
  % 90 seconds...

BINWIDTH = 1; % seconds

% % Entropy calculation:
% BINS_PER_DECADE = 15;
% ORD_H = 1;
% NBOOTS = 10000;

% change significance for categorizing changes in decrH / incrH labels
SIG_PVAL = 0.001 / 2; % to account for two-tailed distributions


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
gaussFiltFR = false;



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
    
% Append data of interest to Tselect
Tfinal = [Tselect, ...
          table(FRpreCorr), table(FRdbsCorr1), table(FRdbsCorr2), table(FRposCorr), ...
          table(timeFRall), table(binFRallCorr)];
      


%% Prepare data for final visualization and analysis

numel(unique(Tfinal.Unit_objectID)) % display current number units       
Tfinal2 = Tfinal;

% Tanalyze = Tfinal;

% Remove any rows that have a firing rate of less than min firing rate:
% filter out any rows that have too-sparse cells, based on 30-second
% portions
isPreRateBelow = Tfinal2{:, 'FRpreCorr'} < SPARSE_HZ;
% isDbs1RateBelow = Tfinal2{:, 'FRdbsCorr1'} < SPARSE_HZ;
isDbs2RateBelow = Tfinal2{:, 'FRdbsCorr2'} < SPARSE_HZ;
% isRemove = isPreRateBelow & isDbs1RateBelow & isDbs2RateBelow;
isRemove = isPreRateBelow & isDbs2RateBelow;

%  isRemove = isPreRateBelow;

Tfinal2(isRemove, :) = [];
numel(unique(Tfinal2.Unit_objectID))


% % Remove any rows that have neurons that were not recorded with enough time
% % as specified by tStard and tEnd (where t = 0 is onset of DBS). Note this
% % is not the same as excluding neurons that were recorded to have a FR of
% % 0.
% tStart = SPKINTERV_PRE(1) + 1;
% tEnd = SPKINTERV_POS(2) - 1;
% timeFRstandard = tStart:BINWIDTH:tEnd;
% nTimeBins = numel(timeFRstandard);
% nRows = height(Tfinal2);
% hasEnoughTimeBins = false(nRows, 1);
% T = Tfinal2;
% for iRow = 1:nRows
%     % check if time bins for this row contains the requisite start and end
%     % times
%     timeBins = T.timeFRall{iRow,1};
%     if (timeBins(1) <= tStart) && (timeBins(end) >= tEnd)
%         hasEnoughTimeBins(iRow) = true;
%         
%     end
%         
% end
% T(~hasEnoughTimeBins,:) = [];
% disp(['Reomved ', num2str(sum(~hasEnoughTimeBins)), ' of ', ...
%        num2str(nRows),' table rows due to not enough time recorded']);
% Tfinal2 = T;


% Remove any neurons that have less than minimum trials 
Tfinal3 = filter_neuronMinimumTrials(Tfinal2, NEURON_MIN_TRIALS);
numel(unique(Tfinal3.Unit_objectID))

Tanalyze = Tfinal3;


% Add in average FR line from -30 to 60 seconds, for each dbsFrequency


% Reassign tStart to the shortest pre-DBS start time in the data, so that
% all data will have equal length




%% Display FR averages, group axes by DBS frequency

Tdisp = Tanalyze;
freqs = [10, 20, 30, 50, 100, 130];
% f1 = figure;
% 
% f1.Position = [1939 80 1833 903];

% if option is checked, gauss-filter all time-series data
if gaussFiltFR % smooth all the data itself
    Tdisp = tab_gaussFiltDataCell(Tdisp, 'binFRallCorr', GAUSSW);

end

% Calculate DBS-on FR as a Z-score relative to distribution of baseline FR
% bins
T = Tdisp;
nRows = height(T);
dbsFRz = zeros(nRows, 1);
for iRow = 1:nRows
    % for row 1:
    % Grab all FR bins pertaining to pre-DBS period
    timeFRall = T.timeFRall{iRow};
    idxPre = (timeFRall >= SPKINTERV_PRE(1)) & (timeFRall <= SPKINTERV_PRE(2));
    binFRallCorr = T.binFRallCorr{iRow};
    binFRpre = binFRallCorr(idxPre);

    % Get average DBS-on FR as a Z-score 
%     dbsFR = T.FRdbsCorr1(iRow);
    dbsFR = T.FRdbsCorr2(iRow);
    
    distrib = [dbsFR, binFRpre];
    Z = zscore(distrib);
    dbsFRz(iRow) = Z(1);
    % figure; histogram(binFRpre); title([num2str(dbsFR)])

end
T = [T, table(dbsFRz)];
Tdisp = T;

% Gather baseline pre-DBS FR values for each single-unit
T = Tdisp;
units = unique(T.Unit_objectID);
nUnits = numel(units);
baseFR = zeros(nUnits, 1);
for iU = 1:nUnits
    % for a given unit, arrange all recordings from earliest to latest in
    % the session
    idxiU = strcmp(units(iU), T.Unit_objectID);
    T_iU = T(idxiU,:);
    sortrows(T_iU, 'blockNum');

    % from the earliers block, grab the average pre-DBS FR as baseline for analysis
    baseFR(iU) = T_iU.FRpreCorr(1);

end
Tdisp = T;


% Group FR values by DBS frequency
% FRvar = 'FRdbsCorr1';
T = Tdisp;
avDbsFRz = zeros(1, 6);
semDbsFRz = zeros(1, 6);
for ifreq = 1:6
    % get sub-selection of table pertaining to current frequency
    idxFreq = T.dbsFrequency == freqs(ifreq);
    T_ifreq = T(idxFreq,:);
    
    avDbsFRz(ifreq) = mean(T_ifreq{:, FRvar});
    semDbsFRz(ifreq) = std(T_ifreq{:, FRvar}, [], 1) ./ sqrt(length(T_ifreq{:, FRvar}));
    
end

% Display FR results as barplots
avBaseFR = mean(baseFR);
semBaseFR = std(baseFR, [], 1) ./ sqrt(length(baseFR));
avFR = [avBaseFR, avDbsFRz];
semFR = [semBaseFR, semDbsFRz];
figure; ax1 = axes;
barwitherr(semFR, avFR)
ax1.XTickLabel = {'base', '10', '20', '30', '50', '100', '130'};
% ax1.YLim = [-1, 1];
tit = title([ppar.subjID, ': FR ']);
xlabel('DBS frequency (Hz)');
ylabel('DBS-on firing rate (spikes/second)');

% ANOVA between group nominal variables

% concatenate all FR values grouped by DBS frequency with baseline too
FRdbs = T{:, FRvar};
grpDouble = T{:, 'dbsFrequency'};
grpdbsStr = string(grpDouble);

FR = [FRdbs; baseFR];

baseStr = strings(length(baseFR), 1);
baseStr(:) = 'base';
grpStr = [grpdbsStr; baseStr];

[p, tbl, stats] = anova1(FR, grpStr);
figure;
[c, m, h, gnames] = multcompare(stats)

% two-sample t-tests comparing each DBS-on FR distrib with baseline distrib
% check for allHz case;
freqs = [10, 20, 30, 50, 100, 130];
for ifreq = 1:6   
    % get sub-selection of table pertaining to current frequency
    idxFreq = T.dbsFrequency == freqs(ifreq);
    T_ifreq = T(idxFreq,:);
    
    avDbsFR = T_ifreq{:, FRvar};
    [h, p, ci, stats] = ttest2(avDbsFR, baseFR);
    
    
    statTable{ifreq,1} = p;
    statTable{ifreq,2} = ci(1);
    statTable{ifreq,3} = ci(2);
    statTable{ifreq,4} = stats.tstat;
    statTable{ifreq,5} = stats.df;
    statTable{ifreq,6} = stats.sd;
    
end



%% Additional rate-increase/decrease classification step for each DBS freq
% does DBS cause a rate change relative to the immediately preceding preDBS
% period? 

T = Tdisp;
nRows = height(T);
pValranksum = zeros(nRows, 1);
  FRsigDiff = false(nRows, 1);
    wasIncrFR = false(nRows, 1);
    wasDecrFR = false(nRows, 1);
   noChangeFR = false(nRows, 1);
   FRchangeLabel = cell(nRows, 1);
for iRow = 1:nRows
    timeFRall = T.timeFRall{iRow};

    % Grab all FR bins pertaining to pre-DBS period
    idxPre = (timeFRall >= SPKINTERV_PRE(1)) & (timeFRall <= SPKINTERV_PRE(2));
    binFRallCorr = T.binFRallCorr{iRow};
    binFRpre = binFRallCorr(idxPre);

    % Grab all FR bins pertaining to DBS-on period
    idxDbs = (timeFRall >= SPKINTERV_DBS_2(1)) & (timeFRall <= SPKINTERV_DBS_2(2));
    binFRallCorr = T.binFRallCorr{iRow};
    binFRdbs = binFRallCorr(idxDbs);

    % use Wilcoxon ranksum test to see if sig diff from predbs to dbson FR for
    % each row
    pValranksum(iRow) = ranksum(binFRpre, binFRdbs);
    if pValranksum(iRow) < SIG_PVAL
        FRsigDiff(iRow) = true;
        
    end

    % also update counts for each response type
    if FRsigDiff(iRow)
        if mean(binFRdbs) - mean(binFRpre) > 0 % if DBS cause increase FR
            wasIncrFR(iRow) = true;
            FRchangeLabel{iRow} = 'wasIncrFR';

        else
            wasDecrFR(iRow) = true;
            FRchangeLabel{iRow} = 'wasDecrFR';
            
        end
             
    else
        noChangeFR(iRow) = true;
        FRchangeLabel{iRow} = 'noChange';
    end
       
end

% add results to table
T = [T, table(pValranksum), table(FRsigDiff), table(wasIncrFR), ...
        table(wasDecrFR), table(noChangeFR), table(FRchangeLabel)];
    
Tdisp = T;

% create table columns of TF values for each type of response

% display relative % proportions of each type of response vs. dbsFrequency
% pre-allocate matrix for counts with each frequency being a row; columns:
% wasIncrFR, noChange, wasDecrFR 
freqs = [10, 20, 30, 50, 100, 130];
nFreqs = numel(freqs);
catCounts = zeros(nFreqs, 3);

% Gather all counts
for iFr = 1:nFreqs
    % get subselection of table rows pertaining to frequency "iFr"
    T_byFreq = Tdisp(Tdisp.dbsFrequency == freqs(iFr),:);
    
    % categorize each trial into "incrH", "decrH", "noCh" (in that order), 
    % keep count
    
    catCounts(iFr,:) = tab_getVariableCounts(T_byFreq, {'wasIncrFR', 'noChangeFR', 'wasDecrFR'});
  
end
 
% transfer to percentages
tot = sum(catCounts, 2);
catPercs = 100 * catCounts ./ tot;
% catPercs = catPercs';

% display results as barplots
figure; ax = axes;
bar(catPercs, 'stacked')
ax.XTickLabel = freqs;
grid on
xlabel('DBS frequency (Hz)')
ylabel('% recorded cells')
lgd = legend('FR incr', 'noChange', 'FR decr', ...
      'Location', 'northeastoutside') ;
tit = title([ppar.subjID, ':  population percentage by 30-60s DBS Entropy change, p < ', num2str(SIG_PVAL*2)]);
ax.YLim = [0, 100];

% % alternate display as lines, noChange category excluded
% figure; ax = axes;
% plot(100 * catPercs(:,1), 'o--'); hold on;
% plot(100 * catPercs(:,2), 'o--');
% ax.YLim = [0, 100];
% ax.XLim = [0.5, 6.5];
% ax.XTickLabel = freqs;
% grid on
% xlabel('DBS frequency (Hz)')
% ylabel('% recorded cells')
% lgd = legend('FR decr', 'FR incr', ...
%       'Location', 'northeastoutside') ;
% tit = title([ppar.subjID, ':  population percentage by 30-60s DBS Entropy change, p < ', num2str(SIG_PVAL*2)]);
% 


% Chi-2 test for independence; add in one more column with group label
Ttest = Tdisp;
[xtab, chi2, p, labels] = crosstab(Ttest.dbsFrequency, Ttest.FRchangeLabel)


% Chi-2 test for independence post-hoc test: FRincr vs. notFRincr
disp('Chi-2 test for FRincr vs. notFRdecr');
Ttest = Tdisp;
FRsigDecr = cell(nRows, 1);
for iRow = 1:nRows
    if Ttest.FRsigDiff(iRow)
        if Ttest.wasDecrFR(iRow)
            FRsigDecr{iRow} = 'FRdecr';
            
        else
            FRsigDecr{iRow} = 'notFRdecr';
            
        end
        
    else
        FRsigDecr{iRow} = 'notFRdecr';
        
    end
       
end
Ttest = [Ttest, table(FRsigDecr)];
[xtab, chi2, p, labels] = crosstab(Ttest.dbsFrequency, Ttest.FRsigDecr)



% Chi-2 test for independence post-hoc test: FRdecr vs. notFRdecr
disp('Chi-2 test for FRincr vs. notFRincr');
Ttest = Tdisp;
FRsigIncr = cell(nRows, 1);
for iRow = 1:nRows
    if Ttest.FRsigDiff(iRow)
        if Ttest.wasIncrFR(iRow)
            FRsigIncr{iRow} = 'FRincr';
            
        else
            FRsigIncr{iRow} = 'notFRincr';
            
        end
        
    else
        FRsigIncr{iRow} = 'notFRincr';
        
    end
       
end
Ttest = [Ttest, table(FRsigIncr)];
[xtab, chi2, p, labels] = crosstab(Ttest.dbsFrequency, Ttest.FRsigIncr)



%% % Gather baseline pre-DBS FR values for each single-unit
% see if pre-DBS FR has a consistent decrease over time trend in most/all
% cases
% update: It doesn't. preDBS rates are all over the frikkin place over time
T = Tdisp;

units = unique(T.Unit_objectID);
nUnits = numel(units);
baseFR = zeros(nUnits, 1);
figure; ax = axes;
for iU = 1:nUnits
    % for a given unit, arrange all recordings from earliest to latest in
    % the session
    idxiU = strcmp(units(iU), T.Unit_objectID);
    T_iU = T(idxiU,:);
    sortrows(T_iU, 'blockNum');

    % from the earliers block, grab the average pre-DBS FR as baseline for analysis
    baseFR(iU) = T_iU.FRpreCorr(1);
    plot(T_iU.FRpreCorr)
    hold on

end
ax.XLim = [0, 7];
ax.XTick = [1,2,3,4,5,6];


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
% 
% function [rowCounts] = tab_getVariableCounts(T, varStrCell)
% % varStrCell must be a cell vector of chars
% 
% nVars = length(varStrCell);
% rowCounts = zeros(1, nVars);
% 
% for iVar = 1:nVars
%     varStr = varStrCell{iVar};
%     iVarTF = T{:, varStr};
%     rowCounts(iVar) = sum(iVarTF);
%     
% end
% 
% totCount = sum(rowCounts);
% 
% if ~(totCount == height(T))
%     error('make sure to include only TF values groups that together sum up to total num of table rows');
% 
% end
% 
% end
% 

