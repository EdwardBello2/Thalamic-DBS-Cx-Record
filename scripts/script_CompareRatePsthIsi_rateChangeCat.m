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

ppar.dispRateCateg = 'excite'; % 'excite' | 'inhib' | 'noChange'
ppar.dispPhaseLck = true; % TF



%% LOAD all relevant tables and MERGE them

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



%% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);



%% Calculate Rate changes for selected rows, detect significant rate-changes

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

Tselect = N;



%% Detect those individual trials RECORDS in which neurons phase-locked 
% phase-locked must have a p-value less than ppar.pAlphaPSTH
% detection gets added as an extra column of information about each RECORD

% Add in phase-locked TF value for each row
nRows = size(Tselect, 1);
PhaseLocked = false(nRows, 1);

for iRec = 1:nRows
    if Tselect.pVal_Hpsth(iRec,1) <= ppar.pAlphaPSTH
        PhaseLocked(iRec,1) = true;
        
    end
    
end

Tselect = [Tselect, table(PhaseLocked)];



%% Select just those rows pertaining to ppar.dbsRateChange and ppar.PhaseLocked
% ppar.dispRateCateg = 'excite'; % 'excite' | 'inhib' | 'noChange'
% ppar.dispPhaseLck = true; % TF

% find idx for rate change choice
isRateChange = strcmp(Tselect{:, 'dbsRateChange'}, ppar.dispRateCateg);

% find idx for phaselock choice
isPhsLckChoice = Tselect{:, 'PhaseLocked'} == ppar.dispPhaseLck;

% whittle table
isData = isRateChange & isPhsLckChoice;

T = Tselect(isData,:);



%% Run thru CompareRatePsthIsi display script with table subselection

% CONSTANTS
PREDBS_TIME = ppar.preDbsTime; % seconds
DBS_TIME = ppar.dbsTime; % seconds
BINS_PER_DECADE = ppar.binsPerDecade;
CURR_FUNC = 'script_CompareRatePsthIsi';


% %% LOAD all relevant tables and MERGE them
% 
% % custom script for selecting and merging all tables for this analysis
% T = mergeTables_Master(CURR_FUNC, ppar);
% 
% % re-order table for most significant PSTH p-value to least:
% if isfield(ppar, 'sortRows')
%     T = sortrows(T, ppar.sortRows.byColumn, ppar.sortRows.sortType);
%     
% end



% MAIN code

% FOR each record, plot the rate over time, the PSTH change, and the
% log-ISI change:


binWidth = 1; %seconds
nRecs = size(T, 1);
for iRec = 1:nRecs
    %% PLOT Instantaneous Rate
    
    % load the variable "bins"
    load([ppar.projRootPath, '\', T.matFileFolder_RateBins{iRec}, '\', ...
        T.matFile_RateBins{iRec}]);
    
    % prep data for inst. rate display
    rates = [bins.pre; bins.dbs; bins.pos];
    t = (0:(numel(rates) - 1)) * binWidth;
    
    
    % prep data for rate-change test
    ratesPRE = bins.pre / binWidth;
    ratesDBS = bins.dbs / binWidth;
    
    [pValRate, h, stats] = ranksum(ratesPRE, ratesDBS);
%     [pVal, h, stats] = ranksum(ratesPRElog, ratesDBSlog);
    
    
    % plot results
    f = figure; 
    f.Position = [1992 305 945 597];

    ax1 = subplot(2, 2, 1:2);
    % ax = axes;
    plot(t, rates, 'k');
    xlabel('time in dbs-trial (seconds)');
    ylabel('inst. spikes/sec');
    title([T.objectID{iRec}, ', ', ...
        num2str(T.dbsFrequency(iRec,1)), ' Hz DBS    ', ...
        'Spks/sec   PreVsDbs: ', ...
        num2str(mean(ratesPRE)), ' vs ', num2str(mean(ratesDBS)), ...
        ', p = ', num2str(pValRate)], 'Interpreter', 'none')

    % also plot a little horizontal line above showing when DBS was on
    nPre = numel(bins.pre);
    nDbs = numel(bins.dbs);
    nPos = numel(bins.pos);
    dbsBeg = nPre + 1;
    dbsEnd = nPre + nDbs + 1;

    hold on;
    ax1.XMinorGrid = 'on';
    ax1.XMinorTick = 'on';

    line([dbsBeg, dbsEnd], [0.9, 0.9] * ax1.YLim(2), 'Color', 'r', 'LineWidth', 2);



    %% PLOT PSTH

    nexFile = readNexFile([ppar.projRootPath, '\', T.nexFileFolder{iRec}, ...
        '\', T.nexFile{iRec}]);
    [spkTimes, stims] = parseNexFile(nexFile);

    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = stims.DBS;
    stimPeriod = median(diff(dbsTimes));


    % get pre-DBS spikes
    isPreDBS = (spkTimes < dbsTimes(1)) & ...
               (spkTimes >= (dbsTimes(1) - PREDBS_TIME));
    spksPRE = spkTimes(isPreDBS);


    % get DBS-on spikes
    isDBSon = (spkTimes >= dbsTimes(1)) & ...
              (spkTimes < (dbsTimes(end) + stimPeriod));
    spksDBS = spkTimes(isDBSon);



    % CALCULATE PSTH of both PRE and DBS conditions
    tBeg = ppar.psthTimeBeg; % seconds
    bw = ppar.psthBinWidth; % seconds
    nBins = ppar.psthNumBins;
    binEdges = tBeg:bw:(bw * nBins);

    psthPRE = psth(spksPRE, stims.VirtPre, binEdges);
    psthDBS = psth(spksDBS, stims.DBS, binEdges);

    % If user specified to remove first and last bins, perform that now:
    if ppar.trimPSTH
        psthDBS = psthDBS(2:(end-1));
        psthPRE = psthPRE(2:(end-1));
        binEdges = binEdges(2:(end-1));
        
    end

    [xSpre, ySpre] = histcounts2stairplot(psthPRE, binEdges);
    [xSdbs, ySdbs] = histcounts2stairplot(psthDBS, binEdges);


    HprePsth = T.H_PREemp_Hpsth(iRec);
    HdbsPsth = T.H_DBSemp_Hpsth(iRec);

    ax2 = subplot(2, 2, 3);
    plot(xSpre, ySpre); hold on;
    plot(xSdbs, ySdbs);
    xlabel('time after DBS pulse (seconds)');
    ylabel('Probability of spike occurrence');
    title(['PSTH, ', ...
    'pVal = ', num2str(T.pVal_Hpsth(iRec)), ', ', ...
    'Hpre = ', num2str(HprePsth), ', ', ...
    'Hdbs = ', num2str(HdbsPsth), ', ', ...
    '%diff: ', num2str(100 * (HdbsPsth - HprePsth) / HprePsth)], ...
    'Interpreter', 'none')
    legend('PRE-dbs', 'DBS-on')



    %% PLOT ISI

    % Load each nexfile one at a time
    nexFile = readNexFile([ppar.projRootPath, '\', T.nexFileFolder{iRec}, ...
        '\', T.nexFile{iRec}]);
    [spkTimes, StimTs] = parseNexFile(nexFile);

    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = StimTs.DBS;
    stimPeriod = median(diff(dbsTimes));

    % get PRE-dbs spikes
    isPreDBS = (spkTimes >= (dbsTimes(1) - PREDBS_TIME)) & ...
              (spkTimes < dbsTimes(1));
    spksPRE = spkTimes(isPreDBS);

    % Define ISIs, and remove any ISIs == 0
    isiPRE = diff(spksPRE);
    isiPRE(isiPRE == 0) = [];

    % log-bin the isi's, function from "Entropy" toolbox
    [isiLogHist, binEdges, ~] = isiLogBinned(isiPRE, BINS_PER_DECADE);
    % convert from histogram counts to histogram probability:
    isiLogProb = isiLogHist ./ sum(isiLogHist); 

    % convert its data to step-Histogram-data
    [xSpre, ySpre] = histcounts2stairplot(isiLogProb, binEdges);




    % get DBS-on spikes
    isDBSon = (spkTimes >= dbsTimes(1)) & ...
              (spkTimes < (dbsTimes(end) + stimPeriod));
    spksDBS = spkTimes(isDBSon);

    % Define ISIs, and remove any ISIs == 0
    isiDBS = diff(spksDBS);
    isiDBS(isiDBS == 0) = [];

    % log-bin the isi's, function from "Entropy" toolbox
    [isiLogHist, binEdges, ~] = isiLogBinned(isiDBS, BINS_PER_DECADE);
    % convert from histogram counts to histogram probability:
    % isiLogProb = 100 * isiLogHist ./ sum(isiLogHist);
    isiLogProb = isiLogHist ./ sum(isiLogHist);


    % convert its data to step-Histogram-data
    [xSdbs, ySdbs] = histcounts2stairplot(isiLogProb, binEdges);


    Hpre = T.H_PREemp_Hisi(iRec);
    Hdbs = T.H_DBSemp_Hisi(iRec);


    % plot both in one figure
    ax3 = subplot(2, 2, 4);
    % f.Position = [680 535 788 443];
    plot(xSpre, ySpre); hold on
    plot(xSdbs, ySdbs);

    ax3.XScale = 'log';
    ax3.XLim = [0.001, 1]; % one ms to one second
    xlabel('log-ISI bins (seconds)');
    ylabel('isi Probability')

    title(['ISI, ', ...
        'Hpre = ', num2str(Hpre), ', ',...
        'Hdbs = ', num2str(Hdbs), ', ',...
        '%diff: ', num2str(100 * (Hdbs - Hpre) / Hpre)], ...
        'Interpreter', 'none')
    legend('PRE-dbs', 'DBS-on')

    
    
    disp([T.objectID{iRec}, ...
          ', file ', num2str(iRec), ' of ', num2str(nRecs)])

    pause()
    close(f)


end




