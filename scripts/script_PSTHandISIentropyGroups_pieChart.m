% script for looking at PSTH-based entropy vs DBS frequency
% NOTE: still need to put in method for Contact-sweep too

% Cited:
% Moran, A., Stein, E., Tischler, H., Belelovsky, K. & Bar-Gad, 
% I. Dynamic Stereotypic Responses of Basal Ganglia Neurons to Subthalamic
% Nucleus High-Frequency Stimulation in the Parkinsonian Primate. 
% Frontiers in Systems Neuroscience 5, 21–21 (2011).

% To-Do:
% - currently the log(entropy) calculations are messed up due to some
% bootstrapped entropuies being "0", giving -inf and NaN problems down the
% road. Must put in step to remove such cases from the distribution...



%% build pipeline-parameter struct "ppar"

% clear; 

% PIPELINE PARAMETERS
% full path on local PC where project folder is (don't include subfolders here,
% that's tracked within the appropriate tables)
ppar.projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record'; % string

% full path on local PC where tables are to be loaded from (or saved to)
% ppar.tablePath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string
ppar.tablePath = 'C:\Users\bello043\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string

ppar.intDataRootFolder = 'DataProcessing\intermediateData';

% controls whether or not intermediate data for this run is simply done
% over again, even if the data already exists. 
ppar.overwriteIntData = false ; % TF: true | false

% For getting sub-selection of data from tables
ppar.neuronType = 'SU';
ppar.subjID = 'bothSubj'; % 'Kramer' | 'Uva' | 'bothSubj'
ppar.trialType = 'FrequencySweep'; % 'FrequencySweep' | 'ContactSweep'
% ppar.groups.dbsFrequency.groupFreqs = {10, 20, 30, 50, 100, 130};
% ppar.groups.dbsFrequency.groupLabels = {'hz10', 'hz20', 'hz30', 'hz50', 'hz100', 'hz130'};

% PSTH parameters
ppar.psthTimeBeg = 0; % seconds
ppar.psthBinWidth = 0.5 / 1000; %seconds
ppar.psthNumBins = 15;
ppar.pAlphaPSTH = 0.05;
% TF indicating whether to remove the first and last bins
ppar.trimPSTH = [1,2]; % [numerical indices] | [] ; note: ONLY trim bins at the edges, NOT any in the middle (will mess up code)


bBeg = ppar.psthTimeBeg;
bw = ppar.psthBinWidth;
nBins = ppar.psthNumBins;
binEdgesPsth = bBeg:bw:(bw * nBins);


% CONSTANTS

% Spike-sort quality:
REFRAC = 1 / 1000; % seconds, neuron refractory period in which spkes shouldn't occur
PERC_REFRACTHRESH = 1.0; % threshold for rejecting a unit with x% spikes within the refractory period


% Filter the data:
NEURON_MIN_TRIALS = 3;
SPARSE_HZ = 2;



% Define time ranges of various peri-DBS events:
SPKINTERV_PRE = [-30, 0]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS_1 = [0, 30]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS_2 = [30, 60]; %seconds, make sure 2nd element > 1st
SPKINTERV_POS = [60, 90]; %seconds, make sure 2nd element > 1st

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

% Initialize columns to add to Tselect for final analysis in this script:
 FRpreCorr = zeros(nNex, 1);
      Hpre_bitspk = zeros(nNex, 1);
FRdbsCorr1 = zeros(nNex, 1);
     Hdbs1 = zeros(nNex, 1);
FRdbsCorr2 = zeros(nNex, 1);
     Hdbs_bitspk = zeros(nNex, 1);
   Hdbs_bitspkZ = zeros(nNex, 1);
  H_DBSemp = zeros(nNex, 1);

    FRposCorr = zeros(nNex, 1);
         Hpos = zeros(nNex, 1);
      Hboot_bitspkAv = zeros(nNex, 1);
    Hboot_bitspkStdv = zeros(nNex, 1);
Hdbs_bitspk_EmpPval = zeros(nNex, 1);

trimIdx = false(length(binEdgesPsth), 1);
trimIdx(ppar.trimPSTH) = true;
HbootDistrib_bitspk = cell(nNex, 1);
psth_Ct = cell(nNex, 1);
psth_nSpk = zeros(nNex, 1);
psthFR = zeros(nNex, 1);

Hdbs_bitsec = zeros(nNex, 1);
Hboot_bitsecAv = zeros(nNex, 1);
Hboot_bitsecStdv = zeros(nNex, 1);
HbootDistrib_bitsec = cell(nNex, 1);
    
Hdbs_bitpulse = zeros(nNex, 1);
HbootDistrib_bitpulse = cell(nNex, 1);
Hboot_bitpulseAv = zeros(nNex, 1);
Hboot_bitpulseStdv = zeros(nNex, 1);
Hboot_bitpulseLogAv = zeros(nNex, 1);
Hboot_bitpulseLogStdv = zeros(nNex, 1);

      
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
    fs = nexFile.freq; % samples/second

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
    formatSpecFn = 'FRpre_%s_%ss-%ss_fromDbsOnset';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_PRE(1)), ...
                              num2str(SPKINTERV_PRE(2)));
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
        save(fullPathFn, 'FRpre');
    
    end
    
    FRpreCorr(iNex,1) = FRpre.rateCorrected;

    

    %% Get DBS-ON period spike rate, both observed and corrected for artifact
    % blanking for first 30 sec

    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'FRdbs';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'FRdbs_%s_%ss-%ss_fromDbsOnset';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_DBS_1(1)), ...
                              num2str(SPKINTERV_DBS_1(2)));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
        
        FRdbs = calcFiringRates(spkTimes, StimTs, SPKINTERV_DBS_1);
        FRdbs.intervalLimits = SPKINTERV_DBS_1;
        save(fullPathFn, 'FRdbs');
    
    end
    
    FRdbsCorr1(iNex,1) = FRdbs.rateCorrected;

    
    
    %% Get DBS-ON period spike rate, both observed and corrected for artifact
    % blanking for SECOND 30 sec

    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'FRdbs';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'FRdbs_%s_%ss-%ss_fromDbsOnset';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_DBS_2(1)), ...
                              num2str(SPKINTERV_DBS_2(2)));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
        
        FRdbs = calcFiringRates(spkTimes, StimTs, SPKINTERV_DBS_2);
        FRdbs.intervalLimits = SPKINTERV_DBS_2;
        save(fullPathFn, 'FRdbs');
    
    end
    
    FRdbsCorr2(iNex,1) = FRdbs.rateCorrected;
    
    
    
    %% Get POST-DBS period spike rate, both observed and corrected for artifact
    % blanking

    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'FRpos';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'FRpos_%s_%ss-%ss_fromDbsOnset';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_POS(1)), ...
                              num2str(SPKINTERV_POS(2)));
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
        save(fullPathFn, 'FRpos');
    
    end
    
    FRposCorr(iNex,1) = FRpos.rateCorrected;


     %% Calculate entropy for onDBS second 30 sec, by bit/spike
   
    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'Hdbs';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'Hdbs_%s_%smsBW_%sBins_%sboots_%s';
    if ppar.trimPSTH
        trimStr = 'Trim';
        
    else
        trimStr = 'noTrim';
        
    end
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(ppar.psthBinWidth*1000), ...
                              num2str(ppar.psthNumBins), ...
                              num2str(NBOOTS), ...
                              trimStr);
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate entropies for this nexfile
    if (0 < exist(fullPathFn, 'file')) && (~ppar.overwriteIntData)
        load(fullPathFn)
        
    else
         % Get spike times occurring within DBS interval 30-60
        refT = StimTs.DBS(1);
        
        % get spike times in relevent interval
        [spkTimesDbs] = getIntervalEvents(spkTimes, refT, SPKINTERV_DBS_2);
        spkTimesDbs = sort(spkTimesDbs); % make sure the seconds are monotonically increasing...
        
        % get DBS times in relevent interval
        dbsTimes = getIntervalEvents(StimTs.DBS, refT, SPKINTERV_DBS_2);
        
        % Calculate empirical psth-Entropy for DBS
        if isempty(spkTimesDbs)
            psth_DBSct = zeros(1, length(binEdgesPsth)-1);
            
        else
            psth_DBSct = psth(spkTimesDbs, dbsTimes, binEdgesPsth, 'count');
            
        end
        % If user specified to remove first bin, perform that now:
        if ppar.trimPSTH % if trimPSTH is '[]', this step is skipped
            psth_DBSct(ppar.trimPSTH) = 0;

        end
        nSpksDBS = sum(psth_DBSct);
        psth_DBSprob = psth_DBSct / nSpksDBS; % norm by bin probability
        psth_DBSspksec = psth_DBSct / (numel(dbsTimes) * bw); % norm by spikes per second
        H_DBSemp = entropyLetter_bitpSpike(psth_DBSprob);
        
        
        % Calculate bootstrapped distribution of psth-Entropies for DBS
        H_DBSboots = zeros(NBOOTS, 1);
        unifRange = binEdgesPsth;
        unifRange(trimIdx) = [];
        A = unifRange(1);
        B = unifRange(end);
        
        if (nSpksDBS < 1) || isempty(nSpksDBS)
            H_DBSboots(:) = NaN;
            
        else
            parfor iBoot = 1:NBOOTS
                % for each bootstrap PSTH: based on current "allowed" bins (trimming?), 
                % random-uniformly sample spike counts for the bins
%                 psth_iBOOTct = zeros(length(binEdgesPsth)-1, 1)
                R = unifrnd(A, B, nSpksDBS, 1);
                [psth_iBOOTct,~] = histcounts(R, binEdgesPsth);
%                 psth_iBOOTct = h.Values;
       
%                 timeRange = dbsTimes(1):(1/fs):(dbsTimes(end));
%                 [dbsTimes_iBoot] = datasample(timeRange, numel(dbsTimes));
% 
%                 % Calculate empirical psth-Entropy for DBS
%                 psth_DBSct = psth(spkTimesDbs, dbsTimes_iBoot, binEdgesPsth, 'count');
%                 % If user specified to remove first bin, perform that now:
%                 if ppar.trimPSTH % if trimPSTH is '[]', this step is skipped
%                     psth_DBSct(ppar.trimPSTH) = 0;
% 
%                 end
                psth_iBOOTprob = psth_iBOOTct / nSpksDBS;

                % gather bootstrapped entropy distribution
                H_DBSboots(iBoot) = entropyLetter_bitpSpike(psth_iBOOTprob);

            
            end

        end
        
            
        
        
        
        
        
        
        
        % Calculate the PSTH of pre-DBS period as well
        % get spike times in relevent interval
        [spkTimesPre] = getIntervalEvents(spkTimes, refT, SPKINTERV_PRE);
        spkTimesPre = sort(spkTimesPre); % make sure the seconds are monotonically increasing...
        
        % get DBS times in relevent interval
        dbsTimesVirt = getIntervalEvents(StimTs.VirtPre, refT, SPKINTERV_PRE);
        
        % Calculate empirical psth-Entropy for DBS
        if isempty(spkTimesPre)
            psthPREct = zeros(1, length(binEdgesPsth)-1);
            
        else
            psthPREct = psth(spkTimesPre, dbsTimesVirt, binEdgesPsth, 'count');
            
        end
        
        % If user specified to remove first bin, perform that now:
        if ppar.trimPSTH % if trimPSTH is '[]', this step is skipped
            psthPREct(ppar.trimPSTH) = 0;

        end
        nSpksPRE = sum(psthPREct);
        psth_PREprob = psthPREct / nSpksPRE;
        psth_PREspksec = psthPREct / (numel(dbsTimesVirt) * bw); % norm by spikes per second
        H_PREemp = entropyLetter_bitpSpike(psth_PREprob);

        
        save(fullPathFn, 'H_DBSemp', 'H_DBSboots', 'H_PREemp', ...
             'psth_DBSprob', 'psth_DBSspksec', 'psth_DBSct', 'psth_PREprob', 'psth_PREspksec');
    
    end

    Hdbs_bitspk(iNex,1) = H_DBSemp;
    Hpre_bitspk(iNex,1) = H_PREemp;
    H_DBSboots(isnan(H_DBSboots)) = []; % remove any NaN values 

    HbootDistrib_bitspk{iNex,1} = H_DBSboots;
    
    Hboot_bitspkAv(iNex,1) = mean(H_DBSboots);
    Hboot_bitspkStdv(iNex,1) = std(H_DBSboots);
    Hdbs_bitspk_EmpPval(iNex,1) = getEmpiricalPval(H_DBSboots, H_DBSemp);
    
    
    % Get z-value versions of the above data:
%     rawH = [H_DBSemp; H_DBSboots];
%     Z = zscore(rawH);
    Hdbs_bitspkZ(iNex,1) = (H_DBSemp - mean(H_DBSboots)) / std(H_DBSboots);
    
    psth_Ct{iNex,1} = psth_DBSct;
    psth_nSpk(iNex,1) = sum(psth_DBSct);
    psthFR(iNex,1) = mean(psth_DBSspksec);
    
    Hdbs_bitsec(iNex,1) = H_DBSemp * mean(psth_DBSspksec);
    HbootDistrib_bitsec{iNex,1} = H_DBSboots * mean(psth_DBSspksec);
    Hboot_bitsecAv(iNex,1) = mean(H_DBSboots * mean(psth_DBSspksec));
    Hboot_bitsecStdv(iNex,1) = std(H_DBSboots * mean(psth_DBSspksec));
    
    Hdbs_bitpulse(iNex,1) = H_DBSemp * mean(psth_DBSspksec) / Tselect.dbsFrequency(iNex);
    Hboots = H_DBSboots * mean(psth_DBSspksec) / Tselect.dbsFrequency(iNex);
    HbootDistrib_bitpulse{iNex,1} = Hboots;
    Hboot_bitpulseAv(iNex,1) = mean(Hboots);
    Hboot_bitpulseStdv(iNex,1) = std(Hboots);
    Hboot_bitpulseLogAv(iNex,1) = mean(log(Hboots));
    Hboot_bitpulseLogStdv(iNex,1) = std(log(Hboots));
      
    
    
    
     %% Calculate PREDBS entropy, bit/spike
   
    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'Hpre';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'Hpre_%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_PRE(1)), ...
                              num2str(SPKINTERV_PRE(2)), ...
                              num2str(ORD_H), ...
                              num2str(BINS_PER_DECADE));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate entropy for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
         % Get spike times occurring within DBS interval 30-60
        refT = StimTs.DBS(1);
        
        preDbsInterval = SPKINTERV_PRE; %s, seconds to gather spikes and virtual stims before DBS onset
        % check and make sure that recording has enough pre-DBS time as requested
        % in "preDbsInterval", crop preDBS interval if necessary
        if abs(preDbsInterval(1)) > StimTs.DBS(1)
            preDbsInterval(1) = -StimTs.DBS(1);

        end
        
        [spkTimesPre] = getIntervalEvents(spkTimes, refT, preDbsInterval);

        % Define ISIs, and remove any ISIs == 0
        isiPRE = diff(spkTimesPre);
        isiPRE(isiPRE == 0) = [];
        nISI = numel(isiPRE);

        if nISI < 2
            H_PRE = NaN;

        else
            H_PRE = entropyISIdirect_multOrder(isiPRE, BINS_PER_DECADE, ORD_H);

        end
    %     entropyDirISI.dbsEmp = H_DBS;


        
        save(fullPathFn, 'H_PRE', 'isiPRE');
    
    end

    Hpre(iNex,1) = H_PRE;
    


    %% Calculate entropy for onDBS first 30 sec,  bit/spike
   
    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'Hdbs';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'Hdbs_%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD_%sboots';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_DBS_1(1)), ...
                              num2str(SPKINTERV_DBS_1(2)), ...
                              num2str(ORD_H), ...
                              num2str(BINS_PER_DECADE), ...
                              num2str(NBOOTS));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
         % Get spike times occurring within DBS interval 30-60
        refT = StimTs.DBS(1);
        [spkTimesDbs] = getIntervalEvents(spkTimes, refT, SPKINTERV_DBS_1);

        % Define ISIs, and remove any ISIs == 0
        isiDBS = diff(spkTimesDbs);
        isiDBS(isiDBS == 0) = [];
        nISI = numel(isiDBS);

        if nISI < 2
            H_DBS = NaN;

        else
            H_DBS = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);

        end
        
        % Generate bootstrapped distribution of Pre-DBS isi's with isi
        % count matching that of the DBS case above
        H_PREboots = zeros(NBOOTS, ORD_H);
        
        if (nISI < 2) || isempty(isiPRE)
            H_PREboots(:) = NaN;
            
        else
            for iBoot = 1:NBOOTS
                % Resample PRE ISIs (with replacement) so that number of ISIs
                % matches the DBS case
                [isiPREresamp] = datasample(isiPRE, numel(isiDBS));

                % Calculate multi-order Direct-Entropy
                H_PREboots(iBoot,:) = entropyISIdirect_multOrder(isiPREresamp, BINS_PER_DECADE, ORD_H);

            end

        end
        
        save(fullPathFn, 'H_DBS', 'H_PREboots', 'isiPRE', 'isiDBS');
    
    end

    H_PREboots(isnan(H_PREboots)) = []; % remove any NaN values 
    
    Hdbs1_preBootAv(iNex,1) = mean(H_PREboots);
              Hdbs1(iNex,1) = H_DBS;
      Hdbs1_EmpPval(iNex,1) = getEmpiricalPval(H_PREboots, H_DBS);
    
   
      
    %% Calculate entropy for onDBS second 30 sec, normalized by bit/spike
   
    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'Hdbs';
    if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
        mkdir([scriptIntermediateDataFolder, '\', subFolder]);
        
    end
    formatSpecFn = 'Hdbs_%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD_%sboots';
    matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
                              num2str(SPKINTERV_DBS_2(1)), ...
                              num2str(SPKINTERV_DBS_2(2)), ...
                              num2str(ORD_H), ...
                              num2str(BINS_PER_DECADE), ...
                              num2str(NBOOTS));
    fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
         % Get spike times occurring within DBS interval 30-60
        refT = StimTs.DBS(1);
        [spkTimesDbs] = getIntervalEvents(spkTimes, refT, SPKINTERV_DBS_2);

        % Define ISIs, and remove any ISIs == 0
        isiDBS = diff(spkTimesDbs);
        isiDBS(isiDBS == 0) = [];
        nISI = numel(isiDBS);

        if nISI < 2
            H_DBS = NaN;

        else
            H_DBS = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);

        end
        
                
        % Generate bootstrapped distribution of Pre-DBS isi's with isi
        % count matching that of the DBS case above
        H_PREboots = zeros(NBOOTS, ORD_H);

        if (nISI < 2) || isempty(isiPRE)
            H_PREboots(:) = NaN;
            
        else
            for iBoot = 1:NBOOTS
                % Resample PRE ISIs (with replacement) so that number of ISIs
                % matches the DBS case
                [isiPREresamp] = datasample(isiPRE, numel(isiDBS));

                % Calculate multi-order Direct-Entropy
                H_PREboots(iBoot,:) = entropyISIdirect_multOrder(isiPREresamp, BINS_PER_DECADE, ORD_H);

            end

        end
              
        save(fullPathFn, 'H_DBS', 'H_PREboots', 'isiPRE', 'isiDBS');
    
    end

    
    H_PREboots(isnan(H_PREboots)) = []; % remove any NaN values 
    
    Hdbs2_preBootAv(iNex,1) = mean(H_PREboots);
              Hdbs2(iNex,1) = H_DBS;
      Hdbs2_EmpPval(iNex,1) = getEmpiricalPval(H_PREboots, H_DBS);

    
    
end
toc



%% Prepare data for final visualization and analysis

% Append data of interest to Tselect
Tfinal = [Tselect, ...
     table(FRpreCorr),    table(FRdbsCorr1),    table(FRdbsCorr2),    table(FRposCorr), ...
     table(Hdbs_bitspk), table(HbootDistrib_bitspk), table(Hboot_bitspkAv), table(Hboot_bitspkStdv), ...
     table(Hdbs_bitsec), table(HbootDistrib_bitsec), table(Hboot_bitsecAv), table(Hboot_bitsecStdv), ...
     table(Hdbs_bitpulse), table(HbootDistrib_bitpulse), table(Hboot_bitpulseAv), table(Hboot_bitpulseStdv) ...
     table(Hboot_bitpulseLogAv), table(Hboot_bitpulseLogStdv), ...
     table(Hdbs_bitspk_EmpPval), table(Hpre_bitspk), table(Hdbs_bitspkZ), table(psth_Ct), ...
     table(psth_nSpk), table(psthFR)];
% 
% Tfinal = N;
% numel(unique(Tfinal.Unit_objectID))    



% Append data of interest to Tselect
Tfinal = [Tfinal, table(Hpre), ...
          table(Hdbs2), table(Hdbs2_preBootAv), table(Hdbs2_EmpPval)];

numel(unique(Tfinal.Unit_objectID)) % display current number units       
Tfinal2 = Tfinal;

% Remove any rows that contain a NaN for entropy
Tfinal2(isnan(Tfinal2.Hpre),:) = [];
% Tfinal2(isnan(Tfinal2.Hdbs1),:) = [];
Tfinal2(isnan(Tfinal2.Hdbs2),:) = [];
% Tfinal2(isnan(Tfinal2.Hpos),:) = [];
numel(unique(Tfinal2.Unit_objectID))


% % Remove any rows that contain a NaN for entropy
% Tfinal(isnan(Tfinal.Hpre),:) = [];
% Tfinal(isnan(Tfinal.Hdbs1),:) = [];
% Tfinal(isnan(Tfinal.Hdbs2),:) = [];
% % Tfinal(isnan(Tfinal.Hpos),:) = [];
% numel(unique(Tfinal.Unit_objectID))

% % Add in FiringRate Index changes
% FRdbs1_idx = calcFRindex(Tfinal.FRpreCorr, Tfinal.FRdbsCorr1);
% FRdbs2_idx = calcFRindex(Tfinal.FRpreCorr, Tfinal.FRdbsCorr2);
%  FRpos_idx = calcFRindex(Tfinal.FRpreCorr, Tfinal.FRposCorr);
 
% % Entropy values normalized to predbs Entropy for each trial
% Hdbs1_BitpSec_norm = Tfinal.Hdbs1_BitpSec ./ Tfinal.Hpre_BitpSec;
% Hdbs2_BitpSec_norm = Tfinal.Hdbs2_BitpSec ./ Tfinal.Hpre_BitpSec;
%  Hpos_BitpSec_norm = Tfinal.Hpos_BitpSec ./ Tfinal.Hpre_BitpSec;
%  Hpre_BitpSec_norm = Tfinal.Hpre_BitpSec ./ mean(Tfinal.Hpre_BitpSec);

% Hdbs1_BitpSPKperc = 100 * ((Tfinal.Hdbs1 - Tfinal.Hpre) ./ Tfinal.Hpre);
% Hdbs2_BitpSPKperc = 100 * ((Tfinal.Hdbs2 - Tfinal.Hpre) ./ Tfinal.Hpre);
%  Hpos_BitpSPKperc = 100 * ((Tfinal.Hpos - Tfinal.Hpre) ./ Tfinal.Hpre);
% 
% Hdbs1_BitpSECperc = 100 * ((Tfinal.Hdbs1_BitpSec - Tfinal.Hpre_BitpSec) ./ Tfinal.Hpre_BitpSec);
% Hdbs2_BitpSECperc = 100 * ((Tfinal.Hdbs2_BitpSec - Tfinal.Hpre_BitpSec) ./ Tfinal.Hpre_BitpSec);
%  Hpos_BitpSECperc = 100 * ((Tfinal.Hpos_BitpSec - Tfinal.Hpre_BitpSec) ./ Tfinal.Hpre_BitpSec);

% Tfinal2 = [Tfinal, table(Hdbs2_BitpSPKperc)];
Tfinal2 = [Tfinal];



% Remove any rows that have a firing rate of less than min firing rate:
% filter out any rows that have too-sparse cells, based on 30-second
% portions
isPreRateBelow = Tfinal2{:, 'FRpreCorr'} < SPARSE_HZ;
% isDbs1RateBelow = Tfinal2{:, 'FRdbsCorr1'} < SPARSE_HZ;
isDbs2RateBelow = Tfinal2{:, 'FRdbsCorr2'} < SPARSE_HZ;
isDbs2PSTHRateBelow = Tfinal2{:, 'psthFR'} < SPARSE_HZ;
% isRemove = isPreRateBelow | isDbs2RateBelow | isDbs2PSTHRateBelow;
isRemove = isPreRateBelow | isDbs2RateBelow;
Tfinal2(isRemove, :) = [];
numel(unique(Tfinal2.Unit_objectID))



% Remove any neurons that have less than minimum trials 
Tfinal3 = filter_neuronMinimumTrials(Tfinal2, NEURON_MIN_TRIALS);
numel(unique(Tfinal3.Unit_objectID))


T = Tfinal3;

% get log(Hbootdistrib) for all boots

% get av and stdv

% nRows = height(T);
% 
% for iRow = 1:nRows
%     
% 
% end
% 




Tanalyze = Tfinal3;



%%

% Gen figure for H bit/pulse, diff between log(Hdbs) - mean(log(Hboots))


% Gen figure for H bit/pulse, log(Hdbs) as Z-score relative to log(Hboots)

%% Gen figure for H bit/spike
f2 = figure;

freqs = Tanalyze{:, 'dbsFrequency'};
Entropy_grpLabel = cellstr(num2str(freqs));
disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);
entropyTitle = 'Entropy 30-60 sec';


ax(1) = subplot(1, 2, 1);
dispDeltaH_Fswp_Boxplot(Tanalyze.Hboot_bitspkAv, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers'); 
t = title([ppar.subjID, ':  BOOTav-dbs ', entropyTitle]);
ylabel('Bits/spike')

ax(2) = subplot(1, 2, 2);
disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);


diffHdbs = (Tanalyze.Hdbs_bitspk);
% diffHdbs = Tanalyze.Hdbs2;

dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
t = title(['REAL-dbs ', entropyTitle]);
ylabel('')
             
f2.Position = FIG_POSITION;

set(ax, 'YLim', [min(ax(1).YLim(1), ax(2).YLim(1)), max(ax(1).YLim(2), ax(2).YLim(2))]);



%% Gen figure for H bit/spike
f2 = figure;

freqs = Tanalyze{:, 'dbsFrequency'};
Entropy_grpLabel = cellstr(num2str(freqs));
% disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);
entropyTitle = 'Entropy 30-60 sec';


diffHdbs = (Tanalyze.Hdbs_bitspk - Tanalyze.Hboot_bitspkAv);

dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
t = title(['REAL-dbs ', entropyTitle]);
ylabel('\DeltaBits/spike')             
f2.Position = FIG_POSITION;

set(ax, 'YLim', [min(ax(1).YLim(1), ax(2).YLim(1)), max(ax(1).YLim(2), ax(2).YLim(2))]);

% Check for any significant difference in groups
[pAnova, tblAnova, statsAnova] = anova1(diffHdbs, Entropy_grpLabel);
figure
multcompare(statsAnova)

% Check t-test results for each individual Frequency group:

% check for allHz case;
hzLabel = str2num(cell2mat(Entropy_grpLabel));

uniqueLab = unique(hzLabel)
nLabs = length(uniqueLab);

for iLab = 1:nLabs
    isHz = hzLabel == uniqueLab(iLab);
    x = diffHdbs(isHz);
    [h, p, ci, stats] = ttest(x);
    
    statTable{iLab,1} = p;
    statTable{iLab,2} = ci(1);
    statTable{iLab,3} = ci(2);
    statTable{iLab,4} = stats.tstat;
    statTable{iLab,5} = stats.df;
    statTable{iLab,6} = stats.sd;
    
end


%% Gen figure for H bit/sec
f2 = figure;

freqs = Tanalyze{:, 'dbsFrequency'};
Entropy_grpLabel = cellstr(num2str(freqs));
disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);
entropyTitle = 'Entropy 30-60 sec';


ax(1) = subplot(1, 2, 1);
dispDeltaH_Fswp_Boxplot(Tanalyze.Hboot_bitsecAv, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers'); 
t = title([ppar.subjID, ':  BOOTav-dbs ', entropyTitle]);
ylabel('Bits/sec')


ax(2) = subplot(1, 2, 2);
diffHdbs = (Tanalyze.Hdbs_bitsec);
dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
t = title(['REAL-dbs ', entropyTitle]);
ylabel('')
    

f2.Position = FIG_POSITION;
set(ax, 'YLim', [min(ax(1).YLim(1), ax(2).YLim(1)), max(ax(1).YLim(2), ax(2).YLim(2))]);



% % Check for any significant difference in groups
% [pAnova, tblAnova, statsAnova] = anova1(diffHdbs, Entropy_grpLabel);
% figure
% multcompare(statsAnova)





%%  Gen figure for H bit/dbs pulse
f2 = figure;

freqs = Tanalyze{:, 'dbsFrequency'};
Entropy_grpLabel = cellstr(num2str(freqs));
disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);
entropyTitle = 'Entropy 30-60 sec';


ax(1) = subplot(1, 2, 1);
dispDeltaH_Fswp_Boxplot(Tanalyze.Hboot_bitpulseAv, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers'); 
t = title([ppar.subjID, ':  BOOTav-dbs ', entropyTitle]);
ylabel('Bits/pulse')


ax(2) = subplot(1, 2, 2);
diffHdbs = (Tanalyze.Hdbs_bitpulse);
dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
t = title(['REAL-dbs ', entropyTitle]);
ylabel('')
    

f2.Position = FIG_POSITION;
set(ax, 'YLim', [min(ax(1).YLim(1), ax(2).YLim(1)), max(ax(1).YLim(2), ax(2).YLim(2))]);



%%  Gen figure for H bit/dbs pulse
f2 = figure;

freqs = Tanalyze{:, 'dbsFrequency'};
Entropy_grpLabel = cellstr(num2str(freqs));
disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);
entropyTitle = 'Entropy 30-60 sec';


ax(1) = subplot(1, 2, 1);
dispDeltaH_Fswp_Boxplot(Tanalyze.Hboot_bitpulseLogAv, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers'); 
t = title([ppar.subjID, ':  BOOTav-dbs ', entropyTitle]);
ylabel('LOG Bits/pulse')


ax(2) = subplot(1, 2, 2);
diffHdbs = (log(Tanalyze.Hdbs_bitpulse));
dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
t = title(['REAL-dbs ', entropyTitle]);
ylabel('')
    

f2.Position = FIG_POSITION;
set(ax, 'YLim', [min(ax(1).YLim(1), ax(2).YLim(1)), max(ax(1).YLim(2), ax(2).YLim(2))]);



%% Gen figure for H bit/dbs pulse
f2 = figure;
ax2 = axes;

freqs = Tanalyze{:, 'dbsFrequency'};
Entropy_grpLabel = cellstr(num2str(freqs));

% disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);

% diffHdbs = Tanalyze.Hdbs1 - Tanalyze.Hpre;
% entropyTitle = 'DBS-on Entropy 0-30 sec';

% diffHdbs = (Tanalyze.Hdbs2) .* Tanalyze.psthFR .* (1 ./ Tanalyze.dbsFrequency);
diffHdbs = (log(Tanalyze.Hdbs_bitpulse) - Tanalyze.Hboot_bitpulseLogAv);

% diffHdbs = Tanalyze.Hdbs2;
entropyTitle = 'DBS-on Entropy 30-60 sec';

dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
ylabel('\Deltalog Bits/pulse')

f2.Position = FIG_POSITION;
t = title([ppar.subjID, ':  ', entropyTitle]);

% Check for any significant difference in groups
[pAnova, tblAnova, statsAnova] = anova1(diffHdbs, Entropy_grpLabel);
figure
multcompare(statsAnova)

% Check t-test results for each individual Frequency group:

% check for allHz case;
hzLabel = str2num(cell2mat(Entropy_grpLabel));

uniqueLab = unique(hzLabel)
nLabs = length(uniqueLab);

for iLab = 1:nLabs
    isHz = hzLabel == uniqueLab(iLab);
    x = diffHdbs(isHz);
    [h, p, ci, stats] = ttest(x);
    
    statTable{iLab,1} = p;
    statTable{iLab,2} = ci(1);
    statTable{iLab,3} = ci(2);
    statTable{iLab,4} = stats.tstat;
    statTable{iLab,5} = stats.df;
    statTable{iLab,6} = stats.sd;
    
end




%% Gen figure for H bit/dbs pulse
f2 = figure;
ax2 = axes;

freqs = Tanalyze{:, 'dbsFrequency'};
Entropy_grpLabel = cellstr(num2str(freqs));

disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);

% diffHdbs = Tanalyze.Hdbs1 - Tanalyze.Hpre;
% entropyTitle = 'DBS-on Entropy 0-30 sec';

% diffHdbs = (Tanalyze.Hdbs2) .* Tanalyze.psthFR .* (1 ./ Tanalyze.dbsFrequency);
diffHdbs = (log(Tanalyze.Hdbs_bitpulse) - Tanalyze.Hboot_bitpulseLogAv) ./ Tanalyze.Hboot_bitpulseLogStdv;
% diffHdbs = Tanalyze.Hdbs2;
entropyTitle = 'DBS-on Entropy 30-60 sec';

dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
ylabel('Z-score Bits/pulse')

f2.Position = FIG_POSITION;
t = title([ppar.subjID, ':  ', entropyTitle]);


HdbsZscore = diffHdbs;

%--------------------------------------------------------------------------
% section for classifying the cell responses as entrained or not entrained:


% get a mod/nonmod category going based on z-score p-value 
nRows = numel(HdbsZscore);
isMod = false(nRows, 1);
isNonMod = false(nRows, 1);

pValZscore = normcdf(HdbsZscore);


for iRow = 1:nRows
    if pValZscore(iRow) < ppar.pAlphaPSTH
        isMod(iRow) = true;
        
    else
        isNonMod(iRow) = true;
        
    end
    
end

% add results to analysis table
Tdisp = [Tanalyze, table(HdbsZscore), table(pValZscore), ...
            table(isMod), table(isNonMod)];

% display categorization data by DBS frequency

% display relative % proportions of each type of response vs. dbsFrequency
% pre-allocate matrix for counts with each frequency being a row; columns:
% wasIncrFR, noChange, wasDecrFR 
freqs = [10, 20, 30, 50, 100, 130];
nFreqs = numel(freqs);
categories = {'isMod', 'isNonMod'};
catCounts = zeros(nFreqs, numel(categories));

% Gather all counts
for iFr = 1:nFreqs
    % get subselection of table rows pertaining to frequency "iFr"
    T_byFreq = Tdisp(Tdisp.dbsFrequency == freqs(iFr),:);
    
    % categorize each trial into "incrH", "decrH", "noCh" (in that order), 
    % keep count
    
    catCounts(iFr,:) = tab_getVariableCounts(T_byFreq, categories);
  
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
lgd = legend('isMod', 'noChange', 'Location', 'northeastoutside') ;
tit = title([ppar.subjID, ':  population percentage by 30-60s DBS Entropy change, p < ', num2str(ppar.pAlphaPSTH)]);
ax.YLim = [0, 100];


[xtab, chi2, p, labels] = crosstab(Tdisp.dbsFrequency, Tdisp.isMod)



% For ISI entropy get categories based on deltaEntropy measure:


dbsPortionStr = 'Hdbs2';

% % change significance
% SIG_PVAL = 0.05 / 2; % to account for two-tailed distributions

% Tdisp = Tanalyze; % just to make typing easier...
nRows = height(Tdisp);

     wasDecrH = false(nRows, 1);
     wasIncrH = false(nRows, 1);
signifChangeH = false(nRows, 1);

for iRow = 1:nRows
    if Tdisp.Hdbs2(iRow) <= Tdisp.Hdbs2_preBootAv(iRow) % entorpy decrease
        wasDecrH(iRow) = true; 
        
    else % entropy increase
        wasIncrH(iRow) = true;
        
    end
    
    if abs(Tdisp.Hdbs2_EmpPval(iRow)) < SIG_PVAL
        signifChangeH(iRow) = true;
        
    end
    
end

Tdisp = [Tdisp, table(wasDecrH), table(wasIncrH), table(signifChangeH)];

% Gather all counts
freqs = [10, 20, 30, 50, 100, 130];
nFreqs = numel(freqs);
catCountsH = zeros(nFreqs, 3);
for iFr = 1:nFreqs
    % get subselection of table rows pertaining to frequency "iFr"
    T_byFreq = Tdisp(Tdisp.dbsFrequency == freqs(iFr),:);
    
    % categorize each trial into "incrH", "decrH", "noCh" (in that order), 
    % keep count
    catCountsH(iFr,:) = getCategoryCounts_table(T_byFreq);
  
end

% Hdbs2: Display the group proportions according to frequency
% pre-allocate matrix for counts with each frequency being a row; columns:
% wasDecrH, wasIncrH, noCh


%% Assign Category label to each row according to the TF values of:

% Columns correspond to whether or not: 
% 1) Phase-locked, 2) isiH increased, 3) isiH decreased
% Rows corresponed to pieLabels below, in order

% cType = [false, false, false;
%          false, false, true;
%          true, false, true;
%          true, false, false;
%          true, true, false;
%          false, true, false;];
% 
% % I've observed 6 possible categories in this data
% pieLabels = {'noChange', 'noPhsLck_Hdecr', 'noPhsLck_Hincr', ...
%              'PhsLckONLY', 'PhsLck_Hdecr', 'PhsLck_Hincr'};
         
cType = [false, false, false;
         false, false, true;
         true, true, false;
         true, false, false;
         true, false, true;
         false, true, false;];

% I've observed 6 possible categories in this data
pieLabels = {'noChange', 'noPhsLck_Hdecr', 'PhsLck_Hdecr', ...
             'PhsLckONLY', 'PhsLck_Hincr', 'noPhsLck_Hincr'};

nRows = size(Tdisp, 1);  
pieCategory = cell(nRows, 1);
for iRow = 1:nRows
    if Tdisp.isMod(iRow) && Tdisp.signifChangeH(iRow) && Tdisp.wasIncrH(iRow) 
        pieCategory{iRow} = 'PhsLck_Hincr';
        
    elseif Tdisp.isMod(iRow) && Tdisp.signifChangeH(iRow) && Tdisp.wasDecrH(iRow) 
        pieCategory{iRow} = 'PhsLck_Hdecr';
        
    elseif Tdisp.isMod(iRow) && ~Tdisp.signifChangeH(iRow)
        pieCategory{iRow} = 'PhsLckONLY';
        
    elseif ~Tdisp.isMod(iRow) && Tdisp.signifChangeH(iRow) && Tdisp.wasIncrH(iRow) 
        pieCategory{iRow} = 'noPhsLck_Hincr';
        
    elseif ~Tdisp.isMod(iRow) && Tdisp.signifChangeH(iRow) && Tdisp.wasDecrH(iRow) 
        pieCategory{iRow} = 'noPhsLck_Hdecr';
        
    else
        pieCategory{iRow} = 'noChange';
        
    end
        
end   
        
 

Tdisp = [Tdisp, table(pieCategory)];


%% TALLY all Trial Observation counts for each category below
% 

% how many groups -- note that '' is the default empty case, to be ignored
% grpLabels = unique(Tselect.group);
% grpLabels(strcmp(grpLabels, '')) = [];
% 
% 
% if isfield(ppar, 'groups')
%     grpLabels = ppar.groups.dbsFrequency.groupLabels;
%     
% else
%     grpLabels = unique(Tselect.group);
%     
% end
% 
%      
% % Gather up data for first group first, then next, etc...
% nFreqs = numel(grpLabels);
% for iGrp = 1:nFreqs
%     isCurrGrp = Tdisp.dbsFrequency, grpLabels{iGrp});
%     TbyGroup = Tselect(isCurrGrp, :);
%     
%     nObs = size(TbyGroup, 1);
%         
%     nPieLabels = numel(pieLabels);
%     for iPiLab = 1:nPieLabels % count num obs for each category
%         pieCatCount{iGrp}(iPiLab) = sum(strcmp(TbyGroup.pieCategory, pieLabels(iPiLab)));
%         pieData{iGrp}(iPiLab) = round(100 * pieCatCount{iGrp}(iPiLab) / nObs);
% 
%     end
%          
% end
% 


%%
freqs = [10, 20, 30, 50, 100, 130];
nFreqs = numel(freqs);
catCounts = zeros(nFreqs, 3);

% Gather all counts
for iFr = 1:nFreqs
    % get subselection of table rows pertaining to frequency "iFr"
    T_byFreq = Tdisp(Tdisp.dbsFrequency == freqs(iFr),:);
    
    
    nObs = size(T_byFreq, 1);
        
    nPieLabels = numel(pieLabels);
    for iPiLab = 1:nPieLabels % count num obs for each category
        pieCatCount{iFr}(iPiLab) = sum(strcmp(T_byFreq.pieCategory, pieLabels(iPiLab)));
%         pieData{iFr}(iPiLab) = round(100 * pieCatCount{iFr}(iPiLab) / nObs);
        pieData{iFr}(iPiLab) = pieCatCount{iFr}(iPiLab);

    end
    
    
    % categorize each trial into "incrH", "decrH", "noCh" (in that order), 
    % keep count
%     catCounts(iFr,:) = getCategoryCounts_table(T_byFreq);
  
end
 

%%

f1 = figure; 
% ax1 = axes;


nFreqs = numel(freqs);
if nFreqs > 1
    
    for iGrp = 1:nFreqs
        ax(iGrp) = subplot(1, nFreqs, iGrp);
        pie1 = pie(pieData{iGrp}, pieLabels);
        title([ppar.subjID, ':   ', ppar.trialType, ', ', num2str(freqs(iGrp))])

    end

else % case where there is only one group, or default case of no group
    pie1 = pie(pieData{1}, pieLabels);
    
end
    




%%

%--------------------------------------------------------------------------

% 
% % Gen figure for Z-score Hdbs (normalized difference from uniform distrib
% % PSTH)
% f2 = figure;
% ax2 = axes;
% 
% freqs = Tanalyze{:, 'dbsFrequency'};
% Entropy_grpLabel = cellstr(num2str(freqs));
% 
% disp_NpointsVsFreq(Entropy_grpLabel)
% Flabels_AllNeu = unique(Entropy_grpLabel);
% 
% % diffHdbs = Tanalyze.Hdbs1 - Tanalyze.Hpre;
% % entropyTitle = 'DBS-on Entropy 0-30 sec';
% 
% diffHdbs = Tanalyze.Hdbs2_z;
% % diffHdbs = Tanalyze.Hdbs2;
% entropyTitle = 'DBS-on Entropy 30-60 sec';
% 
% dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
%                         Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
% ylabel('Normalized DBS-on Entropy (Z-scores)')
% 
% f2.Position = FIG_POSITION;
% t = title([ppar.subjID, ':  ', entropyTitle]);
% % ax2.YLim = [-0.5, 0.5];
% 


%%




% % Check the pre-DBS case as well
% entropyTitle = 'Pre-DBS Entropy 30-60 sec';
% 
% dispDeltaH_Fswp_Boxplot(Tanalyze.Hpre, Entropy_grpLabel, ...
%                         Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
% ylabel('Bits/spike')
% 
% f2.Position = FIG_POSITION;
% t = title([ppar.subjID, ':  ', entropyTitle]);
% % ax2.YLim = [-0.5, 0.5];
% 
% disp('STOP script here')
% 
% % Check for any significant difference in groups
% [pAnova, tblAnova, statsAnova] = anova1(diffHdbs, Entropy_grpLabel);
% figure
% multcompare(statsAnova);





% Check t-test results for each individual Frequency group:

% check for allHz case;
hzLabel = str2num(cell2mat(Entropy_grpLabel));

uniqueLab = unique(hzLabel)
nLabs = length(uniqueLab);

for iLab = 1:nLabs
    isHz = hzLabel == uniqueLab(iLab);
    x = diffHdbs(isHz);
    [h, p, ci, stats] = ttest(x);
    
    statTable{iLab,1} = p;
    statTable{iLab,2} = ci(1);
    statTable{iLab,3} = ci(2);
    statTable{iLab,4} = stats.tstat;
    statTable{iLab,5} = stats.df;
    statTable{iLab,6} = stats.sd;
    
end




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

function [n] = calcFiringRates(spkTimes, StimTs, spkTimeInterval)
% take the "spkTimes" and "StimTs" data from the nexFile and calculate the
% observed spike firing rate and other information occurring within the 
% time interval "spkTimeInterval" (1x2 array, seconds, refenced to DBS 
% onset time).

% CONSTANTS
% assumed to be true for both Uva and Kramer
periStimBlank = 1 / 1000; %s, blanking time around each stim pulse


% CODE
refT = StimTs.DBS(1);

% Calculate observed spike rate for interval
[spkTimesDbs] = getIntervalEvents(spkTimes, refT, spkTimeInterval);
totIntervTime = spkTimeInterval(2) - spkTimeInterval(1);
frDbsObs = numel(spkTimesDbs) / totIntervTime;


% Estimate the "true" firing rate by correcting for DBS artifact blank
% time, as done in Moran et al 2011
stmTimesDbs = getIntervalEvents(StimTs.DBS, refT, spkTimeInterval);
totBlankTime = periStimBlank * numel(stmTimesDbs);
frDbsCorr = frDbsObs * (totIntervTime / (totIntervTime - totBlankTime));


n.refT = refT;
n.spkTimesDbs = spkTimesDbs;
n.stmTimesDbs = stmTimesDbs;
n.rateObserved = frDbsObs; % spikes/second
n.rateCorrected = frDbsCorr; % spikes/second
n.totIntervTime = totIntervTime;
n.totBlankTime = totBlankTime;

end


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








