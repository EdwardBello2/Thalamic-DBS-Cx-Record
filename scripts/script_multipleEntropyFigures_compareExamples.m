% script for looking at example cases of cell Entropy responses 
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
binEdgesPsth2 = [0:0.1:7.5] ./ 1000; % seconds




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
SIG_PVAL = 0.05 / 2; % to account for two-tailed distributions


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

% Initialize columns to add to Tanalyze for final firingRate analysis:
 FRpreCorr = zeros(nNex, 1);
FRdbsCorr1 = zeros(nNex, 1);
FRdbsCorr2 = zeros(nNex, 1);
 FRposCorr = zeros(nNex, 1);

% Initialize columns to add to Tanalyze for final Entropy analysis in this script:
Hpre = zeros(nNex, 1);

          Hdbs1 = zeros(nNex, 1);
Hdbs1_preBootAv = zeros(nNex, 1);
  Hdbs1_EmpPval = zeros(nNex, 1);

          Hdbs2 = zeros(nNex, 1);
Hdbs2_preBootAv = zeros(nNex, 1);
  Hdbs2_EmpPval = zeros(nNex, 1);

Hpos = zeros(nNex, 1);

% Initialize columns to add to Tselect for final analysis in this script:
      Hpre_bitspk = zeros(nNex, 1);
     Hdbs_bitspk = zeros(nNex, 1);
   Hdbs_bitspkZ = zeros(nNex, 1);
  H_DBSemp = zeros(nNex, 1);

      Hboot_bitspkAv = zeros(nNex, 1);
    Hboot_bitspkStdv = zeros(nNex, 1);
Hdbs_bitspk_EmpPval = zeros(nNex, 1);

trimIdx = false(length(binEdgesPsth), 1);
trimIdx(ppar.trimPSTH) = true;
HbootDistrib_bitspk = cell(nNex, 1);
psth_Ct = cell(nNex, 1);
psth_Ct2 = cell(nNex, 1);
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
    if (0 < exist(fullPathFn, 'file')) && (~ppar.overwriteIntData)
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
    if (0 < exist(fullPathFn, 'file')) && (~ppar.overwriteIntData)
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
    if (0 < exist(fullPathFn, 'file')) && (~ppar.overwriteIntData)
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
    if (0 < exist(fullPathFn, 'file')) && (~ppar.overwriteIntData)
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


        %% Calculate PREDBS entropy, normalized by bit/spike
   
    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'HpreISI';
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
    
    % load/calculate firing rates for this nexfile
    if (0 < exist(fullPathFn, 'file')) && (~ppar.overwriteIntData)
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


        
        save(fullPathFn, 'H_PRE');
    
    end

    Hpre(iNex,1) = H_PRE;
    

        
        %% Calculate entropy for onDBS second 30 sec, normalized by bit/spike
   
    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'HdbsISI';
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
    if (0 < exist(fullPathFn, 'file')) && (~ppar.overwriteIntData)
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
      ISIs_PRE(iNex,1) = {isiPRE};
      ISIs_DBS(iNex,1) = {isiDBS};



    %% Calculate entropy for onDBS second 30 sec, by bit/spike
   
    % First check if this calculation has already been done and stored in
    % intermediate data
    subFolder = 'HdbsPSTH';
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
        
        % Calc a 2nd version of the above psth, at higher resolution for
        % display
%         binEdgesPsth2 = [0:0.1:7.5] ./ 1000; % seconds
        if isempty(spkTimesDbs)
            psth_DBSct2 = zeros(1, length(binEdgesPsth2)-1);
            
        else
            psth_DBSct2 = psth(spkTimesDbs, dbsTimes, binEdgesPsth2, 'count');
            
        end
        
        
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
             'psth_DBSprob', 'psth_DBSspksec', 'psth_DBSct', 'psth_DBSct2', ...
             'psth_PREprob', 'psth_PREspksec');
    
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
    psth_Ct2{iNex,1} = psth_DBSct2;
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
     table(Hdbs_bitspk_EmpPval), table(Hpre_bitspk), table(Hdbs_bitspkZ), table(psth_Ct), table(psth_Ct2), ...
     table(psth_nSpk), table(psthFR), table(Hpre), table(Hdbs2), table(Hdbs2_preBootAv), table(Hdbs2_EmpPval), ...
     table(ISIs_PRE), table(ISIs_DBS)];

numel(unique(Tfinal.Unit_objectID))      
Tfinal2 = Tfinal;

% Remove any rows that contain a NaN for entropy
Tfinal2(isnan(Tfinal2.Hpre),:) = [];
% Tfinal2(isnan(Tfinal2.Hdbs1),:) = [];
Tfinal2(isnan(Tfinal2.Hdbs2),:) = [];
% Tfinal2(isnan(Tfinal2.Hpos),:) = [];
numel(unique(Tfinal2.Unit_objectID))


% Remove any rows that have a firing rate of less than min firing rate:
% filter out any rows that have too-sparse cells, based on 30-second
% portions
isPreRateBelow = Tfinal2{:, 'FRpreCorr'} < SPARSE_HZ;
% isDbs1RateBelow = Tfinal2{:, 'FRdbsCorr1'} < SPARSE_HZ;
isDbs2RateBelow = Tfinal2{:, 'FRdbsCorr2'} < SPARSE_HZ;
isPsthRateBelow = Tfinal2{:, 'psthFR'} < SPARSE_HZ;
isRemove = isPreRateBelow | isDbs2RateBelow | isPsthRateBelow;
Tfinal2(isRemove, :) = [];
numel(unique(Tfinal2.Unit_objectID))

% 
% % Remove any neurons that have less than minimum trials 
% Tfinal3 = filter_neuronMinimumTrials(Tfinal2, NEURON_MIN_TRIALS);
% numel(unique(Tfinal3.Unit_objectID))

Tanalyze = Tfinal2;



%% Add in PSTH entorpy zscore and categorization 
% (based on script for PSTH entropy analysis...)

T = Tanalyze;


HpsthZscore = (log(T.Hdbs_bitpulse) - Tanalyze.Hboot_bitpulseLogAv) ./ Tanalyze.Hboot_bitpulseLogStdv;

% get a mod/nonmod category going based on z-score p-value 
nRows = numel(HpsthZscore);
isMod = false(nRows, 1);
isNonMod = false(nRows, 1);

pValZscore = normcdf(HpsthZscore);


for iRow = 1:nRows
    if pValZscore(iRow) < ppar.pAlphaPSTH
        isMod(iRow) = true;
        
    else
        isNonMod(iRow) = true;
        
    end
    
end

% add results to analysis table
Tanalyze = [T, table(HpsthZscore), table(pValZscore), ...
            table(isMod), table(isNonMod)];
        

        
%% Add in delta H and categorization label for isi entropy

T = Tanalyze;
diffHdbs = T.Hdbs2 - T.Hpre;
T = [T, table(diffHdbs)];


Tdisp = T; % just to make typing easier...
nRows = height(Tdisp);

     wasDecrH = false(nRows, 1);
     wasIncrH = false(nRows, 1);
signifChangeH = false(nRows, 1);
HchangeLabel = cell(nRows, 1);

for iRow = 1:nRows
    if Tdisp.Hdbs2(iRow) <= Tdisp.Hdbs2_preBootAv(iRow) % entorpy decrease
        wasDecrH(iRow) = true; 
        
    else % entropy increase
        wasIncrH(iRow) = true;
        
    end
    
    if abs(Tdisp.Hdbs2_EmpPval(iRow)) < SIG_PVAL
        signifChangeH(iRow) = true;
        
    end
    
    if signifChangeH(iRow) && wasDecrH(iRow)
        HchangeLabel(iRow) = {'Hdecr'};
        
    elseif signifChangeH(iRow) && wasIncrH(iRow)
        HchangeLabel(iRow) = {'Hincr'};
        
    else
        HchangeLabel(iRow) = {'nonSignif'};
        
    end
    
end

Tanalyze = [Tdisp, table(wasDecrH), table(wasIncrH), table(signifChangeH), table(HchangeLabel)];




% % display categorization data by DBS frequency
% 
% % display relative % proportions of each type of response vs. dbsFrequency
% % pre-allocate matrix for counts with each frequency being a row; columns:
% % wasIncrFR, noChange, wasDecrFR 
% freqs = [10, 20, 30, 50, 100, 130];
% nFreqs = numel(freqs);
% categories = {'isMod', 'isNonMod'};
% catCounts = zeros(nFreqs, numel(categories));
% 
% % Gather all counts
% for iFr = 1:nFreqs
%     % get subselection of table rows pertaining to frequency "iFr"
%     T_byFreq = Tdisp(Tdisp.dbsFrequency == freqs(iFr),:);
%     
%     % categorize each trial into "incrH", "decrH", "noCh" (in that order), 
%     % keep count
%     
%     catCounts(iFr,:) = tab_getVariableCounts(T_byFreq, categories);
%   
% end
%  
% % transfer to percentages
% tot = sum(catCounts, 2);
% catPercs = 100 * catCounts ./ tot;
% % catPercs = catPercs';
% 
% % display results as barplots
% figure; ax = axes;
% bar(catPercs, 'stacked')
% ax.XTickLabel = freqs;
% grid on
% xlabel('DBS frequency (Hz)')
% ylabel('% recorded cells')
% lgd = legend('isMod', 'noChange', 'Location', 'northeastoutside') ;
% tit = title([ppar.subjID, ':  population percentage by 30-60s DBS Entropy change, p < ', num2str(ppar.pAlphaPSTH)]);
% ax.YLim = [0, 100];


%%

T = Tanalyze; 

% for smoothing PSTH, use gaussian kernel:
w = gausswin(7);
w = w / sum(w);

% subselect and re-order the table according to target variable
isCateg = strcmp(T.HchangeLabel, 'nonSignif');
isPSTH = T.isMod;
isKeep = isCateg & isPSTH;

T = T(isKeep,:);

T = sortrows(T, 'dbsFrequency');


nRows = height(T);
% Generate example figures for each trial, in order
for iRow = 1:nRows
    
    disp(['Row ', num2str(iRow), ' of ', num2str(nRows)]);
    
    Fr = T.dbsFrequency(iRow);
    nStims = Fr * 30;
    f = figure; 
    f.Position = [2039 498 1478 420];
    
    
    % GET FINE-DETAIL PRE-DBS PSTH ----------------------------------------
    nexfn = T.nexFile{iRow,1};
    nexpn = T.nexFileFolder{iRow,1};
    nexFile = readNexFile([ppar.projRootPath, '\', nexpn, '\', nexfn]);
    
    
    % extract spike times and DBS stim times from nexFile struct
    [spkTimes, StimTs] = parseNexFile(nexFile);

    
    % make sure spkTimes are all of time values that are monotonically
    % increasing:
    spkTimes = sort(spkTimes);
    
    
    % Get spike times occurring within 30 sec before DBS onset
    refT = StimTs.DBS(1);
    preDbsInterval = SPKINTERV_PRE; %s, seconds to gather spikes and virtual stims before DBS onset
    % check and make sure that recording has enough pre-DBS time as requested
    % in "preDbsInterval", crop preDBS interval if necessary
    if abs(preDbsInterval(1)) > StimTs.DBS(1)
        preDbsInterval(1) = -StimTs.DBS(1);

    end
    [spkTimesPre] = getIntervalEvents(spkTimes, refT, preDbsInterval);
    
    
     % get DBS times in relevent interval
     dbsTimesVirt = getIntervalEvents(StimTs.VirtPre, refT, SPKINTERV_PRE);
     
     
     % get fine-detal PSTH for pre-DBS condition
     psthPREct2 = psth(spkTimesPre, dbsTimesVirt, binEdgesPsth2, 'count');
     
     
     % 
    
    %----------------------------------------------------------------------
    
    % plot the fine-detail PSTH
    psthDbsDisp = T.psth_Ct2{iRow};
    psthPreDisp = psthPREct2;
    
    psthDbsDisp(1:9) = 0;
    psthPreDisp(1:9) = 0;
    
    psthDbsDisp2 = psthDbsDisp / (0.0001 * nStims);
    psthPreDisp2 = psthPreDisp / (0.0001 * nStims);
    
    psthDbsSm = conv(psthDbsDisp2, w, 'same');
    psthPreSm = conv(psthPreDisp2, w, 'same');

    ax1 = subplot(1,2,1);
    plot(binEdgesPsth2(2:end)*1000, psthDbsSm);
    hold on; 
    plot(binEdgesPsth2(2:end)*1000, psthPreSm);
    
%     plot(ax1.XLim, [T.psthFR(iRow), T.psthFR(iRow)], '--');
    xlabel('Time (ms)');
    ylabel('Spikes/sec');
    
    % Show the pertinent Z-score scaled PSTH entropy of PSTH
    HdbsZ = T.HpsthZscore(iRow);
    
    title(ax1, [num2str(Fr), 'Hz DBS, bits/pulse zscore: ', num2str(HdbsZ), ', isMod: ', num2str(T.isMod(iRow))]);
    
    
    % for this nexfile, show a stairplot of isi histograms, with Entropy
    % info
    [logHistPRE, logBinEdgesPRE] = isiLogBinned(T.ISIs_PRE{iRow}, BINS_PER_DECADE);
    [xSpre, ySpre] = histcounts2stairplot(logHistPRE / sum(logHistPRE), logBinEdgesPRE);
    [logHistDBS, logBinEdgesDBS] = isiLogBinned(T.ISIs_DBS{iRow}, BINS_PER_DECADE);
    [xSdbs, ySdbs] = histcounts2stairplot(logHistDBS / sum(logHistDBS), logBinEdgesDBS);
    
    diffHdbs = T.diffHdbs(iRow);
    
    ax2 = subplot(1,2,2);
    plot(xSpre * 1000, ySpre); hold on
    plot(xSdbs * 1000, ySdbs);
    legend('Pre-DBS', 'DBS-ON', 'Location', 'northeastoutside');
    ax2.XScale = 'log';
%     ax2.XLim = [0, 1000];
    xlabel('ISI, ms')
    ylabel('Probability')
    
    title(ax2, [T.objectID{iRow}, '  \DeltaH: ', num2str(diffHdbs), ', Hchange: ', T.HchangeLabel{iRow}]);
    
    
    pause()
    close(f)

    
end

disp('done viewing all!')
pause()
close all
% Gen figure for DBS on last 30 sec
f2 = figure;
ax2 = axes;

freqs = Tanalyze{:, 'dbsFrequency'};
Entropy_grpLabel = cellstr(num2str(freqs));

disp_NpointsVsFreq(Entropy_grpLabel)
Flabels_AllNeu = unique(Entropy_grpLabel);

% diffHdbs = Tanalyze.Hdbs1 - Tanalyze.Hpre;
% entropyTitle = 'DBS-on Entropy 0-30 sec';

diffHdbs = Tanalyze.Hdbs2 - Tanalyze.Hpre;
entropyTitle = 'DBS-on Entropy 30-60 sec';

dispDeltaH_Fswp_Boxplot(diffHdbs, Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'clean');                        
ylabel('\DeltaBits/spike')

f2.Position = FIG_POSITION;
t = title([ppar.subjID, ':  ', entropyTitle]);
ax2.YLim = [-0.5, 0.5];

disp('STOP script here')

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



%% SUB-FUNCTIONS

function FRidx = calcFRindex(FRbase, FRtest)
% firing rate index. 0 means no change; 1 means cell went from zero
% baseline to "something"; -1 means cell went from some activity to
% complete inhibition.

FRidx = (FRtest - FRbase) ./ (FRbase + FRtest);

end

function [intervEvents] = getIntervalEvents(evTimes, refT, refInterval)
% returns the observed firing rate within the desired time interval.
% "refInterval" indicates the window within which to gather spike times 
% works with the values in "spkTimes" and "StimTs" as inputs

rerefTimes = evTimes - refT;

% get observed spike rate "FRobs" based on count and time duration
idxInInterv = (rerefTimes >= refInterval(1)) & (rerefTimes < refInterval(2));
rerefSubselect = rerefTimes(idxInInterv);
intervEvents = rerefSubselect + refT;



end

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

function [psthSmooth] = smoothGauss(psth, nPoints)
% apply gaussian smoothing by multiplying every bin's value with a gaussian
% window and spreading out its value to its neighboring bins, then adding
% up all the newly-spread bins

% Make psth a 1 x n vector
if size(psth, 1) > size(psth, 2), psth = psth'; end


% Define gaussian window
w = gausswin(nPoints);
w = w / sum(w);
sides = floor(nPoints/2);


% Create matrix of all psth bins spread out
len = length(psth);
psthMat = zeros(len, len);
for i = 1:len
    psthMat(i,i) = psth(i);
    
end

% Multiply each of the rows by gaussian window
psthMat(i,(i-sides):(i+sides)) = w * psth(i)

% Collapse them all for a smoothed psth



end



