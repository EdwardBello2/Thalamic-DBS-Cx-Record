% script for looking at Entropy in bits/second


% Cited:
% Moran, A., Stein, E., Tischler, H., Belelovsky, K. & Bar-Gad, 
% I. Dynamic Stereotypic Responses of Basal Ganglia Neurons to Subthalamic
% Nucleus High-Frequency Stimulation in the Parkinsonian Primate. 
% Frontiers in Systems Neuroscience 5, 21�21 (2011).



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
ppar.trialType = 'ContactSweep'; % 'FrequencySweep' | 'ContactSweep'
% ppar.groups.dbsFrequency.groupFreqs = {10, 20, 30, 50, 100, 130};
% ppar.groups.dbsFrequency.groupLabels = {'hz10', 'hz20', 'hz30', 'hz50', 'hz100', 'hz130'};

NEURON_MIN_TRIALS = 3;
SPARSE_HZ = 2;

% CONSTANTS
FIG_POSITION = [14 59 560 420];
BINS_PER_DECADE = 15;
ORD_H = 1;

SPKINTERV_PRE = [-30, 0]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS_1 = [0, 30]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS_2 = [30, 60]; %seconds, make sure 2nd element > 1st
SPKINTERV_POS = [60, 90]; %seconds, make sure 2nd element > 1st

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

nNex = height(Tselect);

% Initialize columns to add to Tselect for final analysis in this script:
FRpreCorr = zeros(nNex, 1);
Hpre = zeros(nNex, 1);
FRdbsCorr1 = zeros(nNex, 1);
Hdbs1 = zeros(nNex, 1);
FRdbsCorr2 = zeros(nNex, 1);
Hdbs2 = zeros(nNex, 1);
FRposCorr = zeros(nNex, 1);
Hpos = zeros(nNex, 1);


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


%         %% Calculate PREDBS entropy, normalized by bit/spike
%    
%     % First check if this calculation has already been done and stored in
%     % intermediate data
%     subFolder = 'Hpre';
%     if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
%         mkdir([scriptIntermediateDataFolder, '\', subFolder]);
%         
%     end
%     formatSpecFn = 'Hpre_%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD';
%     matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
%                               num2str(SPKINTERV_PRE(1)), ...
%                               num2str(SPKINTERV_PRE(2)), ...
%                               num2str(ORD_H), ...
%                               num2str(BINS_PER_DECADE));
%     fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
%     
%     % load/calculate firing rates for this nexfile
%     if exist(fullPathFn, 'file')
%         load(fullPathFn)
%         
%     else
%          % Get spike times occurring within DBS interval 30-60
%         refT = StimTs.DBS(1);
%         
%         preDbsInterval = SPKINTERV_PRE; %s, seconds to gather spikes and virtual stims before DBS onset
%         % check and make sure that recording has enough pre-DBS time as requested
%         % in "preDbsInterval", crop preDBS interval if necessary
%         if abs(preDbsInterval(1)) > StimTs.DBS(1)
%             preDbsInterval(1) = -StimTs.DBS(1);
% 
%         end
%         
%         [spkTimesDbs] = getIntervalEvents(spkTimes, refT, preDbsInterval);
% 
%         % Define ISIs, and remove any ISIs == 0
%         isiDBS = diff(spkTimesDbs);
%         isiDBS(isiDBS == 0) = [];
%         nISI = numel(isiDBS);
% 
%         if nISI < 2
%             H_PRE = NaN;
% 
%         else
%             H_PRE = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);
% 
%         end
%     %     entropyDirISI.dbsEmp = H_DBS;
% 
% 
%         
%         save(fullPathFn, 'H_PRE');
%     
%     end
% 
%     Hpre(iNex,1) = H_PRE;
%     
% 
% 
%     %% Calculate entropy for onDBS first 30 sec, normalized by bit/spike
%    
%     % First check if this calculation has already been done and stored in
%     % intermediate data
%     subFolder = 'Hdbs';
%     if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
%         mkdir([scriptIntermediateDataFolder, '\', subFolder]);
%         
%     end
%     formatSpecFn = 'Hdbs_%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD';
%     matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
%                               num2str(SPKINTERV_DBS_1(1)), ...
%                               num2str(SPKINTERV_DBS_1(2)), ...
%                               num2str(ORD_H), ...
%                               num2str(BINS_PER_DECADE));
%     fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
%     
%     % load/calculate firing rates for this nexfile
%     if exist(fullPathFn, 'file')
%         load(fullPathFn)
%         
%     else
%          % Get spike times occurring within DBS interval 30-60
%         refT = StimTs.DBS(1);
%         [spkTimesDbs] = getIntervalEvents(spkTimes, refT, SPKINTERV_DBS_1);
% 
%         % Define ISIs, and remove any ISIs == 0
%         isiDBS = diff(spkTimesDbs);
%         isiDBS(isiDBS == 0) = [];
%         nISI = numel(isiDBS);
% 
%         if nISI < 2
%             H_DBS = NaN;
% 
%         else
%             H_DBS = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);
% 
%         end
%     %     entropyDirISI.dbsEmp = H_DBS;
% 
% 
%         
%         save(fullPathFn, 'H_DBS');
%     
%     end
% 
%     Hdbs1(iNex,1) = H_DBS;
%     
%     
%     %% Calculate entropy for onDBS second 30 sec, normalized by bit/spike
%    
%     % First check if this calculation has already been done and stored in
%     % intermediate data
%     subFolder = 'Hdbs';
%     if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
%         mkdir([scriptIntermediateDataFolder, '\', subFolder]);
%         
%     end
%     formatSpecFn = 'Hdbs_%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD';
%     matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
%                               num2str(SPKINTERV_DBS_2(1)), ...
%                               num2str(SPKINTERV_DBS_2(2)), ...
%                               num2str(ORD_H), ...
%                               num2str(BINS_PER_DECADE));
%     fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
%     
%     % load/calculate firing rates for this nexfile
%     if exist(fullPathFn, 'file')
%         load(fullPathFn)
%         
%     else
%          % Get spike times occurring within DBS interval 30-60
%         refT = StimTs.DBS(1);
%         [spkTimesDbs] = getIntervalEvents(spkTimes, refT, SPKINTERV_DBS_2);
% 
%         % Define ISIs, and remove any ISIs == 0
%         isiDBS = diff(spkTimesDbs);
%         isiDBS(isiDBS == 0) = [];
%         nISI = numel(isiDBS);
% 
%         if nISI < 2
%             H_DBS = NaN;
% 
%         else
%             H_DBS = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);
% 
%         end
%     %     entropyDirISI.dbsEmp = H_DBS;
% 
% 
%         
%         save(fullPathFn, 'H_DBS');
%     
%     end
% 
%     Hdbs2(iNex,1) = H_DBS;
%     
%     
%     
%     %% Calculate entropy for Post DBS second 30 sec, normalized by bit/spike
%    
%     % First check if this calculation has already been done and stored in
%     % intermediate data
%     subFolder = 'Hpos';
%     if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
%         mkdir([scriptIntermediateDataFolder, '\', subFolder]);
%         
%     end
%     formatSpecFn = 'Hpos_%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD';
%     matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
%                               num2str(SPKINTERV_POS(1)), ...
%                               num2str(SPKINTERV_POS(2)), ...
%                               num2str(ORD_H), ...
%                               num2str(BINS_PER_DECADE));
%     fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
%     
%     % load/calculate firing rates for this nexfile
%     if exist(fullPathFn, 'file')
%         load(fullPathFn)
%         
%     else
%          % Get spike times occurring within DBS interval 30-60
%         refT = StimTs.DBS(1);
%         [spkTimesDbs] = getIntervalEvents(spkTimes, refT, SPKINTERV_POS);
% 
%         % Define ISIs, and remove any ISIs == 0
%         isiDBS = diff(spkTimesDbs);
%         isiDBS(isiDBS == 0) = [];
%         nISI = numel(isiDBS);
% 
%         if nISI < 2
%             H_DBS = NaN;
% 
%         else
%             H_DBS = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);
% 
%         end
%     %     entropyDirISI.dbsEmp = H_DBS;
% 
% 
%         
%         save(fullPathFn, 'H_DBS');
%     
%     end
% 
%     Hpos(iNex,1) = H_DBS;
%     
%     
%     
end
toc

%% Prepare data for final visualization and analysis

%% Prepare data for final visualization and analysis

% % calculate BitpSecond column
% Hpre_BitpSec = FRpreCorr .* Hpre;
% Hdbs1_BitpSec = FRdbsCorr1 .* Hdbs1;
% Hdbs2_BitpSec = FRdbsCorr2 .* Hdbs2;
% Hpos_BitpSec = FRposCorr .* Hpos;



% Append data of interest to Tselect
N = [Tselect, table(FRpreCorr), table(FRdbsCorr1), table(FRdbsCorr2), ...
     table(FRposCorr)];

Tfinal = N;
numel(unique(Tfinal.Unit_objectID))       

% 
% % Remove any rows that contain a NaN for entropy
% Tfinal(isnan(Tfinal.Hpre),:) = [];
% Tfinal(isnan(Tfinal.Hdbs1),:) = [];
% Tfinal(isnan(Tfinal.Hdbs2),:) = [];
% Tfinal(isnan(Tfinal.Hpos),:) = [];
% numel(unique(Tfinal.Unit_objectID))

% Add in FiringRate Index changes
FRdbs1_idx = calcFRindex(FRpreCorr, FRdbsCorr1);
FRdbs2_idx = calcFRindex(FRpreCorr, FRdbsCorr2);
 FRpos_idx = calcFRindex(FRpreCorr, FRposCorr);
%  
% Entropy values normalized to predbs Entropy for each trial
% Hdbs1_BitpSec_norm = Tfinal.Hdbs1_BitpSec ./ Tfinal.Hpre_BitpSec;
% Hdbs2_BitpSec_norm = Tfinal.Hdbs2_BitpSec ./ Tfinal.Hpre_BitpSec;
% Hpos_BitpSec_norm = Tfinal.Hpos_BitpSec ./ Tfinal.Hpre_BitpSec;
% Hpre_BitpSec_norm = Tfinal.Hpre_BitpSec ./ mean(Tfinal.Hpre_BitpSec);
% % 
% % Hdbs1_BitpSPKperc = 100 * ((Tfinal.Hdbs1 - Tfinal.Hpre) ./ Tfinal.Hpre);
% % Hdbs2_BitpSPKperc = 100 * ((Tfinal.Hdbs2 - Tfinal.Hpre) ./ Tfinal.Hpre);
% % Hpos_BitpSPKperc = 100 * ((Tfinal.Hpos - Tfinal.Hpre) ./ Tfinal.Hpre);
% % 
% % 
% % 
Tfinal2 = [Tfinal, table(FRdbs1_idx), table(FRdbs2_idx), table(FRpos_idx)];



% % Remove any rows that have a firing rate of less than min firing rate:
% % filter out any rows that have too-sparse cells, based on 30-second
% % portions
% isPreRateBelow = Tfinal{:, 'FRpreCorr'} < SPARSE_HZ;
% isDbsRateBelow1 = Tfinal{:, 'FRdbsCorr1'} < SPARSE_HZ;
% isDbsRateBelow2 = Tfinal{:, 'FRdbsCorr2'} < SPARSE_HZ;
% isRemove = isPreRateBelow | isDbsRateBelow1 | isDbsRateBelow2;
% Tfinal(isRemove, :) = [];
% numel(unique(Tfinal.Unit_objectID))


% % Remove any neurons that have less than minimum trials 
% Tfinal = filter_neuronMinimumTrials(Tfinal, NEURON_MIN_TRIALS);
% numel(unique(Tfinal.Unit_objectID))


% 
% 
% % Join the recording electrode table to main results, and get neuron depth
% % relative to cortical surface for each row
% load([ppar.tablePath, '\RecordingElectrodeDepths_Uva']);
% 
% Test = join(Tfinal, RecordingElectrodeDepths, 'LeftKeys', {'YY', 'MM', 'DD'}, ...
%             'RightKeys', {'YY', 'MM', 'DD'});
% 
% % Column labels pertaining to "RecordingElectrodeDepths_Uva" table
% chanLabels = {'ch01', 'ch02', 'ch03', 'ch04', 'ch05', 'ch06', 'ch07', 'ch08', ...
%               'ch09', 'ch10', 'ch11', 'ch12', 'ch13', 'ch14', 'ch15', 'ch16'};        
% nRows = height(Test);
% neuronDepth = zeros(nRows, 1);
% for iRow = 1:nRows    
%     % find which channel this recording corresponds to
%     ch = Test.ChannelRec(iRow);
%    
%     neuronDepth(iRow) = Test.CorticalSurface(iRow) - Test{iRow, chanLabels{ch}};
%    
% end
% 
% Tfinal2 = [Test, table(neuronDepth)];
% 


%% Display scatter of Depth vs. Firing Rate (corrected)
% 
% f = figure;
% f.Position = [21 563 1857 420];
% ax1 = subplot(1,4,1);
% ax2 = subplot(1,4,2);
% ax3 = subplot(1,4,3);
% ax4 = subplot(1,4,4);
% 
% scatter(ax1, Tfinal2{:, 'FRpreCorr'}, Tfinal2.neuronDepth); grid(ax1, 'on');
% scatter(ax2, Tfinal2{:, 'FRdbsCorr1'}, Tfinal2.neuronDepth); grid(ax2, 'on');
% scatter(ax3, Tfinal2{:, 'FRdbsCorr2'}, Tfinal2.neuronDepth); grid(ax3, 'on');
% scatter(ax4, Tfinal2{:, 'FRposCorr'}, Tfinal2.neuronDepth); grid(ax4, 'on');
% 
% ylabel(ax1, 'Depth from cortical Surface (mm)')
% xlabel(ax2, 'spikes/second')
% title(ax1, 'FRpre');
% title(ax2, 'FRdbs1');
% title(ax3, 'FRdbs2');
% title(ax4, 'FRpost');
% 
% % freq sweep Uva
% ax1.XLim = [0, 120];
% ax2.XLim = [0, 120];
% ax3.XLim = [0, 120];
% ax4.XLim = [0, 120];
% 
% % % contact sweep Uva
% % ax1.XLim = [0, 150];
% % ax2.XLim = [0, 150];
% % ax3.XLim = [0, 150];
% % ax4.XLim = [0, 150];
% 
% % all depths
% ax1.YLim = [-2.5, 0];
% ax2.YLim = [-2.5, 0];
% ax3.YLim = [-2.5, 0];
% ax4.YLim = [-2.5, 0];



%% Display scatter of Depth vs. delta Firing Rate (corrected)

f = figure;
f.Position = [21 563 1857 420];
ax1 = subplot(1,3,1);
ax2 = subplot(1,3,2);
ax3 = subplot(1,3,3);

% deltaFRdbsCorr1 = Tfinal2{:, 'FRdbsCorr1'} - Tfinal2{:, 'FRpreCorr'};
% deltaFRdbsCorr2 = Tfinal2{:, 'FRdbsCorr2'} - Tfinal2{:, 'FRpreCorr'};
% deltaFRposCorr = Tfinal2{:, 'FRposCorr'} - Tfinal2{:, 'FRpreCorr'};

scatter(ax1, Tfinal2.FRdbs1_idx, Tfinal2.FRpreCorr); grid(ax1, 'on');
scatter(ax2, Tfinal2.FRdbs2_idx, Tfinal2.FRpreCorr); grid(ax2, 'on');
scatter(ax3, Tfinal2.FRpos_idx, Tfinal2.FRpreCorr); grid(ax3, 'on');

ylabel(ax1, 'Pre-DBS firing rate (baseline) spk/sec')
xlabel(ax2, 'FR change index')
title(ax1, 'FRdbs1');
title(ax2, 'FRdbs2');
title(ax3, 'FRpost');

% freq sweep Uva
ax1.XLim = [-1, 1];
ax2.XLim = [-1, 1];
ax3.XLim = [-1, 1];

% % contact sweep Uva
% ax1.XLim = [-40, 140];
% ax2.XLim = [-40, 140];
% ax3.XLim = [-40, 140];

% all depths
ax1.YLim = [0, 60];
ax2.YLim = [0, 60];
ax3.YLim = [0, 60];


 
%% Display scatterHIST of Depth vs. delta Firing Rate (corrected)

f = figure;
f.Position = [21 563 1857 420];
ax1 = subplot(1,3,1);
ax2 = subplot(1,3,2);
ax3 = subplot(1,3,3);

% deltaFRdbsCorr1 = Tfinal2{:, 'FRdbsCorr1'} - Tfinal2{:, 'FRpreCorr'};
% deltaFRdbsCorr2 = Tfinal2{:, 'FRdbsCorr2'} - Tfinal2{:, 'FRpreCorr'};
% deltaFRposCorr = Tfinal2{:, 'FRposCorr'} - Tfinal2{:, 'FRpreCorr'};

scatter(ax1, Tfinal2.FRdbs1_idx, Tfinal2.FRpreCorr); grid(ax1, 'on');
scatter(ax2, Tfinal2.FRdbs2_idx, Tfinal2.FRpreCorr); grid(ax2, 'on');
scatter(ax3, Tfinal2.FRpos_idx, Tfinal2.FRpreCorr); grid(ax3, 'on');

ylabel(ax1, 'Pre-DBS firing rate (baseline) spk/sec')
xlabel(ax2, 'FR change index')
title(ax1, 'FRdbs1');
title(ax2, 'FRdbs2');
title(ax3, 'FRpost');

% freq sweep Uva
ax1.XLim = [-1, 1];
ax2.XLim = [-1, 1];
ax3.XLim = [-1, 1];

% % contact sweep Uva
% ax1.XLim = [-40, 140];
% ax2.XLim = [-40, 140];
% ax3.XLim = [-40, 140];

% all depths
ax1.YLim = [0, 60];
ax2.YLim = [0, 60];
ax3.YLim = [0, 60];



%% SUB-FUNCTIONS

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

function FRidx = calcFRindex(FRbase, FRtest)
% firing rate index. 0 means no change; 1 means cell went from zero
% baseline to "something"; -1 means cell went from some activity to
% complete inhibition.

FRidx = (FRtest - FRbase) ./ (FRbase + FRtest);

end









