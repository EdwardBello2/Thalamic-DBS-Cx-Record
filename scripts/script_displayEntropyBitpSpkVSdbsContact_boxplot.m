% script for looking delta Entropy vs DBS frequency, showing continuous
% line for each individual neuron
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
ppar.trialType = 'ContactSweep'; % 'FrequencySweep' | 'ContactSweep'
% ppar.groups.dbsFrequency.groupFreqs = {10, 20, 30, 50, 100, 130};
% ppar.groups.dbsFrequency.groupLabels = {'hz10', 'hz20', 'hz30', 'hz50', 'hz100', 'hz130'};



% CONSTANTS

% Spike-sort quality:
REFRAC = 1 / 1000; % seconds, neuron refractory period in which spkes shouldn't occur
PERC_REFRACTHRESH = 1.0; % threshold for rejecting a unit with x% spikes within the refractory period


% Filter the data:
NEURON_MIN_TRIALS = 4;
SPARSE_HZ = 2;



% Define time ranges of various peri-DBS events:
SPKINTERV_PRE = [-30, 0]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS_1 = [0, 30]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS_2 = [30, 60]; %seconds, make sure 2nd element > 1st
SPKINTERV_POS = [60, 90]; %seconds, make sure 2nd element > 1st

% Entropy calculation:
BINS_PER_DECADE = 15;
ORD_H = 1;

ppar.BINS_PER_DECADE = BINS_PER_DECADE;
ppar.ORD_H = ORD_H;

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

%% Get subselection of tableRoot for analysis

% Tselect = filterTable_Master(tableRoot, scriptName, ppar);
tab = tableRoot;

% ppar.neuronType:
tab = filter_neuronType(tab, ppar.neuronType);

% ppar.subjID:
tab = filter_subjID(tab, ppar.subjID);

% ppar.trialType:
tab = filter_trialType(tab, ppar.trialType);


% Remove any units from processing that do not have good-enough sort
% quality:
Tselect = removePoorQualityUnits(tab, REFRAC, PERC_REFRACTHRESH, ppar);



%% Run analysis pipeline on each nexfile specified by table above

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

    

    %% Calculate firing rates for portions of nexfile

    % Get firing rate calculations for various intervals
    FRpre = func_calcFRnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_PRE, ppar);
    FRdbs1 = func_calcFRnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_DBS_1, ppar);
    FRdbs2 = func_calcFRnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_DBS_2, ppar);
    FRpos = func_calcFRnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_POS, ppar);


    % save info in table columns
    FRpreCorr(iNex,1) = FRpre.rateCorrected; 
    FRdbsCorr1(iNex,1) = FRdbs1.rateCorrected;
    FRdbsCorr2(iNex,1) = FRdbs2.rateCorrected;
    FRposCorr(iNex,1) = FRpos.rateCorrected;

    

    %% Calculate entrpies for portions of nexfile

    % Get entropy calculations for various intervals
    H_PRE = func_calcHnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_PRE, ppar);
    H_DBS1 = func_calcHnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_DBS_1, ppar);
    H_DBS2 = func_calcHnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_DBS_2, ppar);
    H_POS = func_calcHnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_POS, ppar);
    
    
    % save info in table columns
    Hpre(iNex,1) = H_PRE; 
    Hdbs1(iNex,1) = H_DBS1; 
    Hdbs2(iNex,1) = H_DBS2; 
    Hpos(iNex,1) = H_POS; 
    
    
    
end
toc


%% Prepare data for final visualization and analysis

% calculate BitpSecond column
Hpre_BitpSec = FRpreCorr .* Hpre;
Hdbs1_BitpSec = FRdbsCorr1 .* Hdbs1;
Hdbs2_BitpSec = FRdbsCorr2 .* Hdbs2;
Hpos_BitpSec = FRposCorr .* Hpos;



% Append data of interest to Tselect
N = [Tselect, ...
     table(FRpreCorr),    table(FRdbsCorr1),    table(FRdbsCorr2),    table(FRposCorr), ...
     table(Hpre),         table(Hdbs1),         table(Hdbs2),         table(Hpos), ...
     table(Hpre_BitpSec), table(Hdbs1_BitpSec), table(Hdbs2_BitpSec), table(Hpos_BitpSec)];

Tfinal = N;
numel(unique(Tfinal.Unit_objectID))       


% Remove any rows that contain a NaN for entropy
Tfinal(isnan(Tfinal.Hpre),:) = [];
Tfinal(isnan(Tfinal.Hdbs1),:) = [];
Tfinal(isnan(Tfinal.Hdbs2),:) = [];
% Tfinal(isnan(Tfinal.Hpos),:) = [];
numel(unique(Tfinal.Unit_objectID))

% Add in FiringRate Index changes
FRdbs1_idx = calcFRindex(Tfinal.FRpreCorr, Tfinal.FRdbsCorr1);
FRdbs2_idx = calcFRindex(Tfinal.FRpreCorr, Tfinal.FRdbsCorr2);
 FRpos_idx = calcFRindex(Tfinal.FRpreCorr, Tfinal.FRposCorr);
 
% Entropy values normalized to predbs Entropy for each trial
Hdbs1_BitpSec_norm = Tfinal.Hdbs1_BitpSec ./ Tfinal.Hpre_BitpSec;
Hdbs2_BitpSec_norm = Tfinal.Hdbs2_BitpSec ./ Tfinal.Hpre_BitpSec;
 Hpos_BitpSec_norm = Tfinal.Hpos_BitpSec ./ Tfinal.Hpre_BitpSec;
 Hpre_BitpSec_norm = Tfinal.Hpre_BitpSec ./ mean(Tfinal.Hpre_BitpSec);

Hdbs1_BitpSPKperc = 100 * ((Tfinal.Hdbs1 - Tfinal.Hpre) ./ Tfinal.Hpre);
Hdbs2_BitpSPKperc = 100 * ((Tfinal.Hdbs2 - Tfinal.Hpre) ./ Tfinal.Hpre);
 Hpos_BitpSPKperc = 100 * ((Tfinal.Hpos - Tfinal.Hpre) ./ Tfinal.Hpre);

Hdbs1_BitpSECperc = 100 * ((Tfinal.Hdbs1_BitpSec - Tfinal.Hpre_BitpSec) ./ Tfinal.Hpre_BitpSec);
Hdbs2_BitpSECperc = 100 * ((Tfinal.Hdbs2_BitpSec - Tfinal.Hpre_BitpSec) ./ Tfinal.Hpre_BitpSec);
 Hpos_BitpSECperc = 100 * ((Tfinal.Hpos_BitpSec - Tfinal.Hpre_BitpSec) ./ Tfinal.Hpre_BitpSec);

Tfinal2 = [Tfinal, ...
           table(FRdbs1_idx), table(FRdbs2_idx), table(FRpos_idx), ...
           table(Hdbs1_BitpSPKperc), table(Hdbs2_BitpSPKperc), table(Hpos_BitpSPKperc), ...
           table(Hdbs1_BitpSECperc), table(Hdbs2_BitpSECperc), table(Hpos_BitpSECperc)];



% Remove any rows that have a firing rate of less than min firing rate:
% filter out any rows that have too-sparse cells, based on 30-second
% portions
isPreRateBelow = Tfinal2{:, 'FRpreCorr'} < SPARSE_HZ;
isDbs1RateBelow = Tfinal2{:, 'FRdbsCorr1'} < SPARSE_HZ;
isDbs2RateBelow = Tfinal2{:, 'FRdbsCorr2'} < SPARSE_HZ;
isRemove = isPreRateBelow | isDbs2RateBelow;
Tfinal2(isRemove, :) = [];
numel(unique(Tfinal2.Unit_objectID))


% Remove any neurons that have less than minimum trials 
Tfinal3 = filter_neuronMinimumTrials(Tfinal2, NEURON_MIN_TRIALS);
numel(unique(Tfinal3.Unit_objectID))

Tanalyze = Tfinal3;



%%

% Gen figure for DBS on last 30 sec
f2 = figure;
ax2 = axes;

Entropy_grpLabel = Tanalyze{:, 'dbsContact'};

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
% hzLabel = str2num(cell2mat(Entropy_grpLabel));
hzLabel = Entropy_grpLabel;


uniqueLab = unique(hzLabel)
nLabs = length(uniqueLab);

for iLab = 1:nLabs
    isHz = strcmp(hzLabel, uniqueLab(iLab));
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

function tabFilt = removePoorQualityUnits(Tselect, REFRAC, PERC_REFRACTHRESH, ppar)
% Check if a given unit has good-enough sort quality to include in
% analysis, remove those that are low quality:
[unitIDs, unitISIs, percInRefrac] = func_assessUnitSpkInRefrac(Tselect, REFRAC, ppar);
isQualityUnit = percInRefrac < PERC_REFRACTHRESH;
tabUnitQ = [table(unitIDs), table(isQualityUnit)];
T = join(Tselect, tabUnitQ, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'unitIDs'); 
T(~T.isQualityUnit,:) = []; % REMOVE THOSE ROWS PERTAINING TO POOR QUALITY UNIT
tabFilt = T;

disp(['Removed ', num2str(sum(~isQualityUnit)), ' of ', ...
       num2str(numel(unitIDs)),' units due to low sort quality']);
  

end

function tabFilt = filter_neuronType(tab, neuronType)

isSU = strcmp(tab.NeuronType, neuronType);
tabFilt = tab(isSU,:);

end

function tabFilt = filter_subjID(tab, subjID)
if strcmp(subjID, 'bothSubj')
    % do nothing (keep data across both subjects)
    tabFilt = tab;
    
elseif strcmp(subjID, 'Kramer') | strcmp(subjID, 'Uva')
    isSubjID = strcmp(tab.subjID, subjID);
    tabFilt = tab(isSubjID,:);
    
else
    error('Oh snap! wrong string for ppar.subjID')

end

end

function tabFilt = filter_trialType(tab, trialType)
% assumes that table has two variables: 1) ContactSweep, 2) FrequencySweep,
% both of which are columns of 1's or 0's. Currently that is how I indicate
% whether a given RECORD represents a trial for Contact or Frequency sweep.

if strcmp(trialType, 'FrequencySweep')
    isTrialSwp = logical(tab.FrequencySweep);
    tabFilt = tab(isTrialSwp,:);

    % Remove any trials where frequency sweep was not done on electrode C0
    isC0 = strcmp(tabFilt.dbsContact, 'C0');
    tabFilt(~isC0,:) = [];
    
elseif strcmp(trialType, 'ContactSweep')
    isTrialSwp = logical(tab.ContactSweep);
    tabFilt = tab(isTrialSwp,:);

    % Remove any trials where Contact sweep was not done with 130Hz freq
    is130Hz = tabFilt.dbsFrequency == 130;
    tabFilt(~is130Hz,:) = [];

else
    error('Wrong values for ppar.TrialType');

end


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










