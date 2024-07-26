% script for looking at Entropy in bits/second


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
ppar.subjID = 'Kramer';
ppar.trialType = 'FrequencySweep';
ppar.groups.dbsFrequency.groupFreqs = {10, 20, 30, 50, 100, 130};
ppar.groups.dbsFrequency.groupLabels = {'hz10', 'hz20', 'hz30', 'hz50', 'hz100', 'hz130'};


% CONSTANTS
FIG_POSITION = [14 59 560 420];
BINS_PER_DECADE = 15;
ORD_H = 1;

SPKINTERV_PRE = [-30, 0]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS = [0, 30]; %seconds, make sure 2nd element > 1st
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

    

%     %% Get DBS-ON period spike rate, both observed and corrected for artifact
%     % blanking
% 
%     % First check if this calculation has already been done and stored in
%     % intermediate data
%     subFolder = 'FRdbs';
%     if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
%         mkdir([scriptIntermediateDataFolder, '\', subFolder]);
%         
%     end
%     formatSpecFn = 'FRdbs_%s_%ss-%ss_fromDbsOnset';
%     matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
%                               num2str(SPKINTERV_DBS(1)), ...
%                               num2str(SPKINTERV_DBS(2)));
%     fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
%     
%     % load/calculate firing rates for this nexfile
%     if exist(fullPathFn, 'file')
%         load(fullPathFn)
%         
%     else
%         
%         FRdbs = calcFiringRates(spkTimes, StimTs, SPKINTERV_DBS);
%         FRdbs.intervalLimits = SPKINTERV_DBS;
%         save(fullPathFn, 'FRdbs');
%     
%     end
%     

    
    
%     %% Get POST-DBS period spike rate, both observed and corrected for artifact
%     % blanking
% 
%     % First check if this calculation has already been done and stored in
%     % intermediate data
%     subFolder = 'FRpos';
%     if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
%         mkdir([scriptIntermediateDataFolder, '\', subFolder]);
%         
%     end
%     formatSpecFn = 'FRpos_%s_%ss-%ss_fromDbsOnset';
%     matfnStr = sprintf(formatSpecFn, Tselect.objectID{iNex}, ...
%                               num2str(SPKINTERV_POS(1)), ...
%                               num2str(SPKINTERV_POS(2)));
%     fullPathFn = [scriptIntermediateDataFolder, '\', subFolder, '\', matfnStr, '.mat'];
%     
%     % load/calculate firing rates for this nexfile
%     if exist(fullPathFn, 'file')
%         load(fullPathFn)
%         
%     else
%         posDbsInterval = SPKINTERV_POS; %s, seconds to gather spikes and virtual stims before DBS onset
%         % check and make sure that recording has enough pre-DBS time as requested
%         % in "preDbsInterval", crop preDBS interval if necessary
%         if posDbsInterval(2) > nexFile.tend
%             posDbsInterval(2) = nexFile.tend;
% 
%         end
%         
%         FRpos = calcFiringRates(spkTimes, StimTs, posDbsInterval);
%         FRpos.intervalLimits = posDbsInterval;
%         save(fullPathFn, 'FRpos');
%     
%     end
%     
    

    %% Calculate PREDBS entropy, normalized by bit/spike
   
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
    
    % load/calculate firing rates for this nexfile
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
        
        [spkTimesDbs] = getIntervalEvents(spkTimes, refT, preDbsInterval);

        % Define ISIs, and remove any ISIs == 0
        isiDBS = diff(spkTimesDbs);
        isiDBS(isiDBS == 0) = [];
        nISI = numel(isiDBS);

        if nISI < 2
            H_PRE = NaN;

        else
            H_PRE = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);

        end
    %     entropyDirISI.dbsEmp = H_DBS;


        
        save(fullPathFn, 'H_PRE');
    
    end

    Hpre(iNex,1) = H_PRE;



end
toc

%% Prepare data for final visualization and analysis

% calculate BitpSecond column
BitpSec = FRpreCorr .* Hpre;

% Append data of interest to Tselect
N = [Tselect, table(FRpreCorr), table(Hpre), table(BitpSec)];

% Remove any rows that contain a NaN for entropy
Tfinal = N;
Tfinal(isnan(Tfinal.Hpre),:) = [];



%% Display boxplots & ANOVA results of ALL cells 

f1 = figure;
ax1 = axes;

% Get Frequency labels of all neurons as array of strings
if isfield(ppar, 'groups')
    if isfield(ppar.groups, 'from')
        % get sub-selection of group labels and entropy values 
        % according to user-grouping
        idxGrp = Tfinal{:, ppar.groups.newLabel{1}};
        freqs = Tfinal{idxGrp, 'dbsFrequency'};

        % Whittle down the Entropy results to user-defined subselection
        EntropyChange_AllNeu = EntropyChange_AllNeu(idxGrp);
        Entropy_grpLabel = cell(size(EntropyChange_AllNeu, 1), 1);
        Entropy_grpLabel(:) = {ppar.groups.newLabel{1}};

    else
        freqs = Tfinal{:, 'dbsFrequency'};
        Entropy_grpLabel = cellstr(num2str(freqs));

    end

end

Flabels_AllNeu = unique(Entropy_grpLabel);

dispDeltaH_Fswp_Boxplot(Tfinal.BitpSec(:), Entropy_grpLabel, ...
                        Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
ylabel('Bits/second')

f1.Position = FIG_POSITION;
entropyTitle = 'DBS-on Entropy';
t = title([ppar.subjID, ':  ', entropyTitle]); 



%% Statistical Test on data

if ~isfield(ppar, 'statTest'), ppar.statTest = 'anova'; end

switch ppar.statTest
    
    case 'anova'
        disp('ANOVA results:')
        [pF, tabF, statsF] = anova1(Tfinal.BitpSec(:), Entropy_grpLabel)

    case 'kruskalwallis'
        disp('KRUSKAL-WALLIS results:')
        [pV, tab, stats] = kruskalwallis(Tfinal.BitpSec(:), Entropy_grpLabel)

    case 'ttest'
        disp('T-TEST results:')
        [h, pV, ci, stats] = ttest(Tfinal.BitpSec(:))
        statTable = table([stats.tstat; stats.df; stats.sd; pV]);
        statTable.Properties.RowNames = {'tstat', 'df', 'sd', 'p'};
        
    otherwise
        error('Oh snap! wrong string input for statTest')
        
end

% Final table summarizing the number of recordings for each condition:


% Final "deltaH_cell" detailing data points according to group order in "label_cell" 
labels_cell= unique(Entropy_grpLabel);
 
nLabs = numel(labels_cell);
for iLab = 1:nLabs
    isLabel = strcmp(labels_cell{iLab}, Entropy_grpLabel);
    deltaH_cell{iLab} = Tfinal.BitpSec(isLabel);
    
end





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








