% script for looking at spike sort quality for each unit as a
% pre-processing step before final analysis




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

% Spike ISI refractory period: 
REFRAC = 1 / 1000; % seconds


% Filter the data:
NEURON_MIN_TRIALS = 3;
SPARSE_HZ = 2;


% Define time ranges of various peri-DBS events: (in seconds, make sure 2nd
% element > 1st element
  SPKINTERV_PRE = [-Inf, 0]; % [-30, 0] | [-Inf, 0]
SPKINTERV_DBS_1 = [0, 30];   % [0, 30]
SPKINTERV_DBS_2 = [30, 60];  % [30, 60]
  SPKINTERV_POS = [60, Inf]; % [60, 90] | [60, Inf]

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

%% String together all isi's for a single unit by grabbing them from
% relevant nexFiles for that Unit:
T = Tselect;

% find all unique units
unitIDs = unique(T.Unit_objectID);
nUnits = numel(unitIDs);
unitISIs = cell(nUnits, 1);

% get combined isi for each unit
for iU = 1:nUnits
    % get chunk of table rows pertaining to iU
    is_iU = strcmp(unitIDs{iU}, T.Unit_objectID);
    TiU = T(is_iU, :);
    
    % load each nexFile and extract ISI, add to current ISI tally for iU
    nNex = height(TiU);
    isi = [];
    for iNex = 1:nNex
        iNexObjectID = TiU.objectID{iNex,1};
               nexfn = TiU.nexFile{iNex,1};
               nexpn = TiU.nexFileFolder{iNex,1};
             nexFile = readNexFile([ppar.projRootPath, '\', nexpn, '\', nexfn]);

        % look at neuron timestamps
        spkTimes = nexFile.neurons{1,1}.timestamps; % seconds
        isi = [isi; diff(spkTimes)];
        
    end
    isi(isi==0) = []; % in case for some reason two timestamps are identical

    unitISIs{iU,1} = isi;
    
end



%% Loop thru each unit ISI to determine % that is less than refractory period

percInRefrac = zeros(nUnits, 1);
histEdges = 0:(REFRAC/10):2;
for iU = 1:nUnits
    percInRefrac(iU) = 100 * sum(unitISIs{iU,1} <= REFRAC) / numel(unitISIs{iU,1});
    
%     f1 = figure; histogram(unitISIs{iU,1}, histEdges);
%     title(['%inRefrac: ', num2str(percInRefrac(iU))])
%         close(f1)

end

isBadUnit = percInRefrac >= 1;





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
plot(ax, [xPoint1, xPoint1], ax.YLim, 'Color', 'r', 'LineWidth', 1.0);
plot(ax, [xPoint2, xPoint2], ax.YLim, 'Color', 'r', 'LineWidth', 1.0);

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









