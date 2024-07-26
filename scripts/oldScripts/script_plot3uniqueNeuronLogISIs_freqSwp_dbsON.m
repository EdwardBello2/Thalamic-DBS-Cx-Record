% This script looks at all Single Unit neurons in Kramer and Uva that had
% responses to all 6 DBS frequencies tried, and plots their ISIs BEFORE DBS
% together in plot3. Cycles thru displaying each Single Unit one at a time.

% Author: Ed Bello
% Created: 06/20/2019

%% Code

clear; 

%% Load relevant tables and join them for analysis

tablepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';


% Load and join all Uva's tables
subjID = 'Uva';

    % load NEXProcfiles table
    load([tablepn, '\', 'NEXprocfiles_', subjID], 'NEXprocfiles');
    N = NEXprocfiles;

    % load SortedUnits table
    load([tablepn, '\', 'SortedUnits_', subjID], 'SortedUnits');
    U = SortedUnits;

    % load SweepAnalysisTrials4Paper2018 table
    load([tablepn, '\', 'SweepAnalysisTrials4Paper2018_', subjID], ...
                                              'SweepAnalysisTrials4Paper2018');
    S = SweepAnalysisTrials4Paper2018;

    % % load RateBins_60sec_1secBins table
    % load([tablepn, '\', 'RateBins_60sec_1secBins_', subjID], 'RateBins');
    % R = RateBins;

    % load EntropyDirectISI_<subjID>_SU_2HzThresh_60pre_60dbs_ordH2_15bins
    load([tablepn, '\', 'EntropyDirectISI_', subjID, '_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots']);
    H = H_DirectISIResults;

    % join tables to NEXprocfiles
    N2 = join(N, U, 'LeftKeys', 2, 'RightKeys', 1);

    N3 = join(N2, S, 'LeftKeys', 3, 'RightKeys', 1);

    % N4 = join(N3, R, 'LeftKeys', 6, 'RightKeys', 1);

    NcomboUva = join(H, N3, 'LeftKeys', 1, 'RightKeys', 1, ...
        'KeepOneCopy', {'Unit_objectID', 'TrialobjectID', 'Filename', 'Pathname', 'dbsFrequency'}); % for now, just take advantage that both tables have identical Primary Keys...
    
% Load and join all Kramer's tables
subjID = 'Kramer';

    % load NEXProcfiles table
    load([tablepn, '\', 'NEXprocfiles_', subjID], 'NEXprocfiles');
    N = NEXprocfiles;

    % load SortedUnits table
    load([tablepn, '\', 'SortedUnits_', subjID], 'SortedUnits');
    U = SortedUnits;

    % load SweepAnalysisTrials4Paper2018 table
    load([tablepn, '\', 'SweepAnalysisTrials4Paper2018_', subjID], ...
                                              'SweepAnalysisTrials4Paper2018');
    S = SweepAnalysisTrials4Paper2018;

    % % load RateBins_60sec_1secBins table
    % load([tablepn, '\', 'RateBins_60sec_1secBins_', subjID], 'RateBins');
    % R = RateBins;

    % load EntropyDirectISI_<subjID>_SU_2HzThresh_60pre_60dbs_ordH2_15bins
    load([tablepn, '\', 'EntropyDirectISI_', subjID, '_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots']);
    H = H_DirectISIResults;

    % join tables to NEXprocfiles
    N2 = join(N, U, 'LeftKeys', 2, 'RightKeys', 1);

    N3 = join(N2, S, 'LeftKeys', 3, 'RightKeys', 1);

    % N4 = join(N3, R, 'LeftKeys', 6, 'RightKeys', 1);

    NcomboKramer = join(H, N3, 'LeftKeys', 1, 'RightKeys', 1, ...
        'KeepOneCopy', {'Unit_objectID', 'TrialobjectID', 'Filename', 'Pathname', 'dbsFrequency'}); % for now, just take advantage that both tables have identical Primary Keys...

    
% Join Uva and Kramer's table
% (eventually it may be useful for put a "subjID" column in the joined
% table to tell the two NHP's data apart

Ncombo = [NcomboKramer; NcomboUva];



%% Specify a selection of the above joined table for final analysis
% based on only: Single-Units spikes & FrequencySweep trials, and those
% neurons that have all 6 frequency responses present


% basically this is like a database query:
isSU = strcmp(Ncombo.NeuronType, 'SU');
isFreqSwp = Ncombo.FrequencySweep == 1;
% isNotSparse = (Ncombo.meanRatePRE >= 2) | (Ncombo.meanRateDBS >= 2);
% Nselect = Ncombo((isSU & isFreqSwp & isNotSparse), :);
Nselect = Ncombo((isSU & isFreqSwp), :);



% Remove any trials where frequency sweep was not done on electrode C0
isC0 = strcmp(Nselect.dbsContact, 'C0');
Nselect(~isC0,:) = [];


% Remove any neurons from analysis that do not have at least 4 trials
% present in the current selection table

% count the number of times each neuron shows up and store in new table
% called "NeuronCounts"
neurons = Nselect.Unit_objectID;
[uNeurons, ~, uNeuIdx] = unique(neurons);

nNeurons = length(uNeurons); % number of unique neurons
neuCount = zeros(nNeurons, 1);
for iNeu = 1:nNeurons
    neuCount(iNeu) = sum(uNeuIdx == iNeu);
    
end

NeuronCounts = [table(uNeurons), table(neuCount)];

% join this table to current selection table
Nselect = join(Nselect, NeuronCounts, 'LeftKeys', 2, 'RightKeys', 1);

% remove those rows with under 4 trials
isPresentforTrials = (Nselect.neuCount == 6);
NselectF = Nselect;
NselectF(~isPresentforTrials,:) = [];



%% Cycle thru displaying each Neuron's isi during DBS of each frequency

PREDBS_TIME = 60;
BINS_PER_DECADE = 15;



% Collect individual unique Neurons table
Unit_objectID = unique(NselectF.Unit_objectID);
uNeurons = table(Unit_objectID);
nNeus = size(uNeurons, 1);

% FOR EACH NEURON
for iNeu = 1:nNeus
    % Select a mini table selecting one Neuron at a time, and reorder
    uNeuronSelect = uNeurons.Unit_objectID{iNeu,1};
    isNeuronSelect = strcmp(uNeuronSelect, NselectF.Unit_objectID);
    tabNeuronSelect = NselectF(isNeuronSelect,:);
    tabNeuronSelect = sortrows(tabNeuronSelect, {'dbsFrequency'}); % order from 10Hz to 130Hz ascending
    
    
    % Collect [xS, yS] of each frequency
    freqs = tabNeuronSelect.dbsFrequency; % should be 10,20,30,50,100,130 hz
    
    nFreqs = numel(freqs);
    for iFreq = 1:nFreqs 
        % for one frequency response at a time...
        nexpn = tabNeuronSelect.Pathname{iFreq,1};
        nexfn = tabNeuronSelect.Filename{iFreq,1};

        nexFile = readNexFile([nexpn, '\', nexfn]);
        [spkTimes, StimTs] = parseNexFile(nexFile);

        % separate the spike times into pre-DBS and DBS-on
        dbsTimes = StimTs.DBS;
        stimPeriod = median(diff(dbsTimes));

        % get DBS-on spikes
        isDBSon = (spkTimes >= dbsTimes(1)) & ...
                  (spkTimes < (dbsTimes(end) + stimPeriod));
        spksDBS = spkTimes(isDBSon);

        % Define ISIs, and remove any ISIs == 0
        isiDBS = diff(spksDBS);
        isiDBS(isiDBS == 0) = [];
        
        % log-bin the isi's, function from "Entropy" toolbox
        [isiLogHist, binEdges, ~] = isiLogBinned(isiDBS, BINS_PER_DECADE);
        
        % convert its data to step-Histogram-data
        [xS, yS] = histcounts2stairplot(isiLogHist, binEdges);
        
        
        % Collect step-Histogram-data into variables for plot3 function
        x{iFreq} = xS; % rescale from ms to seconds
        y{iFreq} = ones(numel(xS), 1) + (iFreq - 1);
        z{iFreq} = yS;
        
    end

        
    % plot3 them together
    f1 = figure; 
    for iFreq = 1:nFreqs
        plot3(x{iFreq}, y{iFreq}, z{iFreq});
        hold on

    end
    ax = gca;
    ax.XScale = 'log';
    xlabel('log-ISI (seconds)')
    
    ylabel('DBS frequency (Hz)')
    ax.YTick = 1:6;
    ax.YTickLabel = freqs;
    
    zlabel('isi binCounts')
    grid on
    title(uNeuronSelect)

    pause()
    close(f1)
end

disp('ALL DONE VIEWING NEURONS!');
