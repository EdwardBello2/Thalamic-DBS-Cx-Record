% % This script looks at all Single Unit neurons in Kramer and Uva that had
% % responses to all 6 DBS frequencies tried, and plots their ISIs during DBS
% % together in plot3. Cycles thru displaying each Single Unit one at a time.
% 
% % Author: Ed Bello
% % Created: 06/20/2019
% 
% %% Code
% 
% clear; 
% 
% %% Load relevant tables and join them for analysis
% 
% tablepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
% 
% 
% % Load and join all Uva's tables
% subjID = 'Uva';
% 
%     % load NEXProcfiles table
%     load([tablepn, '\', 'NEXprocfiles_', subjID], 'NEXprocfiles');
%     N = NEXprocfiles;
% 
%     % load SortedUnits table
%     load([tablepn, '\', 'SortedUnits_', subjID], 'SortedUnits');
%     U = SortedUnits;
% 
%     % load SweepAnalysisTrials4Paper2018 table
%     load([tablepn, '\', 'SweepAnalysisTrials4Paper2018_', subjID], ...
%                                               'SweepAnalysisTrials4Paper2018');
%     S = SweepAnalysisTrials4Paper2018;
% 
%     % % load RateBins_60sec_1secBins table
%     % load([tablepn, '\', 'RateBins_60sec_1secBins_', subjID], 'RateBins');
%     % R = RateBins;
% 
%     % load EntropyDirectISI_<subjID>_SU_2HzThresh_60pre_60dbs_ordH2_15bins
%     load([tablepn, '\', 'EntropyDirectISI_', subjID, '_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots']);
%     H = H_DirectISIResults;
% 
%     % join tables to NEXprocfiles
%     N2 = join(N, U, 'LeftKeys', 2, 'RightKeys', 1);
% 
%     N3 = join(N2, S, 'LeftKeys', 3, 'RightKeys', 1);
% 
%     % N4 = join(N3, R, 'LeftKeys', 6, 'RightKeys', 1);
% 
%     NcomboUva = join(H, N3, 'LeftKeys', 1, 'RightKeys', 1, ...
%         'KeepOneCopy', {'Unit_objectID', 'TrialobjectID', 'Filename', 'Pathname', 'dbsFrequency'}); % for now, just take advantage that both tables have identical Primary Keys...
%     
% % Load and join all Kramer's tables
% subjID = 'Kramer';
% 
%     % load NEXProcfiles table
%     load([tablepn, '\', 'NEXprocfiles_', subjID], 'NEXprocfiles');
%     N = NEXprocfiles;
% 
%     % load SortedUnits table
%     load([tablepn, '\', 'SortedUnits_', subjID], 'SortedUnits');
%     U = SortedUnits;
% 
%     % load SweepAnalysisTrials4Paper2018 table
%     load([tablepn, '\', 'SweepAnalysisTrials4Paper2018_', subjID], ...
%                                               'SweepAnalysisTrials4Paper2018');
%     S = SweepAnalysisTrials4Paper2018;
% 
%     % % load RateBins_60sec_1secBins table
%     % load([tablepn, '\', 'RateBins_60sec_1secBins_', subjID], 'RateBins');
%     % R = RateBins;
% 
%     % load EntropyDirectISI_<subjID>_SU_2HzThresh_60pre_60dbs_ordH2_15bins
%     load([tablepn, '\', 'EntropyDirectISI_', subjID, '_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots']);
%     H = H_DirectISIResults;
% 
%     % join tables to NEXprocfiles
%     N2 = join(N, U, 'LeftKeys', 2, 'RightKeys', 1);
% 
%     N3 = join(N2, S, 'LeftKeys', 3, 'RightKeys', 1);
% 
%     % N4 = join(N3, R, 'LeftKeys', 6, 'RightKeys', 1);
% 
%     NcomboKramer = join(H, N3, 'LeftKeys', 1, 'RightKeys', 1, ...
%         'KeepOneCopy', {'Unit_objectID', 'TrialobjectID', 'Filename', 'Pathname', 'dbsFrequency'}); % for now, just take advantage that both tables have identical Primary Keys...
% 
%     
% % Join Uva and Kramer's table
% % (eventually it may be useful for put a "subjID" column in the joined
% % table to tell the two NHP's data apart
% 
% Ncombo = [NcomboKramer; NcomboUva];
% 
% 
% 
% %% Specify a selection of the above joined table for final analysis
% % based on only: Single-Units spikes & FrequencySweep trials, and those
% % neurons that have all 6 frequency responses present
% 
% 
% % basically this is like a database query:
% isSU = strcmp(Ncombo.NeuronType, 'SU');
% isFreqSwp = Ncombo.FrequencySweep == 1;
% % isNotSparse = (Ncombo.meanRatePRE >= 2) | (Ncombo.meanRateDBS >= 2);
% % Nselect = Ncombo((isSU & isFreqSwp & isNotSparse), :);
% Nselect = Ncombo((isSU & isFreqSwp), :);
% 
% 
% 
% % Remove any trials where frequency sweep was not done on electrode C0
% isC0 = strcmp(Nselect.dbsContact, 'C0');
% Nselect(~isC0,:) = [];
% 
% 
% % Remove any neurons from analysis that do not have at least 4 trials
% % present in the current selection table
% 
% % count the number of times each neuron shows up and store in new table
% % called "NeuronCounts"
% neurons = Nselect.Unit_objectID;
% [uNeurons, ~, uNeuIdx] = unique(neurons);
% 
% nNeurons = length(uNeurons); % number of unique neurons
% neuCount = zeros(nNeurons, 1);
% for iNeu = 1:nNeurons
%     neuCount(iNeu) = sum(uNeuIdx == iNeu);
%     
% end
% 
% NeuronCounts = [table(uNeurons), table(neuCount)];
% 
% % join this table to current selection table
% Nselect = join(Nselect, NeuronCounts, 'LeftKeys', 2, 'RightKeys', 1);
% 
% % remove those rows with under 4 trials
% isPresentforTrials = (Nselect.neuCount == 6);
% NselectF = Nselect;
% NselectF(~isPresentforTrials,:) = [];
% 

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
CURR_FUNC = 'script_plot3uniqueNeuronLogISIs_freqSwp_PREandDBSon'; 



%% LOAD all relevant tables and MERGE them

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



%% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);



%% Cycle thru displaying each Neuron's isi during DBS of each frequency

NselectF = Tselect;

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
        nexpn = tabNeuronSelect.nexFileFolder{iFreq,1};
        nexfn = tabNeuronSelect.nexFile{iFreq,1};

        nexFile = readNexFile([ppar.projRootPath, '\', nexpn, '\', nexfn]);
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
        isiLogProb = 100 * isiLogHist ./ sum(isiLogHist); 
        
        % convert its data to step-Histogram-data
        [xS, yS] = histcounts2stairplot(isiLogProb, binEdges);
        
        
        % Collect step-Histogram-data into variables for plot3 function
        xPRE{iFreq} = xS; % rescale from ms to seconds
        yPRE{iFreq} = ones(numel(xS), 1) + (iFreq - 1);
        zPRE{iFreq} = yS;
        
        
        
        
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
        isiLogProb = 100 * isiLogHist ./ sum(isiLogHist); 
        
        % convert its data to step-Histogram-data
        [xS, yS] = histcounts2stairplot(isiLogProb, binEdges);
        
        
        % Collect step-Histogram-data into variables for plot3 function
        xDBS{iFreq} = xS; % rescale from ms to seconds
        yDBS{iFreq} = ones(numel(xS), 1) + (iFreq - 1);
        zDBS{iFreq} = yS;
        
    end

        
    % plot3 them together
    f1 = figure; 
    f1.Position = [134 558 1106 420];
    ax1 = subplot(1,2,1);
    for iFreq = 1:nFreqs
        plot3(xPRE{iFreq}, yPRE{iFreq}, zPRE{iFreq}, 'LineWidth', 1.5);
        hold on

    end
    ax1.XScale = 'log';
    xlabel('log-ISI (seconds)')
    
    ylabel('DBS frequency (Hz)')
    ax1.YTick = 1:6;
    ax1.YTickLabel = freqs;
    
    zlabel('isi Probability (%)')
    grid on
    title([uNeuronSelect, ' - PRE'])
    
    Zmax = ax1.ZLim(2);
    
    ax2 = subplot(1,2,2);
    for iFreq = 1:nFreqs
        plot3(xDBS{iFreq}, yDBS{iFreq}, zDBS{iFreq}, 'LineWidth', 1.5);
        hold on

    end
    ax2.XScale = 'log';
    xlabel('log-ISI (seconds)')
    
    ylabel('DBS frequency (Hz)')
    ax2.YTick = 1:6;
    ax2.YTickLabel = freqs;
    
    zlabel('isi Probability (%)')
    grid on
    title([uNeuronSelect, ' - DBS'])
    
    Zmax = max(Zmax, ax2.ZLim(2));
    
    ax1.ZLim(2) = Zmax;
    ax2.ZLim(2) = Zmax;

    pause()
    close(f1)
end

disp('ALL DONE VIEWING NEURONS!');
