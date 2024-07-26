% This script looks at all Single Unit neurons in Kramer and Uva and plots
% ISIs for pre and DBS condition, and compares the two

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
% neurons that have at least 4 frequency responses present. Also filtering
% out any rows that have an absolute change of < 5% between pre and dbs
% entropy


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
isPresentforTrials = (Nselect.neuCount >= 4);
NselectF = Nselect;
NselectF(~isPresentforTrials,:) = [];


% add in step calculating % change in entropy, just in empirical ones for
% now:

nRows = size(NselectF, 1);
percChangeH = zeros(nRows, 1);

for iRow = 1:nRows
    Hpre = NselectF.H_PREemp{iRow,1}(1,1);
    Hdbs = NselectF.H_DBS{iRow,1}(1,1);
    
    percChangeH(iRow,1) = 100 * (Hpre - Hdbs) / Hpre;
    
end
isFivePercCh = (abs(percChangeH) >= 5);
P = table(percChangeH);

% append this new column to the main table
NselectF = [NselectF, P];

% filter out any rows that don't have 5% difference in entropy:
NselectF = NselectF(isFivePercCh,:);


%% Cycle thru displaying each trial's Pre-dbs and DBS-on ISI


PREDBS_TIME = 60;
BINS_PER_DECADE = 15;

nNex = size(NselectF, 1);
for iNex = 1:nNex
    % Load each nexfile one at a time
    nexpn = NselectF.Pathname{iNex,1};
    nexfn = NselectF.Filename{iNex,1};

    nexFile = readNexFile([nexpn, '\', nexfn]);
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
    isiLogProb = 100 * isiLogHist ./ sum(isiLogHist);
    
    % convert its data to step-Histogram-data
    [xSdbs, ySdbs] = histcounts2stairplot(isiLogProb, binEdges);


    Hpre = NselectF.H_PREemp{iNex,1}(1,1);
    Hdbs = NselectF.H_DBS{iNex,1}(1,1);


    % plot both in one figure
    f = figure; ax = axes;
    f.Position = [680 535 788 443];
    plot(xSpre, ySpre); hold on
    plot(xSdbs, ySdbs);

    ax.XScale = 'log';
    xlabel('log-ISI bins (seconds)');
    ylabel('isi Probability (%)')

    title([NselectF.objectID{iNex,1}, ', ', ...
        num2str(NselectF.dbsFrequency(iNex,1)), ' Hz DBS, ', ...
        'Hpre = ', num2str(Hpre), ', ',...
        'Hdbs = ', num2str(Hdbs), ', ',...
        '%diff: ', num2str(100 * (Hdbs - Hpre) / Hpre)], ...
        'Interpreter', 'none')
    legend('PRE-dbs', 'DBS-on')


    pause()

    close(f)

end

disp('ALL DONE VIEWING NEURONS!');



% 
% % Collect individual unique Neurons table
% Unit_objectID = unique(NselectF.Unit_objectID);
% uNeurons = table(Unit_objectID);
% nNeus = size(uNeurons, 1);
% 
% % FOR EACH NEURON
% for iNeu = 1:nNeus
%     % Select a mini table selecting one Neuron at a time, and reorder
%     uNeuronSelect = uNeurons.Unit_objectID{iNeu,1};
%     isNeuronSelect = strcmp(uNeuronSelect, NselectF.Unit_objectID);
%     tabNeuronSelect = NselectF(isNeuronSelect,:);
%     tabNeuronSelect = sortrows(tabNeuronSelect, {'dbsFrequency'}); % order from 10Hz to 130Hz ascending
%     
%     
%     % Collect [xS, yS] of each frequency
%     freqs = tabNeuronSelect.dbsFrequency; % should be 10,20,30,50,100,130 hz
%     
%     nFreqs = numel(freqs);
%     for iFreq = 1:nFreqs 
%         % for one frequency response at a time...
%         nexpn = tabNeuronSelect.Pathname{iFreq,1};
%         nexfn = tabNeuronSelect.Filename{iFreq,1};
% 
%         nexFile = readNexFile([nexpn, '\', nexfn]);
%         [spkTimes, StimTs] = parseNexFile(nexFile);
% 
%         % separate the spike times into pre-DBS and DBS-on
%         dbsTimes = StimTs.DBS;
%         stimPeriod = median(diff(dbsTimes));
% 
%         % get PRE-dbs spikes
%         isPreDBS = (spkTimes >= (dbsTimes(1) - PREDBS_TIME)) & ...
%                   (spkTimes < dbsTimes(1));
%         spksPRE = spkTimes(isPreDBS);
% 
%         % Define ISIs, and remove any ISIs == 0
%         isiPRE = diff(spksPRE);
%         isiPRE(isiPRE == 0) = [];
%         
%         % log-bin the isi's, function from "Entropy" toolbox
%         [isiLogHist, binEdges, ~] = isiLogBinned(isiPRE, BINS_PER_DECADE);
%         
%         % convert its data to step-Histogram-data
%         [xS, yS] = histcounts2stairplot(isiLogHist, binEdges);
%         
%         
%         % Collect step-Histogram-data into variables for plot3 function
%         xPRE{iFreq} = xS; % rescale from ms to seconds
%         yPRE{iFreq} = ones(numel(xS), 1) + (iFreq - 1);
%         zPRE{iFreq} = yS;
%         
%         
%         
%         
%         % get DBS-on spikes
%         isDBSon = (spkTimes >= dbsTimes(1)) & ...
%                   (spkTimes < (dbsTimes(end) + stimPeriod));
%         spksDBS = spkTimes(isDBSon);
% 
%         % Define ISIs, and remove any ISIs == 0
%         isiDBS = diff(spksDBS);
%         isiDBS(isiDBS == 0) = [];
%         
%         % log-bin the isi's, function from "Entropy" toolbox
%         [isiLogHist, binEdges, ~] = isiLogBinned(isiDBS, BINS_PER_DECADE);
%         
%         % convert its data to step-Histogram-data
%         [xS, yS] = histcounts2stairplot(isiLogHist, binEdges);
%         
%         
%         % Collect step-Histogram-data into variables for plot3 function
%         xDBS{iFreq} = xS; % rescale from ms to seconds
%         yDBS{iFreq} = ones(numel(xS), 1) + (iFreq - 1);
%         zDBS{iFreq} = yS;
%         
%     end
% 
%         
%     % plot3 them together
%     f1 = figure; 
%     f1.Position = [134 558 1106 420];
%     ax1 = subplot(1,2,1);
%     for iFreq = 1:nFreqs
%         plot3(xPRE{iFreq}, yPRE{iFreq}, zPRE{iFreq});
%         hold on
% 
%     end
%     ax1.XScale = 'log';
%     xlabel('log-ISI (seconds)')
%     
%     ylabel('DBS frequency (Hz)')
%     ax1.YTick = 1:6;
%     ax1.YTickLabel = freqs;
%     
%     zlabel('isi binCounts')
%     grid on
%     title([uNeuronSelect, ' - PRE'])
%     
%     Zmax = ax1.ZLim(2);
%     
%     ax2 = subplot(1,2,2);
%     for iFreq = 1:nFreqs
%         plot3(xDBS{iFreq}, yDBS{iFreq}, zDBS{iFreq});
%         hold on
% 
%     end
%     ax2.XScale = 'log';
%     xlabel('log-ISI (seconds)')
%     
%     ylabel('DBS frequency (Hz)')
%     ax2.YTick = 1:6;
%     ax2.YTickLabel = freqs;
%     
%     zlabel('isi binCounts')
%     grid on
%     title([uNeuronSelect, ' - DBS'])
%     
%     Zmax = max(Zmax, ax2.ZLim(2));
%     
%     ax1.ZLim(2) = Zmax;
%     ax2.ZLim(2) = Zmax;
% 
%     pause()
%     close(f1)
% end

