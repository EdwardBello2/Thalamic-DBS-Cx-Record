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

specID = 'nr17080401-ch12a';



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


    
iNeu = find(strcmp(specID, uNeurons.Unit_objectID));


% FOR EACH NEURON
% for iNeu = 1:nNeus % MAIN FOR LOOP
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
% 
%     pause()
%     close(f1)
    
% end % MAIN FOR LOOP

disp('ALL DONE VIEWING NEURONS!');
