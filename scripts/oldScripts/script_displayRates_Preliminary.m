% script for visualizing all potential files for rate analysis

% open questions for myself: 
% what criteria will exclude some from analysis?
% should I just match the values used in Entropy analysis? 

%%  READ in all NEXfiles for a given subject

clear; 

pipeParams.subjID = 'Uva';

% If a given intermediate table already exists, overwrite it anyway
% WARNING, may take a long time to redo some tables, hours
pipeParams.overwriteTables = false; % true/false

pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 15; % num of bins per decade of log-spaced bins for ISI

% for each neuron observed, exclude it from analysis if it does not stick
% around for at least "neuronMininumTrials" of dbs-trials (i.e. 4 trials)
pipeParams.neuronMinimumTrials = 4; 

pipeParams.entropyType = '%deltaH'; % choose: 'dbsH', 'preH', 'deltaH', 'pctDeltaH'
pipeParams.trialType = 'FrequencySweep'; % choose: 'ContactSweep', 'FrequencySweep'

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = false;



%% LOAD NECESSARY METADATA TABLES
% Table names shortened for each of use in code

% load NEXprocfiles_XXX table, where metadata for individual NEX files is
% stored:
load([pipeParams.tablepn, '\', 'NEXprocfiles_', pipeParams.subjID, '.mat']);
NEX = NEXprocfiles; 

% load Sorted_Units_XXX table, where info related to Spike-Sorting is stored:
load([pipeParams.tablepn, '\', 'SortedUnits_', pipeParams.subjID, '.mat']);
Sort = SortedUnits;

% load SweepAnalysisTrials4Paper2018_XXX table, where info on DBS
% parameters is stored
load([pipeParams.tablepn, '\', 'SweepAnalysisTrials4Paper2018_', pipeParams.subjID, '.mat']);
TrialInfo = SweepAnalysisTrials4Paper2018;



%% SELECT only those Nexfiles that pertain to Single Units

% Load Sorted Units table
isSU = strcmp(Sort.NeuronType, 'SU');
SortSU = Sort(isSU,:);

% MERGE the tables where their primary keys are common
isNexSU = ismember(NEX.Unit_objectID, SortSU.objectID);

NEXsu = NEX(isNexSU,:);





%% LOAD and display each one at a time




% Include code for manipulating nexFiles:
addpath(genpath('C:\Users\bello043\OneDrive\Tools\NEXtools'));

% % Load one NEX file, with has ONE DBS trial of data for ONE unit
% pn = 'D:\PROJECTS\Thalamic DBS Cx Record\DataProcessing\NEX_proc';
% nexfn = '17080401block01-ch01a.nex';


nNex = size(NEXsu, 1);

for iNex = 1:nNex
    close all
    
    nexFile = readNexFile([NEXsu.Pathname{iNex}, '\', NEXsu.Filename{iNex}]);



    % Get the spike times and DBS times
    [spkTimes, stimTimes] = parseNexFile(nexFile);
    dbsTimes = stimTimes.DBS;


    % create all rate bins for entire recording
    tBeg = nexFile.tbeg;
    tEnd = nexFile.tend;
    rateWin = 1;

    binEdges = tBeg:rateWin:tEnd;


    % get rate estimates for each bin by counting spikes that fall within each
    nBins = numel(binEdges) - 1;
    spkCounts = zeros(nBins, 1);

    for i = 1:nBins
        isIn_iBin = spkTimes >= binEdges(i) & ...
                    spkTimes < binEdges(i+1);

        spkCounts(i) = sum(isIn_iBin);    

    end
    spkRates = spkCounts / rateWin;
    
    


    % plot rate estimate with markers for DBS onset and offset
    t = 1:nBins;
    f1 = figure(1); ax = axes;
    plot(t, spkRates);
    title([NEXsu.Filename{iNex},', fileN ', num2str(iNex)]);
    hold on
    plot([dbsTimes(1), dbsTimes(1)], [ax.YLim(1), ax.YLim(2)])
    plot([dbsTimes(end), dbsTimes(end)], [ax.YLim(1), ax.YLim(2)])
    f1.Position = [680 558 560 420];
    
    % Plot histograms of pre and dbs rates for 60 seconds
    nBins = 60;
    rateWin = 1; % seconds

    ratesPRE = getRateBinsPreDBS(spkTimes, dbsTimes(1), 1, 60);
    ratesPRElog = log(ratesPRE + 1);

    ratesDBS = getRateBinsDurDBS(spkTimes, dbsTimes(1), 1, 60);
    ratesDBSlog = log(ratesDBS + 1);


    f2 = figure(2); histogram(ratesPRElog); hold on
    histogram(ratesDBSlog);
    legend('PRE', 'DBS'); title('Histograms of rates');

%     figure; plot(ratesPRElog); hold on;
%     plot(ratesDBSlog);
%     legend('PRE', 'DBS'); title('Rates vs. time')


    [p, h, stats] = ranksum(ratesPRElog, ratesDBSlog)
    f2.Position = [1258 560 560 420];
    
    % Decide if it's DBS-induced inhibition, excitation, or none
    if p < pipeParams.pValAlpha
        if mean(ratesDBSlog) <= mean(ratesPRElog)
            responseType = 'inhib';
            
        else
            responseType = 'excite';
            
        end
        
    else
        responseType = 'none';
        
    end
    title(['p = ', num2str(p),', ', responseType]);



    pause;
    
end


