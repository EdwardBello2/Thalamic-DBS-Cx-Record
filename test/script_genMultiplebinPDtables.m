% script to run a bunch of types of ISI Entropy


%% EntropyDirectISIord1_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation.

clear; 

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 10; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

selectBinsPD = [5, 10, 15];

for i = 1:3
    pipeParams.binsPD = selectBinsPD(i);
    runPipeline('EntropyDirectISIord1_byCell', pipeParams);
    close all
end


% EntropyDirectISIord1_byCell pipeline
% Shows delta-Entropies according to ISI-based (log-binned) Entropy
% estimation.

clear; 

pipeParams.subjID = 'Kramer';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 10; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
% pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = false;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190211_PipelineFigures';

selectBinsPD = [5, 10, 15];

for i = 1:3
    pipeParams.binsPD = selectBinsPD(i);
    try
        runPipeline('EntropyDirectISIord1_byCell', pipeParams);
    catch
        
    end
    
end
