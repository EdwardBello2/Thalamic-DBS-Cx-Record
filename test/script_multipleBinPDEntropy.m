%% EntropyDirect ISI pipeline

clear; 

pipeParams.subjID = 'Uva';
pipeParams.overwriteTables = false; % If a given intermediate table already exists, overwrite it anyway
pipeParams.tablepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
pipeParams.neuTypeFilter = 'SU';
pipeParams.predbsTime = 60; % seconds
pipeParams.dbsTime = 60; % seconds
pipeParams.hzThresh = 2; % Hz
pipeParams.ordH = 2; % num of Entropy orders to fit line to for extrapolation
pipeParams.binsPD = 20; % num of bins per decade of log-spaced bins for ISI

% Bootstrap PSTH parameters
pipeParams.nBoot = 10000;

% Significance
pipeParams.pValAlpha = 0.05;

% save figures
pipeParams.maxdH = -1.5;

pipeParams.finalFigs.save = true;
pipeParams.finalFigs.savepn = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\Reports\Report-190122_Group figures';


binPD = [5, 10, 30];

% for Uva
for iPD = 1:3
    pipeParams.binsPD = binPD(iPD);
    runPipeline('EntropyDirectISI', pipeParams);
    
end

% for Kramer
pipeParams.subjID = 'Uva';
for iPD = 1:3
    pipeParams.binsPD = binPD(iPD);
    runPipeline('EntropyDirectISI', pipeParams);
    
end
