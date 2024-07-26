function runEvokedSpikeLatency(pipeParams)
% Run  Spike-latency analysis to determind if PSTH is antidromic
% excitation, orthodromic excitation, or some kind of phase-locked
% inhibition
%
% This pipeline assumes that the following tables exist and are correct:
%
% NEXprocfiles_subjID.mat
% SortedUnits_sibjID.mat
% SweepAnalysisTrials4Paper2018_subjID.mat
%
% where "subjID" is the nhp name (i.e. 'Uva')

% TO-DO




%% DEFAULT PARAMETERS

DEFAULT.subjID        = 'XXX';
DEFAULT.tablepn       = '\';
DEFAULT.neuTypeFilter = 'SU';
DEFAULT.predbsTime    = 60; % seconds
DEFAULT.dbsTime       = 60; % seconds
DEFAULT.hzThresh      = 2; % Hz
% DEFAULT.ordH          = 2;
% DEFAULT.binsPD        = 20;



%% SET DEFAULT PARAMETERS IF USER HAS NOT SET THEM

% Check for specific User-defined pipeline parameter inputs:
if ~isfield(pipeParams, 'subjID'), pipeParams.subjID               = DEFAULT.subjID; end
if ~isfield(pipeParams, 'tablepn'), pipeParams.tablepn             = DEFAULT.tablepn; end
if ~isfield(pipeParams, 'neuTypeFilter'), pipeParams.neuTypeFilter = DEFAULT.neuTypeFilter; end
if ~isfield(pipeParams, 'predbsTime'), pipeParams.predbsTime       = DEFAULT.predbsTime; end
if ~isfield(pipeParams, 'dbsTime'), pipeParams.dbsTime             = DEFAULT.dbsTime; end
if ~isfield(pipeParams, 'hzThresh'), pipeParams.hzThresh           = DEFAULT.hzThresh; end
% if ~isfield(pipeParams, 'ordH'), pipeParams.ordH                   = DEFAULT.ordH ; end
% if ~isfield(pipeParams, 'binsPD'), pipeParams.binsPD               = DEFAULT.binsPD; end



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



%% CREATE/LOAD a table with all data to be analyzed for entropy

nexLabel = 'NEXprocfiles';


% Create a new name for the intermediate table based on the pipeline params
name = buildDataSelectTableName(pipeParams);
IntTableName = [nexLabel, 'Selected', name];

% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createSelectionTableIfNeeded(NEX, Sort, IntTableName, pipeParams);


% LOAD the specified dataSelect Table: NEX_2analyze
load([pipeParams.tablepn, '\', IntTableName, '.mat' ]) 



%% CREATE TABLE with PSTH-Entropy estimates for both pre- and DBS periods

% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyPsthLetter';
name = buildEntropyPsthLetterTableName(pipeParams);
AnalysisTableName = [baseLabel, name];


% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createEntropyPsthLetterTableIfNeeded(NEX_2analyze, AnalysisTableName, pipeParams);


% LOAD the specified Table: H_letterResults
load([pipeParams.tablepn, '\', AnalysisTableName, '.mat']);


%% CREATE/LOAD a table with all DBS trials and their stim-evoked spk latencies

% Make name for analysis table depending on pipeParams
baseLabel = 'EvokedSpikeLatency';
name = buildEvokedSpikeLatencyTableName(pipeParams);
LatencyTableName = [baseLabel, name];

% SPECIFY a subsection of the H_letterResults table that satisfy 
% p-value < pipeParams.pValAlpha
isSig = H_letterResults.pVal(:) < pipeParams.pValAlpha;
NEX_EntropyLetterSignif = H_letterResults(isSig,:);


% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createEvokedSpikeLatencyTableIfNeeded(H_letterResults, LatencyTableName, pipeParams);


% LOAD the specified Table: EvokedSpikeLatencyResults
load([pipeParams.tablepn, '\', LatencyTableName, '.mat' ]) 



%% CREATE/LOAD add to the previous table with Collision-block Analysis

% Make name for analysis table depending on pipeParams
baseLabel = 'CollisionBlock';
name = buildEvokedSpikeLatencyTableName(pipeParams);
CollisionBlockTableName = [baseLabel, name];


% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createCollisionBlockTableIfNeeded(EvokedSpikeLatencyResults, CollisionBlockTableName, pipeParams);


% LOAD the specified Table: CollisionBlockResults
load([pipeParams.tablepn, '\', CollisionBlockTableName, '.mat' ]) 




%% VIEW all mean and std stim-evoked spike-latencies

% look only at those recordings proven to be phase-locked
isSig = CollisionBlockResults.pVal < pipeParams.pValAlpha;
CollisionBlockResultsSignif = CollisionBlockResults(isSig,:);

% E = EvokedSpikeLatencyResults;

E = CollisionBlockResultsSignif;
latMeans = E.latency_mean; 
latStds = E.latency_std;

figure; 
scatter(latMeans, latStds);
xlabel('Mean evoked spike latency (seconds)');
ylabel('Stdv evoked spike latency (seconds)');
title([pipeParams.subjID, ': Stim-Evoked first-spike Latency of all recordings']);


% Try looking at conduction-block %
tcBlockPerc = E.nMissDueToTcMinus ./ E.nTcMinus;
figure;
histogram(tcBlockPerc);


% Try looking at Spike Per Pulse
spp = E.nPhsLock ./ E.nStims;
figure;
histogram(spp);


% Try looking at SPP vs. DBS entropy
spp = E.nPhsLock ./ E.nStims;
Hdbs = E.H_DBS;
figure; 
scatter(spp, Hdbs);



nRows = size(E, 1);
HPreAv = zeros(nRows, 1);
for iRow = 1:nRows
    HPreAv(iRow,1) = mean(E.H_PREbootdistr{iRow,1});
    
end
% Try looking at SPP vs. Standard Deviation vs. Letter-Entropy
spp = E.nPhsLock ./ E.nStims;
stdv = E.latency_std;
av = E.latency_mean;
deltaH = E.H_DBS - HPreAv;
Hdbs = E.H_DBS;

figure
scatter3(spp, av, deltaH)
xlabel('spp'); 
ylabel('av');
zlabel('delta-H');

% use kmeans to define 2 clusters
data = [spp, stdv, av, deltaH];
idx = kmeans(data,2);

figure;
gscatter(spp, stdv, idx)

figure; 
gscatter(spp, deltaH, idx)

figure;
gscatter(stdv, deltaH, idx)

figure; 
gscatter(av, deltaH, idx)

isAnti = idx == 1;
AntiTable = E(isAnti,:);












%% SELECT candidate antidromics and orthodromics

threshMeanLat = 4.7 / 1000; % seconds 
threshStdLat = 2.4 / 1000; % seconds

% Select all NEX files whose phase-locked neurons meet the above thresholds
isCandidate = (E.latency_mean(:) < threshMeanLat) & ... 
              (E.latency_std < threshStdLat);
          
CandidateNEX = E(isCandidate,:);
C = CandidateNEX;

latMeans = C.latency_mean(:); 
latStds = C.latency_std(:);

figure; 
scatter(latMeans, latStds);
xlabel('Mean evoked spike latency (seconds)');
ylabel('Stdv evoked spike latency (seconds)');
title([pipeParams.subjID, ': Stim-Evoked first-spike Latency of all recordings']);


%% CREATE/LOAD table of selected candidates and their pulse-following results





%% VIEW  

























end

%% SUB-FUNCTIONS

