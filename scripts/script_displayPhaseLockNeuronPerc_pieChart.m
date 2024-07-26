% Display the population percentage of phase-locked cells, according to
% dbs-trial condition (% at each frequency | % at each Contact)

% Author: Ed Bello
% Created: 2019/07/04

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
CURR_FUNC = 'script_displayPhaseLockNeuronPerc_pieChart'; 



%% LOAD all relevant tables and MERGE them

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



%% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);



%% Detect those neurons that phase-locked vs those that did not, label rows
% phase-locked must have a p-value less than ppar.pAlphaPSTH

% Add in phase-locked TF value for each row
nRows = size(Tselect, 1);
PhaseLocked = false(nRows, 1);

for iRec = 1:nRows
    if Tselect.pVal_Hpsth(iRec,1) <= ppar.pAlphaPSTH
        PhaseLocked(iRec,1) = true;
        
    end
    
end

Tselect = [Tselect, table(PhaseLocked)];



%% DETECT which neurons had at least one phase-lock reaction to DBS 

% Tally list of unique neurons with at least one phase-locking in new table
isPhsLck = Tselect.PhaseLocked;
P = Tselect(isPhsLck,:);
neuPhsLck = unique(P.Unit_objectID);

 isPhsLckNeuron = true(size(neuPhsLck, 1), 1);
phsLckTab = [table(neuPhsLck), table(isPhsLckNeuron)];

NphsLckInfo = outerjoin(Tselect, phsLckTab, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'neuPhsLck');

isPhsLckNeur = NphsLckInfo.isPhsLckNeuron;


%% TALLY total number of unique neurons vs those with at least one phsLck

allNeurons = unique(Tselect.Unit_objectID);
nNeurons = numel(allNeurons);
nPhsLck = numel(neuPhsLck);
nNonLck = nNeurons - nPhsLck;

% pieData = [nNeuronsPhsLck / nNeuronsPshLck, (nNeurons - nNeuronsPhsLck)];
pieData(1) = round(100 * nPhsLck / nNeurons); % percentage
pieData(2) = round(100 * nNonLck / nNeurons); % percentage
pieLabels = {['(', num2str(pieData(1)), '%) Phase-Locked'], ...
             ['(', num2str(pieData(2)), '%) non-PhaseLocked']};


% display result in a pie chart:
f1 = figure; ax1 = axes;
pie(pieData, pieLabels);

title([ppar.subjID, ': Proportion of Phase-locking neurons (', ...
       num2str(nPhsLck), '/', num2str(nNeurons), ')']);


