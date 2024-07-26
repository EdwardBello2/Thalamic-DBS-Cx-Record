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
CURR_FUNC = 'script_displayPhaseLockPercVsdbsTrialCond_barplot'; 



%% LOAD all relevant tables and MERGE them

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



%% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);



%% Detect those neuron trial-observations that phase-locked vs those that did not, label rows
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



%% PLOT stacked barplot population % that are phase locked


 % Gather all unique DBS Frequencies
freqs = unique(Tselect.dbsFrequency);
nFreqs = numel(freqs);


% Gather up all rateChange types coutns for each frequency
Fcount = zeros(nFreqs, 2);
for i = 1:nFreqs
    isHz = Tselect.dbsFrequency == freqs(i);
    HzR = Tselect(isHz,:);
    nPhsLck = sum(HzR.PhaseLocked);
    notPhsLck = sum(~HzR.PhaseLocked);

    Fcount(i,1:2) = [nPhsLck, notPhsLck];

end

Fpercent = 100 * (Fcount ./ sum(Fcount, 2));


figure; ax = axes;
b = bar(Fpercent, 'stacked');
b(1).FaceColor = [0.5, 0.5, 0.5]; 
b(2).FaceColor = [1.0, 1.0, 1.0];
% b(3).FaceColor = [0.0, 0.0, 0.0];
set(gca, 'xticklabel', freqs)
ax.YLim = [0, 100];

legend('PhaseLocked', 'notPhaseLocked', 'Location', 'northeastoutside')

title(['% of PhaseLocking Neurons in ', ppar.subjID, ' for ', ppar.trialType]);
ylabel('% neurons recorded')
xlabel('DBS Frequency (Hz)')



%% Tally up neuron recording counts in a final table according to:
% 1) Phase-locked
% 2) NON-locked
% 3) TotNeurons Obs

labelStr = cellstr(num2str(freqs));

for i = 1:numel(labelStr)
    noSpaceStr = labelStr{i}(~isspace(labelStr{i})) ;
    labelStr{i} = ['hz', noSpaceStr];
    
end

tab = Fcount';
tot = sum(tab, 1);
tab = [tab; tot];


sz = [3, 6];
varTypes = {'double', 'double', 'double', 'double', 'double', 'double'};
NeuronCounts = table('Size', sz, 'VariableTypes', varTypes);
NeuronCounts.Variables = tab;
NeuronCounts.Properties.VariableNames = labelStr;
NeuronCounts.Properties.RowNames = {'PhsLock', 'NonLock', 'TotNeurons'}









