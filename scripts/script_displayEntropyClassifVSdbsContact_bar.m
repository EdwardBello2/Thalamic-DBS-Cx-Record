% script for looking delta Entropy vs DBS frequency, showing continuous
% line for each individual neuron
% NOTE: still need to put in method for Contact-sweep too

% Cited:
% Moran, A., Stein, E., Tischler, H., Belelovsky, K. & Bar-Gad, 
% I. Dynamic Stereotypic Responses of Basal Ganglia Neurons to Subthalamic
% Nucleus High-Frequency Stimulation in the Parkinsonian Primate. 
% Frontiers in Systems Neuroscience 5, 21–21 (2011).



%% build pipeline-parameter struct "ppar"

clear; 

% PIPELINE PARAMETERS
% full path on local PC where project folder is (don't include subfolders here,
% that's tracked within the appropriate tables)
ppar.projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record'; % string

% full path on local PC where tables are to be loaded from (or saved to)
% ppar.tablePath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string
ppar.tablePath = 'C:\Users\bello043\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string

ppar.intDataRootFolder = 'DataProcessing\intermediateData';

% For getting sub-selection of data from tables
ppar.neuronType = 'SU';
ppar.subjID = 'Uva'; % 'Kramer' | 'Uva' | 'bothSubj'
ppar.trialType = 'ContactSweep'; % 'FrequencySweep' | 'ContactSweep'
% ppar.groups.dbsFrequency.groupFreqs = {10, 20, 30, 50, 100, 130};
% ppar.groups.dbsFrequency.groupLabels = {'hz10', 'hz20', 'hz30', 'hz50', 'hz100', 'hz130'};



% CONSTANTS

% Spike-sort quality:
REFRAC = 1 / 1000; % seconds, neuron refractory period in which spkes shouldn't occur
PERC_REFRACTHRESH = 1.0; % threshold for rejecting a unit with x% spikes within the refractory period


% Filter the data:
NEURON_MIN_TRIALS = 4;
SPARSE_HZ = 2;



% Define time ranges of various peri-DBS events:
SPKINTERV_PRE = [-30, 0]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS_1 = [0, 30]; %seconds, make sure 2nd element > 1st
SPKINTERV_DBS_2 = [30, 60]; %seconds, make sure 2nd element > 1st
SPKINTERV_POS = [60, 90]; %seconds, make sure 2nd element > 1st

% Entropy calculation:
BINS_PER_DECADE = 15;
ORD_H = 1;

ppar.BINS_PER_DECADE = BINS_PER_DECADE;
ppar.ORD_H = ORD_H;


% Bootrapped Entropy parameters:
NBOOTS = 10000;
% change significance for categorizing changes in decrH / incrH labels
SIG_PVAL = 0.0001 / 2; % to account for two-tailed distributions


% Final figure:
FIG_POSITION = [14 59 560 420];

t = now;
d = datetime(t, 'ConvertFrom', 'datenum');
YYYYMMDD = num2str(yyyymmdd(d));
YYMMDD = YYYYMMDD(3:end);



%%

% Get the name of the currently running script:
[scriptDirectoryFullPath, scriptName] = fileparts(mfilename('fullpath'));

% First make sure that this script has a folder within the project folder's
% intermediate data section
scriptIntermediateDataFolder = [ppar.projRootPath, '\', ppar.intDataRootFolder, '\', scriptName];
if ~exist(scriptIntermediateDataFolder, 'dir')
    mkdir(scriptIntermediateDataFolder)
    
end

% Perform the pipeline for every Row in tableRoot
% load root table "tableRoot"
load([ppar.tablePath, '\', 'tableRoot'], 'tableRoot');

%% Get subselection of tableRoot for analysis

% Tselect = filterTable_Master(tableRoot, scriptName, ppar);
tab = tableRoot;

% ppar.neuronType:
tab = filter_neuronType(tab, ppar.neuronType);

% ppar.subjID:
tab = filter_subjID(tab, ppar.subjID);

% ppar.trialType:
tab = filter_trialType(tab, ppar.trialType);


% Remove any units from processing that do not have good-enough sort
% quality:
Tselect = removePoorQualityUnits(tab, REFRAC, PERC_REFRACTHRESH, ppar);



%% Run analysis pipeline on each nexfile specified by table above

nNex = height(Tselect);

% Initialize columns to add to Tselect for final analysis in this script:
 FRpreCorr = zeros(nNex, 1);
      Hpre = zeros(nNex, 1);
FRdbsCorr1 = zeros(nNex, 1);
     Hdbs1 = zeros(nNex, 1);
FRdbsCorr2 = zeros(nNex, 1);
     Hdbs2 = zeros(nNex, 1);
 FRposCorr = zeros(nNex, 1);
      Hpos = zeros(nNex, 1);
Hdbs2_preBootAv = zeros(nNex, 1);
  Hdbs2_EmpPval = zeros(nNex, 1);


% For each nexfile contained within each record-row in Tselect, calculate
% the corrected Firing Rate and entropies for pre, dbs1, dbs2, and post
tic
for iNex = 1:nNex

    %% load nexfile of interest
%     iNex
    iNexObjectID = Tselect.objectID{iNex,1};
    nexfn = Tselect.nexFile{iNex,1};
    nexpn = Tselect.nexFileFolder{iNex,1};
    nexFile = readNexFile([ppar.projRootPath, '\', nexpn, '\', nexfn]);

    

    %% Calculate firing rates for portions of nexfile

    % Get firing rate calculations for various intervals
    FRpre = func_calcFRnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_PRE, ppar);
    FRdbs1 = func_calcFRnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_DBS_1, ppar);
    FRdbs2 = func_calcFRnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_DBS_2, ppar);
    FRpos = func_calcFRnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_POS, ppar);


    % save info in table columns
    FRpreCorr(iNex,1) = FRpre.rateCorrected; 
    FRdbsCorr1(iNex,1) = FRdbs1.rateCorrected;
    FRdbsCorr2(iNex,1) = FRdbs2.rateCorrected;
    FRposCorr(iNex,1) = FRpos.rateCorrected;

    

    %% Calculate entrpies for portions of nexfile

    % Get entropy calculations for various intervals
    H_PRE = func_calcHnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_PRE, ppar);
    H_DBS1 = func_calcHnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_DBS_1, ppar);
    H_DBS2 = func_calcHnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_DBS_2, ppar);
    H_POS = func_calcHnex(nexFile, Tselect.objectID{iNex}, SPKINTERV_POS, ppar);
    
    
    % Calculate bootstrapped distrib of pre-DBS entropies with ISI's
    % composed of spike events matching the DBS2 interval
    % Generate bootstrapped distribution of Pre-DBS isi's with isi
    % count matching that of the DBS case above 
    H_PREboots = func_calcHnex_bootstrap(nexFile, Tselect.objectID{iNex}, ...
        SPKINTERV_PRE, NBOOTS, SPKINTERV_DBS_2, ppar);
    H_PREboots(isnan(H_PREboots)) = []; % remove any NaN values 
    

    
    % save info in table columns

    Hdbs2_preBootAv(iNex,1) = mean(H_PREboots);
    Hdbs2_EmpPval(iNex,1) = getEmpiricalPval(H_PREboots, H_DBS2);
    
    
    Hpre(iNex,1) = H_PRE; 
    Hdbs1(iNex,1) = H_DBS1; 
    Hdbs2(iNex,1) = H_DBS2; 
    Hpos(iNex,1) = H_POS; 
    
    
    
end
toc



%% Prepare data for final visualization and analysis

% Append data of interest to Tselect
Tfinal = [Tselect, table(Hpre), ...
          table(Hdbs1), table(Hdbs2), table(Hdbs2_preBootAv), table(Hdbs2_EmpPval), ...
          table(Hpos), ...
          table(FRpreCorr), table(FRdbsCorr1), table(FRdbsCorr2), table(FRposCorr)];

numel(unique(Tfinal.Unit_objectID)) % display current number units       
Tfinal2 = Tfinal;

% Remove any rows that contain a NaN for entropy
Tfinal2(isnan(Tfinal2.Hpre),:) = [];
% Tfinal2(isnan(Tfinal2.Hdbs1),:) = [];
Tfinal2(isnan(Tfinal2.Hdbs2),:) = [];
% Tfinal2(isnan(Tfinal2.Hpos),:) = [];
numel(unique(Tfinal2.Unit_objectID))



% Remove any rows that have a firing rate of less than min firing rate:
% filter out any rows that have too-sparse cells, based on 30-second
% portions
isPreRateBelow = Tfinal2{:, 'FRpreCorr'} < SPARSE_HZ;
isDbs1RateBelow = Tfinal2{:, 'FRdbsCorr1'} < SPARSE_HZ;
isDbs2RateBelow = Tfinal2{:, 'FRdbsCorr2'} < SPARSE_HZ;
isRemove = isPreRateBelow | isDbs2RateBelow;
Tfinal2(isRemove, :) = [];
numel(unique(Tfinal2.Unit_objectID))


% Remove any neurons that have less than minimum trials 
Tfinal3 = filter_neuronMinimumTrials(Tfinal2, NEURON_MIN_TRIALS);
numel(unique(Tfinal3.Unit_objectID))

Tanalyze = Tfinal3;



%% Hdbs2: For each trial, categorize as Entorpy incr, decr, and significance
% based on bootstrapped distributions of pre-DBS entropy

dbsPortionStr = 'Hdbs2';

% % change significance
% SIG_PVAL = 0.05 / 2; % to account for two-tailed distributions

% display numbers of cells for each frequency in command line:
contactsADD = Tanalyze{:, 'dbsContact'};
Entropy_grpLabelADD = cellstr(contactsADD);
disp_NpointsVsFreq(Entropy_grpLabelADD)


Tdisp = Tanalyze; % just to make typing easier...






% CATEGORIZE ALL ENTROPY CHANGES
nRows = height(Tdisp);

     wasDecrH = false(nRows, 1);
     wasIncrH = false(nRows, 1);
signifChangeH = false(nRows, 1);

for iRow = 1:nRows
    if Tdisp.Hdbs2(iRow) <= Tdisp.Hdbs2_preBootAv(iRow) % entorpy decrease
        wasDecrH(iRow) = true; 
        
    else % entropy increase
        wasIncrH(iRow) = true;
        
    end
    
    if abs(Tdisp.Hdbs2_EmpPval(iRow)) < SIG_PVAL
        signifChangeH(iRow) = true;
        
    end
    
end

Tdisp = [Tdisp, table(wasDecrH), table(wasIncrH), table(signifChangeH)];


% Hdbs2: Display the group proportions according to frequency
% pre-allocate matrix for counts with each frequency being a row; columns:
% wasDecrH, wasIncrH, noCh
contacts = {'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};
nConts = numel(contacts);
catCountsH = zeros(nConts, 3);

% Gather all counts
for iCnt = 1:nConts
    % get subselection of table rows pertaining to frequency "iFr"
    T_byContact = Tdisp(strcmp(Tdisp.dbsContact, contacts{iCnt}),:);
    
    % categorize each trial into "incrH", "decrH", "noCh" (in that order), 
    % keep count
    catCountsH(iCnt,:) = getCategoryCounts_table(T_byContact);
  
end




% DISPLAY CATEGORIZATION RESULTS AS BARPLOTS
 
% transfer to percentages
tot = sum(catCountsH, 2);
catPercs = catCountsH ./ tot;
catPercs(:,3) = []; % remove the noChange category from display
% display results as barplots
figure; ax = axes;
bar(catPercs * 100)
ax.XTickLabel = contacts;
grid on
xlabel('DBS Contact')
ylabel('% recorded cells')
lgd = legend('H increase', 'H decrease', ...
      'Location', 'northeastoutside') ;
tit = title([ppar.subjID, ':  population percentage by ', dbsPortionStr,' Entropy change']);
ax.YLim = [0, 60];

% alternate display as lines, noChange category excluded
figure; ax = axes;
plot(100 * catPercs(:,1), 'o--'); hold on;
plot(100 * catPercs(:,2), 'o--');
ax.YLim = [0, 60];
ax.XLim = [0.5, 8.5];
ax.XTickLabel = contacts;
grid on
xlabel('DBS contact')
ylabel('% recorded cells')
lgd = legend('H increase', 'H decrease', ...
      'Location', 'northeastoutside') ;
tit = title([ppar.subjID, ':  population percentage by ', dbsPortionStr,' Entropy change']);


% EXTRA VISUALIZATION: SEE WHAT THE COMBINED NUMBERS OF HDECR AND HINCR ARE
% LIKE AS PERCENTAGES


% DISPLAY CATEGORIZATION RESULTS AS BARPLOTS
 
% transfer to percentages
catCountsHco = [sum(catCountsH(:,1:2), 2), catCountsH(:,3)];
tot = sum(catCountsHco, 2);
catPercs = catCountsHco ./ tot;
catPercs(:,end) = []; % remove the noChange category from display
% display results as barplots
figure; ax = axes;
bar(catPercs * 100)
ax.XTickLabel = contacts;
grid on
xlabel('DBS Contact')
ylabel('% recorded cells')
lgd = legend('H increase', 'H decrease', ...
      'Location', 'northeastoutside') ;
tit = title([ppar.subjID, ':  population ALL CHANGE percentage by ', dbsPortionStr,' Entropy change']);
ax.YLim = [0, 60];

% alternate display as lines, noChange category excluded
figure; ax = axes;
plot(100 * catPercs(:,1), 'o--'); hold on;
% plot(100 * catPercs(:,2), 'o--');
ax.YLim = [0, 60];
ax.XLim = [0.5, 8.5];
ax.XTickLabel = contacts;
grid on
xlabel('DBS contact')
ylabel('% recorded cells')
lgd = legend('H increase', 'H decrease', ...
      'Location', 'northeastoutside') ;
tit = title([ppar.subjID, ':  population ALL CHANGE percentage by ', dbsPortionStr,' Entropy change']);




% STATS TESTS FOR PERCENTAGE CATEGORIZATION
% Chi-squared test for independence for each group:
% 
% % noChange:
% sz = [6, 2];
% varTypes = {'double', 'double'};
% tab = table('Size', sz, 'VariableTypes', varTypes);
% tab.Properties.RowNames = {'10', '20', '30', '50', '100', '130'};
% tab.Properties.VariableNames = {'Change', 'noChange'};
% Ch = sum(catCounts(:,1:2), 2);
% noCh = catCounts(:,3);
% tab{:,1} = Ch;
% tab{:,2} = noCh;
% 
% [h, p, stats] = fishertest(tab)


% try crosstab for chi2 test for independence of all groupings:
% first create label column in table:
Hchange = cell(nRows, 1);

for iRow = 1:nRows
    if Tdisp.signifChangeH(iRow)
        if Tdisp.wasIncrH(iRow)
            Hchange{iRow} = 'Hincr';
            
        elseif Tdisp.wasDecrH(iRow)
            Hchange{iRow} = 'Hdecr';
            
        else
            Hchange{iRow} = 'other';
            
        end
        
    else
        Hchange{iRow} = 'noChange';
        
    end
       
end
Ttest = [Tdisp, table(Hchange)];

[xtab, chi2, p, labels] = crosstab(Ttest.dbsContact, Ttest.Hchange)


% try crosstab for chi2 test for independence further refined to just Hdecr:
% first create label column in table:
Hchange = cell(nRows, 1);

for iRow = 1:nRows
    if Tdisp.signifChangeH(iRow)
        if Tdisp.wasDecrH(iRow)
            Hchange{iRow} = 'Hdecr';
            
        else
            Hchange{iRow} = 'other';
            
        end
        
    else
        Hchange{iRow} = 'other';
        
    end
       
end
Ttest = [Tdisp, table(Hchange)];

[xtab, chi2, p, labels] = crosstab(Ttest.dbsContact, Ttest.Hchange)


% try crosstab for chi2 test for independence further refined to just Hincr:
% first create label column in table:
Hchange = cell(nRows, 1);

for iRow = 1:nRows
    if Tdisp.signifChangeH(iRow)
        if Tdisp.wasIncrH(iRow)
            Hchange{iRow} = 'Hincr';
            
        else
            Hchange{iRow} = 'other';
            
        end
        
    else
        Hchange{iRow} = 'other';
        
    end
       
end
Ttest = [Tdisp, table(Hchange)];

[xtab, chi2, p, labels] = crosstab(Ttest.dbsContact, Ttest.Hchange)



%% SUB-FUNCTIONS

function [rowCounts] = getCategoryCounts_table(T)

% Gather all counts
TnoCh = T(~T.signifChangeH,:);
numNoCh = height(TnoCh);

Tchange = T(T.signifChangeH,:);

  TdecrH = Tchange(Tchange.wasDecrH,:);
numDecrH = height(TdecrH);

  TincrH = Tchange(Tchange.wasIncrH,:);
numIncrH = height(TincrH);

rowCounts = [numIncrH, numDecrH, numNoCh];

end

function [EmpPval] = getEmpiricalPval(normDistrib, testValue)
% outputs a - or + number that tells how far "value" is from the mean of
% "normDistrib" much like how p-value works. Assumes that the randomly
% distributed data in "normDistrib" is normally distributed, but this is
% not strictly necessary. Interpret with caution. 

isLess = normDistrib < testValue;

Pval = sum(isLess) / numel(normDistrib);

if Pval > 0.5
    EmpPval = 1 - Pval;
    
else
    EmpPval = -Pval;
    
end


end

function tabFilt = removePoorQualityUnits(Tselect, REFRAC, PERC_REFRACTHRESH, ppar)
% Check if a given unit has good-enough sort quality to include in
% analysis, remove those that are low quality:
[unitIDs, unitISIs, percInRefrac] = func_assessUnitSpkInRefrac(Tselect, REFRAC, ppar);
isQualityUnit = percInRefrac < PERC_REFRACTHRESH;
tabUnitQ = [table(unitIDs), table(isQualityUnit)];
T = join(Tselect, tabUnitQ, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'unitIDs'); 
T(~T.isQualityUnit,:) = []; % REMOVE THOSE ROWS PERTAINING TO POOR QUALITY UNIT
tabFilt = T;

disp(['Removed ', num2str(sum(~isQualityUnit)), ' of ', ...
       num2str(numel(unitIDs)),' units due to low sort quality']);
  

end

function tabFilt = filter_neuronType(tab, neuronType)

isSU = strcmp(tab.NeuronType, neuronType);
tabFilt = tab(isSU,:);

end

function tabFilt = filter_subjID(tab, subjID)
if strcmp(subjID, 'bothSubj')
    % do nothing (keep data across both subjects)
    tabFilt = tab;
    
elseif strcmp(subjID, 'Kramer') | strcmp(subjID, 'Uva')
    isSubjID = strcmp(tab.subjID, subjID);
    tabFilt = tab(isSubjID,:);
    
else
    error('Oh snap! wrong string for ppar.subjID')

end

end

function tabFilt = filter_trialType(tab, trialType)
% assumes that table has two variables: 1) ContactSweep, 2) FrequencySweep,
% both of which are columns of 1's or 0's. Currently that is how I indicate
% whether a given RECORD represents a trial for Contact or Frequency sweep.

if strcmp(trialType, 'FrequencySweep')
    isTrialSwp = logical(tab.FrequencySweep);
    tabFilt = tab(isTrialSwp,:);

    % Remove any trials where frequency sweep was not done on electrode C0
    isC0 = strcmp(tabFilt.dbsContact, 'C0');
    tabFilt(~isC0,:) = [];
    
elseif strcmp(trialType, 'ContactSweep')
    isTrialSwp = logical(tab.ContactSweep);
    tabFilt = tab(isTrialSwp,:);

    % Remove any trials where Contact sweep was not done with 130Hz freq
    is130Hz = tabFilt.dbsFrequency == 130;
    tabFilt(~is130Hz,:) = [];

else
    error('Wrong values for ppar.TrialType');

end


end

function [intervEvents] = getIntervalEvents(evTimes, refT, refInterval)
% returns the observed firing rate within the desired time interval.
% "refInterval" indicates the window within which to gather spike times 
% works with the values in "spkTimes" and "StimTs" as inputs

rerefTimes = evTimes - refT;

% get observed spike rate "FRobs" based on count and time duration
idxInInterv = (rerefTimes >= refInterval(1)) & (rerefTimes < refInterval(2));
rerefSubselect = rerefTimes(idxInInterv);
intervEvents = rerefSubselect + refT;



end

function [n] = calcFiringRates(spkTimes, StimTs, spkTimeInterval)
% take the "spkTimes" and "StimTs" data from the nexFile and calculate the
% observed spike firing rate and other information occurring within the 
% time interval "spkTimeInterval" (1x2 array, seconds, refenced to DBS 
% onset time).

% CONSTANTS
% assumed to be true for both Uva and Kramer
periStimBlank = 1 / 1000; %s, blanking time around each stim pulse


% CODE
refT = StimTs.DBS(1);

% Calculate observed spike rate for interval
[spkTimesDbs] = getIntervalEvents(spkTimes, refT, spkTimeInterval);
totIntervTime = spkTimeInterval(2) - spkTimeInterval(1);
frDbsObs = numel(spkTimesDbs) / totIntervTime;


% Estimate the "true" firing rate by correcting for DBS artifact blank
% time, as done in Moran et al 2011
stmTimesDbs = getIntervalEvents(StimTs.DBS, refT, spkTimeInterval);
totBlankTime = periStimBlank * numel(stmTimesDbs);
frDbsCorr = frDbsObs * (totIntervTime / (totIntervTime - totBlankTime));


n.refT = refT;
n.spkTimesDbs = spkTimesDbs;
n.stmTimesDbs = stmTimesDbs;
n.rateObserved = frDbsObs; % spikes/second
n.rateCorrected = frDbsCorr; % spikes/second
n.totIntervTime = totIntervTime;
n.totBlankTime = totBlankTime;

end


function tabFilt = filter_neuronMinimumTrials(tab, neuronMinimumTrials)
% Remove any neurons from analysis that do not have at least 4 trials
% present in the current selection table

% count the number of times each neuron shows up and store in new table
% called "NeuronCounts"
neurons = tab.Unit_objectID;
[uNeurons, ~, uNeuIdx] = unique(neurons);

nNeurons = length(uNeurons); % number of unique neurons
neuCount = zeros(nNeurons, 1);
for iNeu = 1:nNeurons
    neuCount(iNeu) = sum(uNeuIdx == iNeu);

end

NeuronCounts = [table(uNeurons), table(neuCount)];

% join this table to current selection table
% leftKey = find(strcmp(Nselect.Properties.VariableNames, 'Unit_objectID'));
tab = join(tab, NeuronCounts, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'uNeurons');

% remove those rows with under "n" trials
isPresentforTrials = (tab.neuCount >= neuronMinimumTrials);
tabFilt = tab;
tabFilt(~isPresentforTrials,:) = [];

end

function FRidx = calcFRindex(FRbase, FRtest)
% firing rate index. 0 means no change; 1 means cell went from zero
% baseline to "something"; -1 means cell went from some activity to
% complete inhibition.

FRidx = (FRtest - FRbase) ./ (FRbase + FRtest);

end

function disp_NpointsVsFreq(labels)
% show the numbers of neuron observations for each grouping-label
uniqueLabels = sort(unique(labels));
numLabels = numel(uniqueLabels);

% prefill 

for iLab = 1:numel(uniqueLabels)
    iLabStr = uniqueLabels{iLab};
    num_iLab = sum(strcmp(iLabStr, labels));
    
    disp([iLabStr, ': ', num2str(num_iLab)])
    
    
end


end










