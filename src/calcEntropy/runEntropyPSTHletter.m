function runEntropyPSTHletter(pipeParams)
% Run PSTH-based Entropy Letter-Method pipeline
%
% This is based off of Joe Xiao's work on the "letter" method of entropy
% estimation, see [Xiao et al paper].
% 
% This pipeline assumes that the following tables exist and are correct:
%
% NEXprocfiles_subjID.mat
% SortedUnits_sibjID.mat
% SweepAnalysisTrials4Paper2018_subjID.mat
%
% where "subjID" is the nhp name (i.e. 'Uva')

%% DEFAULT PARAMETERS

DEFAULT.subjID        = 'Uva';
DEFAULT.tablepn       = 'K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
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
H_PSTHtableName = [baseLabel, name];


% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createEntropyPsthLetterTableIfNeeded(NEX_2analyze, H_PSTHtableName, pipeParams);


% LOAD the specified table: H_letterResults
load([pipeParams.tablepn, '\', H_PSTHtableName, '.mat' ]);



%% GET boolean indices of all trials according to labels and modulations

R = H_letterResults; % simplify variable name
alpha = pipeParams.pValAlpha;


% Divide Analyze table into Contact-sweep and Current-sweep
isCswp = strcmp(R.sweepType(:), 'Contact');
R_Cswp = R(isCswp, :);

isFswp = strcmp(R.sweepType(:), 'Frequency');
R_Fswp = R(isFswp, :);


% Get labels and idx's where they occur for electrode Contact-sweeps
[Clabels, isRowClabel] = getAlldbsElectrodeLabels(R_Cswp);
isModCswp = R_Cswp.pVal(:) < alpha;
isRowClabelMod = isRowClabel & isModCswp;
isRowClabelnonMod = isRowClabel & ~isModCswp;


% Get labels and idx's where they occur for Frequency-sweeps
[Flabels, isRowFlabel] = getAlldbsFrequencyLabels(R_Fswp);
isModFswp = R_Fswp.pVal(:) < alpha;
isRowFlabelMod = isRowFlabel & isModFswp;
isRowFlabelnonMod = isRowFlabel & ~isModFswp;



%% Calculate Delta-Entropies for PSTH-based Entropy of all trials

% For Contact-sweep:
nRows = size(R_Cswp, 1);
H_PreAv = zeros(nRows, 1);
for iRow = 1:nRows
    H_PreAv(iRow,1) = mean(R_Cswp.H_PREbootdistr{iRow,1});
    
end
deltaH_Cswp = R_Cswp.H_DBS(:) - H_PreAv;


% For Frequency-sweep:
nRows = size(R_Fswp, 1);
H_PreAv = zeros(nRows, 1);
for iRow = 1:nRows
    H_PreAv(iRow,1) = mean(R_Fswp.H_PREbootdistr{iRow,1});
    
end
deltaH_Fswp = R_Fswp.H_DBS(:) - H_PreAv;



%% PLOT overall % of trials that show phase-locking vs none


% Display both sweeps % modulated results:
f1 = figure;
dispModPercent_Csweep(Clabels, isRowClabel, isRowClabelMod, alpha)
title(['Phase-locked Contact-sweep trial-population (p<', num2str(alpha), ')']);
xlabel('% of trials')
ylabel('DBS electrode')
legend('phase-locked', 'non P-L');


f2 = figure;
dispModPercent_Fsweep(Flabels, isRowFlabel, isRowFlabelMod, alpha)
title(['Phase-locked Frequency-sweep trial-population (p<', num2str(alpha), ')']);
xlabel('DBS Frequency')
ylabel('% of trials')
legend('phase-locked', 'non P-L');



%% PLOT boxplots of delta Entropy for modulated and non-modulated

% For Contact-sweep:
rowClabel = R_Cswp.dbsElectrode(:);

rowClabelMod = rowClabel(isModCswp);
deltaH_CswpMod = deltaH_Cswp(isModCswp);

rowClabelNonMod = rowClabel(~isModCswp);
deltaH_CswpNonMod = deltaH_Cswp(~isModCswp);

f3 = figure;
dispDeltaH_Cswp_Boxplot(deltaH_CswpMod, rowClabelMod, Clabels);
title('Phase-lock Entropy change in Modulated Contact-Sweep trials');

f4 = figure; 
dispDeltaH_Cswp_Boxplot(deltaH_CswpNonMod, rowClabelNonMod, Clabels);
title('Phase-lock Entropy change in NON-Modulated Contact-Sweep trials');


% For Frequency-sweep:
rowFrequency = R_Fswp.dbsFrequency(:);
nRows = size(rowFrequency, 1);
rowFlabel = cell(nRows, 1);
for i = 1:nRows, rowFlabel{i,1} = num2str(rowFrequency(i,1)); end


rowFlabelMod = rowFlabel(isModFswp);
deltaH_FswpMod = deltaH_Fswp(isModFswp);

rowFlabelNonMod = rowFlabel(~isModFswp);
deltaH_FswpNonMod = deltaH_Fswp(~isModFswp);

f5 = figure;
dispDeltaH_Fswp_Boxplot(deltaH_FswpMod, rowFlabelMod, Flabels);
title('Phase-lock Entropy change in Modulated Frequency-Sweep trials');

f6 = figure; 
dispDeltaH_Fswp_Boxplot(deltaH_FswpNonMod, rowFlabelNonMod, Flabels);
title('Phase-lock Entropy change in NON-Modulated Frequency-Sweep trials');



%% Final optional saving of figures

if pipeParams.finalFigs.save == true
    savfPn = pipeParams.finalFigs.savepn;
    
    % save as .eps files
    saveas(f1, [savfPn, '\percModpsthTrials_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f2, [savfPn, '\percModpsthTrials_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');
    saveas(f3, [savfPn, '\dHletterAll_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f4, [savfPn, '\dHletterAll_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');

    % save as .jpg files
    saveas(f1, [savfPn, '\percModpsthTrials_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f2, [savfPn, '\percModpsthTrials_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');
    saveas(f3, [savfPn, '\dHletterAll_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f4, [savfPn, '\dHletterAll_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');

end

end % END function



%% SUB-FUNCTIONS

function [labels, isLab, isLabMod] = findCsweepModLabels(R_Cswp, alpha)



% Get idx's of modulated (phase-locked) trials
isModCswp = R_Cswp.pVal(:) <= alpha;


% Get positions of each label in Table
labels = unique(R_Cswp.dbsElectrode);
labelTable = R_Cswp.dbsElectrode(:);


% Get total number of trials for each label
nTableRows = size(R_Cswp, 1);
nLabs = numel(labels);

isLab = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLab(:,iLab) = strcmp(labels(iLab), labelTable);

end



% Get total number of modulated trials for each label
isLabMod = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLabMod(:,iLab) = (isLab(:,iLab) & isModCswp);
    
end



end % END function

function [labels, isLab, isLabMod] = findFsweepModLabels(R_Fswp, alpha)



% Get idx's of modulated (phase-locked) trials
isModFswp = R_Fswp.pVal(:) <= alpha;


% Get positions of each label in Table
labels = unique(R_Fswp.dbsFrequency);
labelTable = R_Fswp.dbsFrequency(:);


% Get total number of trials for each label
nTableRows = size(R_Fswp, 1);
nLabs = numel(labels);

isLab = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLab(:,iLab) = (labelTable == labels(iLab));

end



% Get total number of modulated trials for each label
isLabMod = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLabMod(:,iLab) = (isLab(:,iLab) & isModFswp);
    
end










end

% function dispModPercent_Csweep(Clabels, isClab, isClabMod, alpha)
% 
% Ctot = sum(isClab)';
% Cmod = sum(isClabMod)';
% 
% Cpercent(:,1) = 100 * (Cmod ./ Ctot);
% Cpercent(:,2) = 100 * (1 - (Cmod ./ Ctot));
% 
% 
% % figure; 
% b = barh(Cpercent, 'stacked');
% b(1).FaceColor = [0.5, 0.5, 0.5];
% b(2).FaceColor = [1, 1, 1];
% set(gca,'yticklabel',Clabels)
% 
% title(['Phase-locked Contact-sweep trial-population (p<', num2str(alpha), ')']);
% xlabel('% of trials')
% ylabel('DBS electrode')
% legend('phase-locked', 'non P-L');
% 
% end

% function dispModPercent_Fsweep(Flabels, isFlab, isFlabMod, alpha)
% 
% 
% Ftot = sum(isFlab)';
% Fmod = sum(isFlabMod)';
% 
% Fpercent(:,1) = 100 * (Fmod ./ Ftot);
% Fpercent(:,2) = 100 * (1 - (Fmod ./ Ftot));
% 
% 
% % figure; 
% b = bar(Fpercent, 'stacked');
% b(1).FaceColor = [0.5, 0.5, 0.5];
% b(2).FaceColor = [1, 1, 1];
% set(gca,'xticklabel',Flabels)
% 
% title(['Phase-locked Frequency-sweep trial-population (p<', num2str(alpha), ')']);
% xlabel('DBS Frequency')
% ylabel('% of trials')
% legend('phase-locked', 'non P-L');
% 
% 
% 
% end

function boxplotDeltaHMod_Cswp(R_Cswp, Clabels, isClabMod)
% Get Csweep deltaH values for a given label
nLab = numel(Clabels);
boxDeltaH = [];
boxGroup = {};
for iLab = 1:nLab

    lab = Clabels{iLab};

    % First get Average Hpre from the bootstrapped Hpre estimates
    Hpreboot = R_Cswp.H_PREbootdistr(isClabMod(:,iLab));

    nHpre = size(Hpreboot, 1);
    HpreAv = zeros(nHpre, 1);

    for iHpre = 1:nHpre
        HpreAv(iHpre, 1) = mean(Hpreboot{iHpre,1});

    end

    % Next get empirically measured Hdbs
    HdbsEmp = R_Cswp.H_DBS(isClabMod(:,iLab));

    % Finally get deltaH numbers for given label
    deltaH = HpreAv - HdbsEmp;
    labelGroup = cell(nHpre, 1); for i = 1:nHpre, labelGroup{i,1} = lab; end

    
    %update the running tally of deltaH's and respective labels
    boxDeltaH = [boxDeltaH; deltaH];
    boxGroup = [boxGroup; labelGroup];
    
end
% Display Csweep modulated delta H


% figure;
boxplot(boxDeltaH, boxGroup);
title('Entropy change in Modulated Contact-Sweep trials');
ylabel('delta-H (bits/spike)');
xlabel('Electrode Contact (all at 130Hz)');
hLine = refline(0,0); 
hLine.LineStyle = '--';
hLine.Color = [0,0,0];



end

function boxplotDeltaHMod_Fswp(R_Fswp, Flabels, isFlabMod)
% Get Csweep deltaH values for a given label
nLab = numel(Flabels);
boxDeltaH = [];
boxGroup = {};
for iLab = 1:nLab

    lab = Flabels(iLab);

    % First get Average Hpre from the bootstrapped Hpre estimates
    Hpreboot = R_Fswp.H_PREbootdistr(isFlabMod(:,iLab));

    nHpre = size(Hpreboot, 1);
    HpreAv = zeros(nHpre, 1);

    for iHpre = 1:nHpre
        HpreAv(iHpre, 1) = mean(Hpreboot{iHpre,1});

    end

    % Next get empirically measured Hdbs
    HdbsEmp = R_Fswp.H_DBS(isFlabMod(:,iLab));

    % Finally get deltaH numbers for given label
    deltaH = HpreAv - HdbsEmp;
    labelGroup = cell(nHpre, 1); for i = 1:nHpre, labelGroup{i,1} = lab; end

    
    %update the running tally of deltaH's and respective labels
    boxDeltaH = [boxDeltaH; deltaH];
    boxGroup = [boxGroup; labelGroup];
    
end
% Display Csweep modulated delta H


% figure;
boxplot(boxDeltaH, boxGroup);
title('Entropy change in Modulated Frequency-Sweep trials');
ylabel('delta-H (bits/spike)');
xlabel('DBS Frequency (all at C0)');
hLine = refline(0,0); 
hLine.LineStyle = '--';
hLine.Color = [0,0,0];



end
