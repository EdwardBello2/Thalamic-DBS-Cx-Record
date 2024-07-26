function runEntropyPSTH_percentPopPhaseLocked(pipeParams)
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
% TrialInfo = SweepAnalysisTrials4Paper2018;



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


% Get labels and idx's where they occur for Frequency-sweeps
[Flabels, isRowFlabel] = getAlldbsFrequencyLabels(R_Fswp);
isModFswp = R_Fswp.pVal(:) < alpha;
isRowFlabelMod = isRowFlabel & isModFswp;



%% Calculate Delta-Entropies for PSTH-based Entropy of all trials

% For Contact-sweep:
nRows = size(R_Cswp, 1);
H_PreAv = zeros(nRows, 1);
for iRow = 1:nRows
    H_PreAv(iRow,1) = mean(R_Cswp.H_PREbootdistr{iRow,1});
    
end


% For Frequency-sweep:
nRows = size(R_Fswp, 1);
H_PreAv = zeros(nRows, 1);
for iRow = 1:nRows
    H_PreAv(iRow,1) = mean(R_Fswp.H_PREbootdistr{iRow,1});
    
end


%% PERFORM Fisher's exact test on the data

% Contact-sweep:

% Build a table of Csweep counts

C_All = sum(isRowClabel)';
C_PhsLck = sum(isRowClabelMod)';
C_NonPhsLck = C_All - C_PhsLck;

% P = table(PhsLck);
% NP = table(NonPhsLck);
% 
% Cross = [P, NP];
% Cross.Properties.RowNames = Clabels
% 
% Cpercent(:,2) = 100 * (1 - (Cmod ./ Ctot));



% 
% 
% % figure; 
% b = barh(Cpercent);
% b(1).FaceColor = [0.5, 0.5, 0.5];
% % b(2).FaceColor = [1, 1, 1];
% set(gca, 'yticklabel', Clabels)
% set(gca, 'XLim', pipeParams.PhsLckPercentXLim);




% Frequeny-sweep:
F_All = sum(isRowFlabel)';
F_PhsLck = sum(isRowFlabelMod)';
F_NonPhsLck = F_All - F_PhsLck;




%% PLOT overall % of trials that show phase-locking vs none


% Display both sweeps % modulated results:
f1 = figure;
% dispModPercent_Csweep(Clabels, isRowClabel, isRowClabelMod, alpha)
dispPhsLckPercent_Csweep(Clabels, isRowClabel, isRowClabelMod, pipeParams);
title(['Phase-locked Contact-sweep trial-population (p<', num2str(alpha), ')']);
xlabel('% of trials')
ylabel('DBS electrode')
% legend('phase-locked', 'non P-L');


f2 = figure;
% dispModPercent_Fsweep(Flabels, isRowFlabel, isRowFlabelMod, alpha)
dispPhsLckPercent_Fsweep(Flabels, isRowFlabel, isRowFlabelMod, pipeParams)

title(['Phase-locked Frequency-sweep trial-population (p<', num2str(alpha), ')']);
xlabel('DBS Frequency')
ylabel('% of trials')
% legend('phase-locked', 'non P-L');



%% Final optional saving of figures

if pipeParams.finalFigs.save == true
    savfPn = pipeParams.finalFigs.savepn;
    
    % save as .eps files
    saveas(f1, [savfPn, '\percModpsthTrials_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f2, [savfPn, '\percModpsthTrials_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');


    % save as .jpg files
    saveas(f1, [savfPn, '\percModpsthTrials_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f2, [savfPn, '\percModpsthTrials_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');


end

end % END function


