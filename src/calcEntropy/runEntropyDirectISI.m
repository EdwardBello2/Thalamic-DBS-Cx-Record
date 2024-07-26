function runEntropyDirectISI(pipeParams)
% Run ISI-based Entropy Estimation Pipeline
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
DEFAULT.ordH          = 2;
DEFAULT.binsPD        = 20;



%% SET DEFAULT PARAMETERS IF USER HAS NOT SET THEM

% Check for specific User-defined pipeline parameter inputs:
if ~isfield(pipeParams, 'subjID'), pipeParams.subjID               = DEFAULT.subjID; end
if ~isfield(pipeParams, 'tablepn'), pipeParams.tablepn             = DEFAULT.tablepn; end
if ~isfield(pipeParams, 'neuTypeFilter'), pipeParams.neuTypeFilter = DEFAULT.neuTypeFilter; end
if ~isfield(pipeParams, 'predbsTime'), pipeParams.predbsTime       = DEFAULT.predbsTime; end
if ~isfield(pipeParams, 'dbsTime'), pipeParams.dbsTime             = DEFAULT.dbsTime; end
if ~isfield(pipeParams, 'hzThresh'), pipeParams.hzThresh           = DEFAULT.hzThresh; end
if ~isfield(pipeParams, 'ordH'), pipeParams.ordH                   = DEFAULT.ordH ; end
if ~isfield(pipeParams, 'binsPD'), pipeParams.binsPD               = DEFAULT.binsPD; end



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



%% CREATE/LOAD a table with Direct-Entropy estimates for both pre- and DBS periods

% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyDirectISI';
name = buildEntropyDirectISITableName(pipeParams);
AnalysisTableName = [baseLabel, name];


% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table
createEntropyDirectISITableIfNeeded(NEX_2analyze, AnalysisTableName, pipeParams);


% LOAD the specified Table: H_DirectISIResults
load([pipeParams.tablepn, '\', AnalysisTableName, '.mat']);



%% PLOT boxplots of entropy results for both sweeps-types

% comment out the below for now:


R = H_DirectISIResults; % simplify variable name
alpha = pipeParams.pValAlpha;


% Divide Analyze table into Contact-sweep and Current-sweep
isCswp = strcmp(R.sweepType(:), 'Contact');
R_Cswp = R(isCswp, :);

isFswp = strcmp(R.sweepType(:), 'Frequency');
R_Fswp = R(isFswp, :);


% Get labels and idx's where they occur for C-sweeps
[Clabels, isClab, isClabMod] = findCsweepModLabels(R_Cswp, alpha);

% Get labels and idx's where they occur for F-sweeps
[Flabels, isFlab, isFlabMod] = findFsweepModLabels(R_Fswp, alpha);




% Display both sweeps % modulated results:
f1 = figure;
dispModPercent_Csweep(Clabels, isClab, isClabMod, alpha)
title(['ISI pattern-change Contact-sweep trial-population (p<', num2str(alpha), ')']);
xlabel('% of trials')
ylabel('DBS electrode')
legend('pattern-mod', 'non P-M');


f2 = figure;
dispModPercent_Fsweep(Flabels, isFlab, isFlabMod, alpha)
title(['ISI pattern-change Frequency-sweep trial-population (p<', num2str(alpha), ')']);
xlabel('DBS Frequency')
ylabel('% of trials')
legend('pattern-mod', 'non P-M');

% % 
% % % Plot the deltaHs of both sweep-types
% % deltaH_boxplots(H_Results)
% 
% f1 = figure;
% boxplot_deltaH_Fswp(H_Results)
% 
% 
% f2 = figure;
% boxplot_deltaH_Fswp(H_Results)
% 
% 
% 
% 
% 
% f3 = figure;
% ax = axes;
% barplot_deltaH_Cswp(H_Results);
% ax.YGrid = 'on';
% ax.GridLineStyle = '-';
% yLim = get(gca, 'YLim');
% yLim(1) = pipeParams.maxdH;
% set(gca, 'YLim', yLim);
% 
% f4 = figure; 
% ax = axes;
% barplot_deltaH_Fswp(H_Results);
% ax.YGrid = 'on';
% ax.GridLineStyle = '-';
% yLim = get(gca, 'YLim');
% yLim(1) = pipeParams.maxdH;
% set(gca, 'YLim', yLim);
% 
% 
% % Get paired-sample T-test results for each condition in Contact-sweep
% 
% pValCsweep = deltaH_ttest(H_Results, 'Contact')
% 
% pValFsweep = deltaH_ttest(H_Results, 'Frequency')
% 
% 
% 
% %% Final optional saving of figures
% 
% if pipeParams.finalFigs.save == true
%     savfPn = pipeParams.finalFigs.savepn;
%     
%     % save as .eps files
%     saveas(f1, [savfPn, '\dHisiAllBox_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
%     saveas(f2, [savfPn, '\dHisiAllBox_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');
%     saveas(f3, [savfPn, '\dHisiAllBAR_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
%     saveas(f4, [savfPn, '\dHisiAllBAR_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');
% 
%     % save as .jpg files
%     saveas(f1, [savfPn, '\dHisiAllBox_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
%     saveas(f2, [savfPn, '\dHisiAllBox_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');
%     saveas(f3, [savfPn, '\dHisiAllBAR_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
%     saveas(f4, [savfPn, '\dHisiAllBAR_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');
% 
% end
% 
% 



end % END function