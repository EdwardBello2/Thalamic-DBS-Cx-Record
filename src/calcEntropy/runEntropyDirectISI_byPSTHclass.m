function runEntropyDirectISI_byPSTHclass(pipeParams)
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



%% Create a table with all data to be analyzed for entropy

nexLabel = 'NEXprocfiles';


% Create a new name for the intermediate table based on the pipeline params
name = buildDataSelectTableName(pipeParams);
IntTableName = [nexLabel, 'Selected', name];


% Check if the intermediate table to be created already exists. If it does
% not exist, runPipeline will proceed to create the table

if ~exist([pipeParams.tablepn, '\', IntTableName, '.mat']) || ...
        pipeParams.overwriteTables
    
    disp(['Did not find ', IntTableName, '.mat'])
    disp('Creating...')

    
    % Create a table of selected rows
    Selected = createSelectedDataTable(NEX, Sort, pipeParams);

    % Match each row with DBS Contact-Sweep or Freq-Sweep info, and remove any
    % rows that pertain to neither
    NEX_2analyze = updateNEXprocfiles_subjID_4isiEntropyAnalysis(Selected, TrialInfo);


    % SAVE final result, as both .mat and .xlsx
    save([pipeParams.tablepn, '\', IntTableName, '.mat'], 'NEX_2analyze');

    writetable(NEX_2analyze, [pipeParams.tablepn, '\', IntTableName, '.xlsx'])



    disp('CREATED')
    
else
    disp(['FOUND ', IntTableName, '.mat'])

end % END if ~exist



%% CREATE a table with Direct-Entropy estimates for both pre- and DBS periods

% OPEN the specified dataSelect Table:
load([pipeParams.tablepn, '\', IntTableName, '.mat' ]) 


% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyDirectISI';
name = buildEntropyAnalyzeTableName(pipeParams);
EntroyDirectTableName = [baseLabel, name];


% Check if such a table already exists
if ~exist([pipeParams.tablepn, '\', EntroyDirectTableName, '.mat']) || ...
        pipeParams.overwriteTables
    
    disp(['Did not find ', EntroyDirectTableName, '.mat'])
    disp('Creating...')
    
    % Perform Direct-Entropy estimate on each row
    [H_Results] = calcDirectEntropyISI_batch(NEX_2analyze, pipeParams); %<---------


    % SAVE final result, as both .mat and .xlsx
    % save the table
    save([pipeParams.tablepn, '\', EntroyDirectTableName, '.mat'], 'H_Results');

    writetable(H_Results, [pipeParams.tablepn, '\', EntroyDirectTableName, '.xlsx']);
    

    disp('CREATED');

else
    disp(['FOUND ', EntroyDirectTableName, '.mat'])

end % END if ~exist



%% CREATE TABLE with PSTH-Entropy estimates for both pre- and DBS periods

% OPEN the specified dataSelect Table:
load([pipeParams.tablepn, '\', IntTableName, '.mat' ])


% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyPsthLetter';

name = buildEntropyPsthLetterTableName(pipeParams);
EntroyLetterTableName = [baseLabel, name];


% Check if such a table already exists
if ~exist([pipeParams.tablepn, '\', EntroyLetterTableName, '.mat']) || ...
        pipeParams.overwriteTables
    
    disp(['Did not find ', EntroyLetterTableName, '.mat'])
    disp('Creating...')
    
    % Perform PSTH-Entropy estimate on each row
    [H_letterResults] = calcPSTHletterEntropy_batch(NEX_2analyze, pipeParams); %<---------


    % SAVE final result, as both .mat and .xlsx
    % save the table
    save([pipeParams.tablepn, '\', EntroyLetterTableName, '.mat'], 'H_letterResults');

%     writetable(H_letterResults, [pipeParams.tablepn, '\', AnalysisTableName, '.xlsx']);
    

    disp('CREATED');

else
    disp(['FOUND ', EntroyLetterTableName, '.mat'])

end % END if ~exist



%% CREATE TABLE of PSTH-Entropy estimates WITH added column of psthTypes


% OPEN the specified Entropy Results table
load([pipeParams.tablepn, '\', EntroyLetterTableName, '.mat' ]);
R = H_letterResults; % simplify variable name


% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyPsthLetter';

name = buildEntropyPsthLetterTableName(pipeParams);
EntropyLetterPsthTableName = [baseLabel, name, '_psthType'];



% Check if such a table already exists
if ~exist([pipeParams.tablepn, '\', EntropyLetterPsthTableName, '.mat']) || ...
        pipeParams.overwriteTables
    
    disp(['Did not find ', EntropyLetterPsthTableName, '.mat'])
    disp('Creating...')
    
    % Load in data with PSTH classifications for each trial
    UnitPSTH = readtable([pipeParams.tablepn, '\', 'Units_PSTHlabels_', pipeParams.subjID, '.xlsx']);


    % Clean up string entries; Add column indicating Sweep Type
    switch pipeParams.subjID
        case 'Uva'
            UnitPSTHclean = cleanupUnitPSTH_Uva(UnitPSTH);
            
        case 'Kramer'
            UnitPSTHclean = cleanupUnitPSTH_Kramer(UnitPSTH);
            
        otherwise
            error(['Need to code up a Units_PSTHlabels read-in strategy for subjID: ', pipeParams.subjID])
    end

    % Add column of psthType labels to H_letterResults
    H_letterResPSTH = addPSTHcol2Results(H_letterResults, UnitPSTHclean);


    % SAVE final result as .mat table only
    save([pipeParams.tablepn, '\', EntropyLetterPsthTableName, '.mat'], 'H_letterResPSTH');


    disp('CREATED');

else
    disp(['FOUND ', EntropyLetterPsthTableName, '.mat'])

end % END if ~exist



%%  OPEN the specified ISI-Entropy Results table and Phase-lock Entropy results table
% Combine info from them

% Loads as "H_Results"
load([pipeParams.tablepn, '\', EntroyDirectTableName, '.mat' ]);
H_ResD = H_Results;
% Loads as "H_letterResPSTH"
load([pipeParams.tablepn, '\', EntropyLetterPsthTableName, '.mat' ]);
H_ResL = H_letterResPSTH;


% ADD psthType colum to ISI-entropy table from corresponding rows in PSTH
% table:

% find row in L that corresponds to row in D
nRowsD = size(H_ResD, 1);
psthType = cell(nRowsD, 1);
for iRowD = 1:nRowsD
    isRowDinL = strcmp(H_ResD.objectID{iRowD}, H_ResL.objectID(:));
    psthType{iRowD,1} = H_ResL.psthType{isRowDinL};
    
end
PSTH = table(psthType);

H_ResDpsth = [H_ResD, PSTH];


% Divide table into groupings based on PSTH type
[isAnti, isOrtho, isInhib, isOther] = findPsthType(H_ResDpsth);
H_ResDanti = H_ResDpsth(isAnti,:);
H_ResDortho = H_ResDpsth(isOrtho,:);
H_ResDinhib = H_ResDpsth(isInhib,:);
H_ResDother = H_ResDpsth(isOther,:);






%% PLOT boxplots of entropy results for both sweeps-types

% 
% % Plot the deltaHs of both sweep-types
% deltaH_boxplots(H_Results)

f1 = figure;
f1.Position = [2044 520 771 420];
barplot_deltaH_Cswp_allPsthType(H_ResDpsth)
yLim = get(gca, 'YLim');
yLim(1) = pipeParams.maxdH;
set(gca, 'YLim', yLim);

f2 = figure;
f2.Position = [2044 520 771 420];
barplot_deltaH_Fswp_allPsthType(H_ResDpsth)
yLim = get(gca, 'YLim');
yLim(1) = pipeParams.maxdH;
set(gca, 'YLim', yLim);

% switch pipeParams.dispPsthType
%     
%     case 'anti'   
%         f1 = figure;
%         boxplot_deltaH_Cswp(H_ResDanti)
% 
%         f2 = figure;
%         boxplot_deltaH_Fswp(H_ResDanti)
% 
%         f3 = figure;
%         barplot_deltaH_Cswp(H_ResDanti);
% 
%         f4 = figure; 
%         barplot_deltaH_Fswp(H_ResDanti);
%         
%     case 'ortho'
%         f1 = figure;
%         boxplot_deltaH_Cswp(H_ResDortho)
% 
%         f2 = figure;
%         boxplot_deltaH_Fswp(H_ResDortho)
% 
%         f3 = figure;
%         barplot_deltaH_Cswp(H_ResDortho);
% 
%         f4 = figure; 
%         barplot_deltaH_Fswp(H_ResDortho);
%         
%     case 'inhib'
%         f1 = figure;
%         boxplot_deltaH_Cswp(H_ResDinhib)
% 
%         f2 = figure;
%         boxplot_deltaH_Fswp(H_ResDinhib)
% 
%         f3 = figure;
%         barplot_deltaH_Cswp(H_ResDinhib);
% 
%         f4 = figure; 
%         barplot_deltaH_Fswp(H_ResDinhib);
%         
%     case 'other'
%         f1 = figure;
%         boxplot_deltaH_Cswp(H_ResDother)
% 
%         f2 = figure;
%         boxplot_deltaH_Fswp(H_ResDother)
% 
%         f3 = figure;
%         barplot_deltaH_Cswp(H_ResDother);
% 
%         f4 = figure; 
%         barplot_deltaH_Fswp(H_ResDother);
%         
%     otherwise
%         error('pipeParams.dispPsthType must be one of these strings: anti, ortho, inhib, other');
%         
% end



% % Get paired-sample T-test results for each condition in Contact-sweep
% 
% pValCsweep = deltaH_ttest(H_Results, 'Contact')
% 
% pValFsweep = deltaH_ttest(H_Results, 'Frequency')



%% Final optional saving of figures

if pipeParams.finalFigs.save == true
    savfPn = pipeParams.finalFigs.savepn;
    
    % save as .eps files
%     saveas(f1, [savfPn, '\dHisiAllBox_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
%     saveas(f2, [savfPn, '\dHisiAllBox_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');
    saveas(f1, [savfPn, '\dHisiPsthBAR_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f2, [savfPn, '\dHisiPsthBAR_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');

    % save as .jpg files
%     saveas(f1, [savfPn, '\dHisiAllBox_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
%     saveas(f2, [savfPn, '\dHisiAllBox_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');
    saveas(f1, [savfPn, '\dHisiPsthBAR_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f2, [savfPn, '\dHisiPsthBAR_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');

end





end % END function