function runEntropyPSTHletter_byPSTHclass(pipeParams);
% Run PSTH-based Entropy Letter-Method by-PSTH classification pipeline
%
% This is like runEntropyPSTHletter, except taking PSTH type into account,
% such as antidromic, orthodromic...
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



%% CREATE TABLE with all data to be analyzed for entropy

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




%% CREATE TABLE with PSTH-Entropy estimates for both pre- and DBS periods

% OPEN the specified dataSelect Table:
load([pipeParams.tablepn, '\', IntTableName, '.mat' ])


% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyPsthLetter';

name = buildEntropyPsthLetterTableName(pipeParams);
AnalysisTableName = [baseLabel, name];


% Check if such a table already exists
if ~exist([pipeParams.tablepn, '\', AnalysisTableName, '.mat']) || ... 
        pipeParams.overwriteTables
    
    disp(['Did not find ', AnalysisTableName, '.mat'])
    disp('Creating...')
    
    % Perform PSTH-Entropy estimate on each row
    [H_letterResults] = calcPSTHletterEntropy_batch(NEX_2analyze, pipeParams); %<---------


    % SAVE final result, as both .mat and .xlsx
    % save the table
    save([pipeParams.tablepn, '\', AnalysisTableName, '.mat'], 'H_letterResults');

%     writetable(H_letterResults, [pipeParams.tablepn, '\', AnalysisTableName, '.xlsx']);
    

    disp('CREATED');

else
    disp(['FOUND ', AnalysisTableName, '.mat'])

end % END if ~exist



%% CREATE TABLE of PSTH-Entropy estimates WITH added column of psthTypes


% OPEN the specified Entropy Results table
load([pipeParams.tablepn, '\', AnalysisTableName, '.mat' ]);
R = H_letterResults; % simplify variable name


% Make name for analysis table depending on pipeParams
baseLabel = 'EntropyPsthLetter';

name = buildEntropyPsthLetterTableName(pipeParams);
AnalPsthTypeName = [baseLabel, name, '_psthType'];



% Check if such a table already exists
if ~exist([pipeParams.tablepn, '\', AnalPsthTypeName, '.mat']) || ... 
        pipeParams.overwriteTables
    
    disp(['Did not find ', AnalPsthTypeName, '.mat'])
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
    save([pipeParams.tablepn, '\', AnalPsthTypeName, '.mat'], 'H_letterResPSTH');


    disp('CREATED');

else
    disp(['FOUND ', AnalPsthTypeName, '.mat'])

end % END if ~exist



%% PLOT stacked bar of %psthType vs DBS condition (for MODULATED trials)

alpha = pipeParams.pValAlpha;

% OPEN the specified Entropy Results table
load([pipeParams.tablepn, '\', AnalPsthTypeName, '.mat' ]);
R = H_letterResPSTH; % simplify variable name


% Get only rows that correspond to modulated trials
isMod = getModIdx(R, alpha);
R_mod = R(isMod, :);


% Divide table into Csweep and Fsweep tables
[R_modCswp, R_modFswp] = divideBySweepType(R_mod);


% Csweep: Get row indices for pshtTypes
[isAntiC, isOrthoC, isInhibC, isOtherC] = findPsthType(R_modCswp);
[Clabels, isClabel] = findDBSlabels_Cswp(R_modCswp);


% Prepare Csweep data for stacked bar
ncLabs = numel(Clabels);
psthType = {'ANTI', 'ORTHO', 'other', 'INHIB'};
for iLab = 1:ncLabs
    % Get numbers of SPECIFIC psth occurrences for each label:
    Cpsth(iLab,1) = sum(isClabel(:,iLab) & isAntiC);
    Cpsth(iLab,2) = sum(isClabel(:,iLab) & isOrthoC);
    Cpsth(iLab,3) = sum(isClabel(:,iLab) & isOtherC);
    Cpsth(iLab,4) = sum(isClabel(:,iLab) & isInhibC);
    
end
CtotPsth = sum(Cpsth,2);
Cpsth_perc = Cpsth ./ CtotPsth;


% Plot stacked bar of psthType vs DBS condition
f1 = figure; 
stackbarCswpPsths(Cpsth_perc, Clabels, psthType)



% Fsweep: Get row indices for pshtTypes
[isAntiF, isOrthoF, isInhibF, isOtherF] = findPsthType(R_modFswp);
[Flabels, isFlabel] = findDBSlabels_Fswp(R_modFswp);


% Prepare Fsweep data for stacked bar
psthType = {'ANTI', 'ORTHO', 'other', 'INHIB'};
nfLabs = numel(Flabels);
for iLab = 1:nfLabs
    % Get numbers of SPECIFIC psth occurrences for each label:
    Fpsth(iLab,1) = sum(isFlabel(:,iLab) & isAntiF);
    Fpsth(iLab,2) = sum(isFlabel(:,iLab) & isOrthoF);
    Fpsth(iLab,3) = sum(isFlabel(:,iLab) & isOtherF);
    Fpsth(iLab,4) = sum(isFlabel(:,iLab) & isInhibF);
    
end
FtotPsth = sum(Fpsth,2);
Fpsth_perc = Fpsth ./ FtotPsth;


% Plot stacked bar of psthType vs DBS condition
f2 = figure; 
stackbarFswpPsths(Fpsth_perc, Flabels, psthType)



%% Final optional saving of figures

if pipeParams.finalFigs.save == true
    savfPn = pipeParams.finalFigs.savepn;
    
    % save as .eps files
    saveas(f1, [savfPn, '\percPSTHtype_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
    saveas(f2, [savfPn, '\percPSTHtype_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');
%     saveas(f3, [savfPn, '\dHletterAll_vs_DBScond_Cswp_', pipeParams.subjID ],'epsc');
%     saveas(f4, [savfPn, '\dHletterAll_vs_DBScond_Fswp_', pipeParams.subjID ],'epsc');

    % save as .jpg files
    saveas(f1, [savfPn, '\percPSTHtype_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
    saveas(f2, [savfPn, '\percPSTHtype_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');
%     saveas(f3, [savfPn, '\dHletterAll_vs_DBScond_Cswp_', pipeParams.subjID ],'jpg');
%     saveas(f4, [savfPn, '\dHletterAll_vs_DBScond_Fswp_', pipeParams.subjID ],'jpg');

end



end % END function



%% SUB-FUNCTIONS

function stackbarFswpPsths(Fpsth_perc, Flabels, psthType)


bar(Fpsth_perc, 'stacked');
set(gca,'xticklabel',Flabels)

title(['Modulated Freq-Sweep trials: % observed PSTH types']);
ylabel('% of modulated trials')
xlabel('DBS Freqeuncy')
legend(psthType);


end

function stackbarCswpPsths(Cpsth_perc, Clabels, psthType)

barh(Cpsth_perc, 'stacked');
set(gca,'yticklabel',Clabels)

title(['Modulated Contact-Sweep trials: % observed PSTH types']);
xlabel('% of modulated trials')
ylabel('DBS electrode')
legend(psthType);


end

function [labels, isLab] = findDBSlabels_Cswp(R_Cswp);



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




end

function [labels, isLab] = findDBSlabels_Fswp(R_Fswp);



% Get positions of each label in Table
labels = unique(R_Fswp.dbsFrequency);
labelTable = R_Fswp.dbsFrequency(:);


% Get total number of trials for each label
nTableRows = size(R_Fswp, 1);
nLabs = numel(labels);

isLab = false(nTableRows, nLabs);
for iLab = 1:nLabs
    isLab(:,iLab) = (labels(iLab) == labelTable);

end



end

function [R_Cswp, R_Fswp] = divideBySweepType(R)

isCswp = strcmp('Contact', R.sweepType(:));
R_Cswp = R(isCswp, :);

isFswp = strcmp('Frequency', R.sweepType(:));
R_Fswp = R(isFswp, :);



end

function [isMod] = getModIdx(R, alpha)
% Get 


% Get idx's of modulated (phase-locked) trials
isMod = R.pVal(:) <= alpha;



end % END function

function [isAnti, isOrtho, isInhib, isOther] = findPsthType(R);

isAnti = strcmp('anti', R.psthType(:)) | ...
         strcmp('anti-ortho', R.psthType(:));

isOrtho = strcmp('ortho', R.psthType(:));

isInhib = strcmp('inhib', R.psthType(:));

isOther = ~(isAnti | isOrtho | isInhib);



end

function Rp = addPSTHcol2Results(R, P)
% Now for given row in REsults (R), find its matching row in UnitPSTHclean 
% (P) and get psthlabel

nRows = size(R, 1);
psthType = cell(nRows, 1);
for iRow = 1:nRows
    iUnit = R.Unit_objectID{iRow,1};
    iElec = R.dbsElectrode{iRow,1};
    iFreq = R.dbsFrequency(iRow,1);
    iSwp = R.sweepType{iRow,1};
    
    isUnitinP = strcmp(iUnit, P.UnitID(:));
    isElecinP = strcmp(iElec, P.DBScontact(:));
    isFreqinP = (P.DBSfrequency(:) == iFreq);
    isSwpinP  = strcmp(iSwp, P.sweepType(:));
    
    % Choose the ONE row in P for which ALL the below is TRUE
    iRowinP = isUnitinP & isElecinP & isFreqinP & isSwpinP;
%     if sum(iRowinP) > 1
%         error(['Row ', num2str(iRow), ' in R is mapped to more than one row in P! Shit!'])
%     end
    
    psthType{iRow,1} = P.psth{iRowinP,1};
    
        
end

Ptype = table(psthType);
Rp = [R, Ptype];

end

function UnitPSTHclean = cleanupUnitPSTH_Kramer(UnitPSTH)
% Remove excess quotes from strings from UnitID and sweepType; 

U = UnitPSTH; % rename for simplicity
Ucl = U; % initialize "clean" table
nRows = size(U, 1);

% Remove excess '' markes from strings
for iRow = 1:nRows
    Ucl.UnitID{iRow,1} = U.UnitID{iRow,1}(2:end-1);

    Ucl.DBScontact{iRow,1} =  U.DBScontact{iRow,1}(2:end-1);
    Ucl.sweepType{iRow,1} =  U.sweepType{iRow,1}(2:end-1);
   
end


UnitPSTHclean = Ucl;



end

function UnitPSTHclean = cleanupUnitPSTH_Uva(UnitPSTH)
% Remove excess quotes from strings; Add column indicating Sweep Type

U = UnitPSTH; % rename for simplicity
Ucl = U; % initialize "clean" table
nRows = size(U, 1);

% Remove excess '' markes from strings
for iRow = 1:nRows
    Ucl.UnitID{iRow,1} = U.UnitID{iRow,1}(2:end-1);
    Ucl.NeuronType{iRow,1} =  U.NeuronType{iRow,1}(2:end-1);
    Ucl.DBScontact{iRow,1} =  U.DBScontact{iRow,1}(2:end-1);
    
    
    
end

% assign proper sweepType labels to rows
sweepType = cell(nRows, 1);

for iRow = 1:nRows
    % if a row has anything but C0 as DBScontact, it's Csweep
    if ~strcmp(Ucl.DBScontact{iRow,1}, 'C0')
        sweepType{iRow,1} = 'Contact';
    
    % if a row has anything but 130Hz as DBSfrequency, it's Fsweep
    elseif Ucl.DBSfrequency(iRow,1) ~= 130
        sweepType{iRow,1} = 'Frequency';
        
    % if a row has C0 and 130Hz assum Csweep, unless...   
    elseif (strcmp(Ucl.DBScontact{iRow,1},'C0')) && ...
           (Ucl.DBSfrequency(iRow,1) == 130)
        sweepType{iRow,1} = 'Contact';
        
        % ...in case of C0 130 belonging to Fsweep, see if previous row was Csweep:
        if (strcmp(Ucl.DBScontact{iRow-1,1},'C0')) && ...
           (Ucl.DBSfrequency(iRow-1,1) == 130)
            
            sweepType{iRow,1} = 'Frequency';
            
        end
      
       
    else
        sweepType{iRow,1} = 'IGNORE';
        
    end
    
end % END for-loop



Sw = table(sweepType);

Ucl = [Ucl, Sw];

UnitPSTHclean = Ucl;



end % END function

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

function dispModPercent_Csweep(Clabels, isClab, isClabMod, alpha)

Ctot = sum(isClab)';
Cmod = sum(isClabMod)';

Cpercent(:,1) = 100 * (Cmod ./ Ctot);
Cpercent(:,2) = 100 * (1 - (Cmod ./ Ctot));


% figure; 
barh(Cpercent, 'stacked');
set(gca,'yticklabel',Clabels)

title(['Percent modulated Contact-sweep trial-population (p<', num2str(alpha), ')']);
xlabel('% of trials')
ylabel('DBS electrode')
legend('Modulated', 'non-Mod');

end

function dispModPercent_Fsweep(Flabels, isFlab, isFlabMod, alpha)


Ftot = sum(isFlab)';
Fmod = sum(isFlabMod)';

Fpercent(:,1) = 100 * (Fmod ./ Ftot);
Fpercent(:,2) = 100 * (1 - (Fmod ./ Ftot));


% figure; 
bar(Fpercent, 'stacked');
set(gca,'xticklabel',Flabels)

title(['Percent modulated Frequency-sweep trial-population (p<', num2str(alpha), ')']);
xlabel('DBS Frequency')
ylabel('% of trials')
legend('Modulated', 'non-Mod');



end

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

