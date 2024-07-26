% script to analyze Direct-ISI entropy changes according to Phase-locking
% behavior

% Created: 2019/11/22
% based on script_displayPhaseLockEntropyChangeVsdbsCond_boxplots.m, except
% calling the table
% "EntropyDirectISI_Analysis_SU_2hzThresh_30pre_30dbs_ordH1_15binsPD_1000boots.mat"
% instead of
% "EntropyDirectISI_Analysis_SU_2HzThresh_60pre_60dbs_ordH1_15binsPD_10000boots.mat"
% for entropy values to include


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
% 
% % Final display parameters
% ppar.dispPhaseLocked = 'phaseLck'; % 'phaseLck' | 'nonPhaseLck' | 'both'
clear; 

% Initialize fields in ppar struct in an accompanying script:
script_pipelineParams

 % CONSTANTS
CURR_FUNC = 'script_displayPhaseLockEntropyChangeVsdbsCond_boxplots_v30s'; 
FIG_POSITION = [14 59 560 420];


%% LOAD all relevant tables and MERGE them

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



%% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);



%% Detect those individual trials RECORDS in which neurons phase-locked 
% phase-locked must have a p-value less than ppar.pAlphaPSTH
% detection gets added as an extra column of information about each RECORD

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

% Tally unique list of neurons with phase-locking in new table
isPhsLck = Tselect.PhaseLocked;
P = Tselect(isPhsLck,:);
neuPhsLck = unique(P.Unit_objectID);

 isPhsLckNeuron = true(size(neuPhsLck, 1), 1);
phsLckTab = [table(neuPhsLck), table(isPhsLckNeuron)];

NphsLckInfo = outerjoin(Tselect, phsLckTab, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'neuPhsLck');

isPhsLckNeur = NphsLckInfo.isPhsLckNeuron;


% For final results, track number of data points 



%% Create User-defined Group column, with TF indices for which rows have data that belongs to said group

% ppar.groups.from.vars = {'dbsFrequency'};
% ppar.groups.from.select = {10};
% ppar.groups.newLabel = {'Hz10'};
if isfield(ppar, 'groups')
    
    if isfield(ppar.groups, 'from') & isfield(ppar.groups, 'newLabel')
        % Get data from user-specified existing column:
        varsData = Tselect{:, ppar.groups.from.vars{1}};

        % find user-defined selection from vars in table

        isUserGrp = false(size(varsData, 1), 1);
        for i = 1:numel(ppar.groups.from.select)
            isUserGrp = isUserGrp | (varsData == ppar.groups.from.select{i});

        end

        Grp = table(isUserGrp);
        Grp.Properties.VariableNames{1} = ppar.groups.newLabel{1};


        Tselect = [Tselect, Grp];

    end

end


%% CALCULATE and PLOT phase-locked neurons' Direct-ISI entropy changes

switch ppar.dispPhaseLocked
    
    case 'phaseLck'
        Tfinal = Tselect(isPhsLckNeur,:);
        patternTitle = 'phaselockers';
        
    case 'nonPhaseLck'
        Tfinal = Tselect(~isPhsLckNeur,:);
        patternTitle = 'non-phaselockers';

    case 'both'
        Tfinal = Tselect;
        patternTitle = 'all patterns';
        
    otherwise
        error('wrong string input for ppar.dispPhaseLocked');
 
        
end



%% Apply log-transform according to ppar.logTransform
% The +1 operation is to protect against the bias that would happen for any
% entropy values less that 1.0 (log(n<1) = "a really negative number")
% 'H_DBSemp_Hisi' | 'H_PREbootAv' | 'finalEntropy'
% If logTransform is not a field (not specified), then skip this step
if isfield(ppar, 'logTransform')
    if any(strcmp(ppar.logTransform, 'H_DBSemp_Hisi'))
        Tfinal.H_DBSemp_Hisi(:) = log(Tfinal.H_DBSemp_Hisi(:));
        
    end
        
    if any(strcmp(ppar.logTransform, 'H_PREbootAv'))
        Tfinal.H_PREbootAv(:) = log(Tfinal.H_PREbootAv(:));
        
    end        
    
end



%% Calculate final Entropy analysis based on ppar.entropyType

% Collect Entropy values for comparison
switch ppar.entropyType
    case 'isiDbsH'
%             EntropyChange_AllNeu = TbyGroup.H_DBSemp_Hisi;
        EntropyChange_AllNeu = Tfinal.H_DBSemp_Hisi;
        entropyTitle = 'DBS-Entropy';
        entropyAxisLabel = 'bits/spike';
        
     case 'isiPreH'
%             EntropyChange_AllNeu = TbyGroup.H_DBSemp_Hisi;
        EntropyChange_AllNeu = Tfinal.H_PREbootAv;
        entropyTitle = 'baseline-Entropy';
        entropyAxisLabel = 'bits/spike';

    case 'isiDiffH'
%             EntropyChange_AllNeu = TbyGroup.H_DBSemp_Hisi - TbyGroup.H_PREbootAv;
        EntropyChange_AllNeu = Tfinal.H_DBSemp_Hisi - Tfinal.H_PREbootAv;
        entropyTitle = '\DeltaEntropy';
        entropyAxisLabel = '\Delta bits/spike';

    case 'isiPercDiffH'
%             EntropyChange_AllNeu = (TbyGroup.H_DBSemp_Hisi - TbyGroup.H_PREbootAv) ./ ...
%                         TbyGroup.H_PREbootAv;
        EntropyChange_AllNeu = (Tfinal.H_DBSemp_Hisi - Tfinal.H_PREbootAv_Hisi) ./ ...
                    Tfinal.H_PREbootAv_Hisi;
        entropyTitle = '%\DeltaEntropy';
        entropyAxisLabel = '%\Delta Entropy';

    otherwise
        error('wrong input for ppar.entropyType');

end

% log-transform final entropy data results if specified in
% ppar.logTransform:
if isfield(ppar, 'logTransform')
    if any(strcmp(ppar.logTransform, 'finalH'))
        % control for negative difference cases
        nRows = size(EntropyChange_AllNeu, 1);
        
        for iRow = 1:nRows
%             disp('old value: ')
%             EntropyChange_AllNeu(iRow)
            
            if EntropyChange_AllNeu(iRow) < 0
                EntropyChange_AllNeu(iRow) = -log(abs(EntropyChange_AllNeu(iRow)) + 1);
                
            elseif EntropyChange_AllNeu(iRow) == 0
                EntropyChange_AllNeu(iRow) = 0;
                
            else
                EntropyChange_AllNeu(iRow) = log(EntropyChange_AllNeu(iRow) + 1);
            
            end
            
%             disp('new value: ')
%             EntropyChange_AllNeu(iRow)
            
        end
        
    end
    
end

% H1dbs = Tfinal.H_DBSemp_Hisi;
% H1pre = Tfinal.H_PREbootAv;




% % Show deltaH normalized by preDBS Entropy:
% EntropyChange_AllNeu = 100*((H1dbs - H1pre) ./ H1pre);  



%% Display boxplots & ANOVA results of ALL cells 

f1 = figure;
ax1 = axes;
switch ppar.trialType
    
    case 'FrequencySweep'
        % Get Frequency labels of all neurons as array of strings
        if isfield(ppar, 'groups')
            if isfield(ppar.groups, 'from')
                % get sub-selection of group labels and entropy values 
                % according to user-grouping
                idxGrp = Tfinal{:, ppar.groups.newLabel{1}};
                freqs = Tfinal{idxGrp, 'dbsFrequency'};

                % Whittle down the Entropy results to user-defined subselection
                EntropyChange_AllNeu = EntropyChange_AllNeu(idxGrp);
                Entropy_grpLabel = cell(size(EntropyChange_AllNeu, 1), 1);
                Entropy_grpLabel(:) = {ppar.groups.newLabel{1}};
            
            else
                freqs = Tfinal{:, 'dbsFrequency'};
                Entropy_grpLabel = cellstr(num2str(freqs));
            
            end
 
        end
        
        Flabels_AllNeu = unique(Entropy_grpLabel);
        
        dispDeltaH_Fswp_Boxplot(EntropyChange_AllNeu, Entropy_grpLabel, ...
                                Flabels_AllNeu, ppar, 'Pres', 'cleanWoutliers');                        
        ylabel(entropyAxisLabel)

    case 'ContactSweep'
        % Get Electrode Contact labels for all neurons as array of strings
        if isfield(ppar, 'groups')
            % get sub-selection of group labels and entropy values 
            % according to user-grouping
            idxGrp = Tfinal{:, ppar.groups.newLabel{1}};
            contacts = Tfinal{idxGrp, 'dbsContact'};
            
            % Whittle down the Entropy results to user-defined subselection
            EntropyChange_AllNeu = EntropyChange_AllNeu(idxGrp);
                       
        else
            contacts = Tfinal{:, 'dbsContact'};
             
        end
        
        Entropy_grpLabel = contacts;
        Clabels_AllNeu = unique(contacts);
        
        dispDeltaH_Cswp_Boxplot(EntropyChange_AllNeu, Entropy_grpLabel, ...
                                Clabels_AllNeu, 'Pres', 'cleanWoutliers');
        xlabel(entropyAxisLabel)
        
    otherwise
        error('Oh snap! wrong trialType')
        
end

f1.Position = FIG_POSITION;
t = title([ppar.subjID, ':  ', entropyTitle, ' for ', ...
            patternTitle]); 



%% Statistical Test on data

if ~isfield(ppar, 'statTest'), ppar.statTest = 'anova'; end

switch ppar.statTest
    
    case 'anova'
        disp('ANOVA results:')
        [pF, tabF, statsF] = anova1(EntropyChange_AllNeu, Entropy_grpLabel)

    case 'kruskalwallis'
        disp('KRUSKAL-WALLIS results:')
        [pV, tab, stats] = kruskalwallis(EntropyChange_AllNeu, Entropy_grpLabel)

    case 'ttest'
        disp('T-TEST results:')
        [h, pV, ci, stats] = ttest(EntropyChange_AllNeu)
        statTable = table([stats.tstat; stats.df; stats.sd; pV]);
        statTable.Properties.RowNames = {'tstat', 'df', 'sd', 'p'};
        
    otherwise
        error('Oh snap! wrong string input for statTest')
        
end

% Final table summarizing the number of recordings for each condition:


% Final "deltaH_cell" detailing data points according to group order in "label_cell" 
labels_cell= unique(Entropy_grpLabel);
 
nLabs = numel(labels_cell);
for iLab = 1:nLabs
    isLabel = strcmp(labels_cell{iLab}, Entropy_grpLabel);
    deltaH_cell{iLab} = EntropyChange_AllNeu(isLabel);
    
end



% %% Now group LFS and HFS
% 
% switch ppar.trialType
%     
%     case 'FrequencySweep'
%         % Choose the two DBS frequency labels from the data for comparison
%         % strings: 10 | 20 | 30 | 50 | 100 | 130
%         LFSfreq1 = ' 10'; 
%         LFSfreq2 = ' 20';
% 
%         HFSfreq1 = '100';
%         HFSfreq2 = '130';
% 
% 
%         % Sift out all frequencies but the two choices
%         isLFS1 = strcmp(LFSfreq1, Entropy_grpLabel);
%         isLFS2 = strcmp(LFSfreq2, Entropy_grpLabel);
% 
%         isHFS1 = strcmp(HFSfreq1, Entropy_grpLabel);
%         isHFS2 = strcmp(HFSfreq2, Entropy_grpLabel);
% 
%         % isInclude = isLFS1 | isLFS2 | isHFS1 | isHFS2;
%         % 
%         % H_include = Entropy_AllNeu(isInclude);
%         % Freq_include = Entropy_grpLabel(isInclude);
% 
% 
%         % combine lows and combine highs
%         isLFS = isLFS1 | isLFS2;
%         isHFS = isHFS1 | isHFS2;
% 
% 
%         % Divide Entropy data between the two groups
%         Hlfs = EntropyChange_AllNeu(isLFS);
%         Hhfs = EntropyChange_AllNeu(isHFS);
% 
%         [hTtest, pTtest, ci, statsTtest] = ttest2(Hlfs, Hhfs)
% 
%         [pRanksum, hRanksum, statsRanksum] = ranksum(Hlfs, Hhfs)
%         
%     otherwise % do nothing, I don't have this figured out for ContactSweep case...
%         
% end
% 

