% This is the function that all my scripts rely on for paring down tables
% into sub-selections for final analysis. Each case below shows what
% filtering steps are undertaken by each script. Since many of my scripts
% filter down their input tables in similar ways, it made sense to put all
% the individual steps in one place for easy code maintenance. 

% Author: Ed Bello
% Created: 2019/07/04

function FilteredTable = filterTable_Master(tab, callingFunc, ppar)

switch callingFunc
      
    case {'script_displayPhaseLockEntropyChangeVsdbsCond_boxplots', ...
          'script_displayPhaseLockNeuronPerc_pieChart', ...
          'script_displayPhaseLockPercVsdbsTrialCond_barplot', ...
          'script_displayRateVsEntropy_scatter', ...
          'script_displayPercPop_isiAndPhaseLckCategory_pieChart', ...
          'script_displayPopRates_PREandDBS_barplot', ...
          'script_displayRateChangeCategory_barplotStack', ...
          'script_plot3uniqueNeuronLogISIs_freqSwp_PREandDBSon', ...
          'script_displayPopRates_PREandDBS_barplot', ...
          'script_displayPhaseLockEntropyChangeVsdbsCond_boxplots_v30s'}
      
        % ppar.neuronType:
        tab = filter_neuronType(tab, ppar.neuronType);
        
        % ppar.subjID:
        tab = filter_subjID(tab, ppar.subjID);
        
        % ppar.trialType:
        tab = filter_trialType(tab, ppar.trialType);

               % ppar.neuronMinimumTrials = 4; % integer
        tab = filter_neuronMinimumTrials(tab, ppar.neuronMinimumTrials);
        
        tab = addColumn_group(tab, ppar);

    case {'script_displayAverageNeuronEntropyByCondGrp_boxplots', ...
          'script_displayEntropyBitpSecVsFreq_boxplot', ...
          'script_displayEntropyBitpSecVsFreq_PreDBS_boxplot', ...
          'script_displayEntropyBitpSecVsFreq_PosDBS_boxplot', ...
          'script_displayEntropyBitpSecVsFreq_onDBS_boxplot', ...
          'script_displayEntropyBitpSecVsFreqSpecific_boxplot', ...
          'script_displaySpkDataVsDepth', ...
          'script_displayFRindexVsFRbaseline_scatter', ...
          'script_displayEntropyVsLFSorHFSchoose_boxplot', ...
          'script_displayEntropyVsCondition_byNeuron_lines', ...
          'script_displayFRvsFrequencyLogTrans_boxplotAndLines', ...
          'script_displayEntropyBitpSpkVSdbsFreq_boxplot', ...
          'script_displayPercPop_groups', ...
          'script_displayFRvsDBSfrequency_variousFigs', ...
          'script_assessUnitSortQuality', ...
          'script_displayFRvsDBSfrequency_Zscores', ...
          'script_displayFRvsDBSfrequency_FRdbson_bar', ...
          'script_displayEntropyPsthVSdbsFreq_boxplot', ...
          'script_multipleEntropyFigures_compareExamples', ...
          'script_PSTHandISIentropyGroups_pieChart', ...
          'script_displayFRvsDBScontact_FRdbson_bar_v1'}
        
        % ppar.neuronType:
        tab = filter_neuronType(tab, ppar.neuronType);
        
        % ppar.subjID:
        tab = filter_subjID(tab, ppar.subjID);
        
        % ppar.trialType:
        tab = filter_trialType(tab, ppar.trialType);

%                % ppar.neuronMinimumTrials = 4; % integer
%         tab = filter_neuronMinimumTrials(tab, ppar.neuronMinimumTrials);
%         
        tab = addColumn_group(tab, ppar);
   
    otherwise
        error('Calling Function not recognized! Peek inside code for acceptable cases.')
        
        
end

FilteredTable = tab;


end

%% SUB-FUNCTIONS

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

function tabFilt = addColumn_group(tab, ppar)
% add a column of user-defined grouping id according to ppar.groups and its
% sub-fields. If the user has defined ppar.groups, the "groups" column will
% assign user-defined labels to each row element; any row element that was
% not specified in ppar.groups gets called ''. If ppar.groups was not
% defined, then every row element gets called "default". Default behavior
% is that there is no special grouping. 

% so far this is just for the Frequency-sweep case, need to work in a
% Contact-sweep implementation too... for now it's default behavior

group = cell(size(tab, 1), 1);
group(:) = {''};


if isfield(ppar, 'groups')
    switch ppar.trialType
        
        case 'FrequencySweep'
            groupFreqs = ppar.groups.dbsFrequency.groupFreqs;
            groupLabels = ppar.groups.dbsFrequency.groupLabels;

            % assign group 1 label to group column
            nGroups = numel(groupFreqs);

            for iGrp = 1:nGroups
                currGroup = groupFreqs{iGrp};
                isCurrGrp = false(size(group, 1), 1);

                nElements = numel(currGroup);
                for iEl = 1:nElements
                    currGrpElement = currGroup(iEl);
                    isCurrGrp = isCurrGrp | (tab.dbsFrequency == currGrpElement);

                end

                group(isCurrGrp) = groupLabels(iGrp);

            end
            
        case 'ContactSweep'
                group(:) = {'default'};

        otherwise
            error('wrong input for ppar.trialType')
            
    end
    

else
    group(:) = {'default'};
    
end

tabFilt = [tab, table(group)];






end
