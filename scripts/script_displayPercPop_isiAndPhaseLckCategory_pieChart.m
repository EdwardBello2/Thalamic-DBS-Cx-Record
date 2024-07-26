% Display the population percentage of phase-locked cells, according to
% dbs-trial condition (% at each frequency | % at each Contact)

% Author: Ed Bello
% Created: 2019/07/11

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
CURR_FUNC = 'script_displayPercPop_isiAndPhaseLckCategory_pieChart'; 



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



%% Detect those neuron trial-observations that either had entropy increase or decrease

% Add in entropy increase/decrease TF values for each row
nRows = size(Tselect, 1);
isiH_incr = false(nRows, 1);
isiH_decr = false(nRows, 1);

for iRec = 1:nRows
    if Tselect.pVal_Hisi(iRec) <= ppar.pAlphaISI
        
        if Tselect.H_DBSemp_Hisi(iRec) > Tselect.H_PREbootAv(iRec)
            isiH_incr(iRec) = true;
            
        elseif Tselect.H_DBSemp_Hisi(iRec) < Tselect.H_PREbootAv(iRec)
            isiH_decr(iRec) = true;
        
        else % in this case, both columns will have a false
            
        end
        
    end % if not below pAlphaISI, both columns will have a false
    
end

Tselect = [Tselect, table(isiH_incr), table(isiH_decr)];



%% Assign Category label to each row according to the TF values of:

% Columns correspond to whether or not: 
% 1) Phase-locked, 2) isiH increased, 3) isiH decreased
% Rows corresponed to pieLabels below, in order

% cType = [false, false, false;
%          false, false, true;
%          true, false, true;
%          true, false, false;
%          true, true, false;
%          false, true, false;];
% 
% % I've observed 6 possible categories in this data
% pieLabels = {'noChange', 'noPhsLck_Hdecr', 'noPhsLck_Hincr', ...
%              'PhsLckONLY', 'PhsLck_Hdecr', 'PhsLck_Hincr'};
         
cType = [false, false, false;
         false, false, true;
         true, true, false;
         true, false, false;
         true, false, true;
         false, true, false;];

% I've observed 6 possible categories in this data
pieLabels = {'noChange', 'noPhsLck_Hdecr', 'PhsLck_Hdecr', ...
             'PhsLckONLY', 'PhsLck_Hincr', 'noPhsLck_Hincr'};

nRows = size(Tselect, 1);  
pieCategory = cell(nRows, 1);
for iRow = 1:nRows
    if Tselect.PhaseLocked(iRow) % PhsLckONLY, PhsLck_Hdecr, PhsLck_Hincr
        
        if Tselect.isiH_incr(iRow)
            pieCategory{iRow} = 'PhsLck_Hincr';
            
        elseif Tselect.isiH_decr(iRow)
            pieCategory{iRow} = 'PhsLck_Hdecr';
            
        else
            pieCategory{iRow} = 'PhsLckONLY';
            
        end
        
        
    else % noChange, noPhsLck_Hdecr, noPhsLck_Hincr
       
        if Tselect.isiH_incr(iRow)
            pieCategory{iRow} = 'noPhsLck_Hincr';
            
        elseif Tselect.isiH_decr(iRow)
            pieCategory{iRow} = 'noPhsLck_Hdecr';
            
        else
            pieCategory{iRow} = 'noChange';
            
        end 
        
    end
    
end

Tselect = [Tselect, table(pieCategory)];


%% TALLY all Trial Observation counts for each category below
% 

% how many groups -- note that '' is the default empty case, to be ignored
% grpLabels = unique(Tselect.group);
% grpLabels(strcmp(grpLabels, '')) = [];


if isfield(ppar, 'groups')
    grpLabels = ppar.groups.dbsFrequency.groupLabels;
    
else
    grpLabels = unique(Tselect.group);
    
end

     
% Gather up data for first group first, then next, etc...
nGroups = numel(grpLabels);
for iGrp = 1:nGroups
    isCurrGrp = strcmp(Tselect.group, grpLabels{iGrp});
    TbyGroup = Tselect(isCurrGrp, :);
    
    nObs = size(TbyGroup, 1);
        
    nPieLabels = numel(pieLabels);
    for iPiLab = 1:nPieLabels % count num obs for each category
        catCount{iGrp}(iPiLab) = sum(strcmp(TbyGroup.pieCategory, pieLabels(iPiLab)));
        pieData{iGrp}(iPiLab) = round(100 * catCount{iGrp}(iPiLab) / nObs);

    end
         
end




%% PLOT piecharts according to user-defined groupings

f1 = figure; 
% ax1 = axes;


nGroups = numel(grpLabels);
if nGroups > 1
    
    for iGrp = 1:nGroups
        ax(iGrp) = subplot(1, nGroups, iGrp);
        pie1 = pie(pieData{iGrp}, pieLabels);
        title([ppar.subjID, ':   ', ppar.trialType, ', ', grpLabels{iGrp}])

    end

else % case where there is only one group, or default case of no group
    pie1 = pie(pieData{1}, pieLabels);
    
end
    



%% DISPLAY table with counts according to user-defined grouping
sz = [2, 6];
varTypes = {'double', 'double', 'double', 'double', 'double', 'double'};

% create a cell array of N tables
grpTables = cell(nGroups, 1);
grpTablesLabel = grpLabels;

for iGrp = 1:nGroups
    NeuronCounts = table('Size', sz, 'VariableTypes', varTypes);
    NeuronCounts.Variables = [catCount{iGrp}; pieData{iGrp}];
    NeuronCounts.Properties.VariableNames = pieLabels;
    NeuronCounts.Properties.RowNames = {'Count', 'Percent'};

    totNum = sum(catCount{iGrp});
    totalObservations = [totNum; 100];

    NeuronCounts = [NeuronCounts, table(totalObservations)];
    grpTables{iGrp} = NeuronCounts;

    % display in command line:
    grpTablesLabel{iGrp}
    grpTables{iGrp}
       
    
end


