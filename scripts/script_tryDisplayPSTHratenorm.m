% test script for an example cell EPF over frequencies





CURR_FUNC = 'script_displayPopRates_PREandDBS_barplot';
%% LOAD all relevant tables and MERGE them

close all

script_pipelineParams

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);

isFiveNeurons = Tselect.neuCount == 5;
Tselect = Tselect(isFiveNeurons,:);

% Detect those individual trials RECORDS in which neurons phase-locked 
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



% DETECT which neurons had at least one phase-lock reaction to DBS 

% go thru each detected phase-lock, find its Unit_objectID, mark every row
% for that unit at PhaseLockNeuron = 1

idxPhsLck = find(Tselect.PhaseLocked);
nPhsLckRows = numel(idxPhsLck);
isPhsLckNeuron = false(size(Tselect, 1) , 1);

for i = 1:nPhsLckRows 
    iRow = idxPhsLck(i);
        
    % find the Unit_objectID belonging to this row
    currentNeuron = Tselect.Unit_objectID{iRow};
    
    % find all rows that pertain to the current neuron
    idxCurrNeuron = strcmp(currentNeuron, Tselect.Unit_objectID);
    
    % update column tracking all neurons showing at least one phase-lock
    isPhsLckNeuron = isPhsLckNeuron | idxCurrNeuron;
    
end

Tselect = [Tselect, table(isPhsLckNeuron)];

% % Tally unique list of neurons with phase-locking in new table
% isPhsLck = Tselect.PhaseLocked;
% P = Tselect(isPhsLck,:);
% neuPhsLck = unique(P.Unit_objectID);
% 
%  isPhsLckNeuron = true(size(neuPhsLck, 1), 1);
% phsLckTab = [table(neuPhsLck), table(isPhsLckNeuron)];
% 
% NphsLckInfo = outerjoin(Tselect, phsLckTab, 'LeftKeys', 'Unit_objectID', 'RightKeys', 'neuPhsLck');
% 
% isPhsLckNeur = NphsLckInfo.isPhsLckNeuron;
% 
% Tselect = [Tselect, table(isPhsLckNeur)];
% 
% Tselect = Tselect(isPhsLckNeur,:);
% % For final results, track number of data points 

Tselect(~Tselect.isPhsLckNeuron,:) = [];

unique(Tselect.Unit_objectID)


%% load table subselection for nr17080401-ch12a
isN = strcmp(Tselect.Unit_objectID, 'nr17080401-ch03a');
Tfinal = Tselect(isN,:);


nRows = size(Tfinal, 1);
for iRec = 1:nRows

% % for one NEX file:

nexFile = readNexFile([ppar.projRootPath, '\', Tfinal.nexFileFolder{iRec}, ...
        '\', Tfinal.nexFile{iRec}]);



disp([num2str(Tfinal.dbsFrequency(iRec)), ' Hz'])
[eEPF, pfs, pfsb] = calcEEPF(nexFile)
title([num2str(Tfinal.dbsFrequency(iRec)), ' Hz, eEPF: ', num2str(eEPF)])

end