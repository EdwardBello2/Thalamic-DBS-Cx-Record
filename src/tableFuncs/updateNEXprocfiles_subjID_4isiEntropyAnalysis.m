function Updated = updateNEXprocfiles_subjID_4isiEntropyAnalysis(NexTable, TrialTable)
% Adds dbs parameter details to each entry, using the second input table. 
% Also removes any rows that do not correspond to either a Contact-sweep or
% Frequency-sweep trial

% Input: NexTable, TrialTable


N = NexTable;
nRows = size(N,1);

T = TrialTable;




NtrialIDs = N.Trial_objectID(:);
allTrialIDs = T.objectID(:);

% Create the three columns of interest first, then merge them to table
dbsElectrode = cell(nRows, 1);
dbsFrequency = zeros(nRows, 1);
sweepType = cell(nRows, 1);


for iRow = 1:nRows
    
    % find the trial in T that corresponds to the current row in N
    isTrialiRow = strcmp(NtrialIDs{iRow,1}, allTrialIDs);

    
    dbsElectrode{iRow,1} = T.dbsContact{isTrialiRow,1};
    
    dbsFrequency(iRow,1) = T.dbsFrequency(isTrialiRow,1);
    
    
    % sort out which sweep type it is -- a 1 or 0 in the T table indicates
    % identity, not one of my better ideas...
    if (T.ContactSweep(isTrialiRow,1)) && ...
      ~(T.FrequencySweep(isTrialiRow,1))
  
        sweepType{iRow,1} = 'Contact';
        
    elseif ~(T.ContactSweep(isTrialiRow,1)) && ...
            (T.FrequencySweep(isTrialiRow,1))
        
        sweepType{iRow,1} = 'Frequency';
        
    else
        sweepType{iRow,1} = 'REMOVE';
        
    end
    
end


dbsE = table(dbsElectrode);
dbsF = table(dbsFrequency);
swpT = table(sweepType);

Updated = [N, dbsE, dbsF, swpT];


% final check for rows that should be removed because they do not have
% either label
isRemove = strcmp('REMOVE', sweepType);


Updated(isRemove,:) = [];




end