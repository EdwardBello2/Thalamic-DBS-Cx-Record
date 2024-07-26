function pValTable = deltaH_ttest(ResultsTable, sweepLabel)
% Get paired t-test p-values for each individual DBS condition in the 
% specified Sweep-Type for delta-H measures

% Check to make sure sweepLable is either 'Contact' or 'Frequency'
if ~xor(strcmp(sweepLabel, 'Contact'), strcmp(sweepLabel, 'Frequency'))
    error('Invalid string for sweepLabel');
end

% Choose the rows of the specified sweepType
tableSweeps = ResultsTable.sweepType(:);
isSwp = strcmp(tableSweeps, sweepLabel);
R_swp = ResultsTable(isSwp,:);



switch sweepLabel 
    case 'Contact'
        pValTable = deltaH_ttest_Cswp(R_swp);

    case 'Frequency'
        pValTable = deltaH_ttest_Fswp(R_swp);
        
    otherwise
        disp('No p-values calculated')

end






end

% SUB-FUNCTIONS
function pValTable = deltaH_ttest_Cswp(R_swp)


% FIND unique dbs Contact labels
dbsElec = R_swp.dbsElectrode(:);
sortedDbsElec = sort(dbsElec);
Contact = unique(sortedDbsElec);

pValue = zeros(numel(Contact), 1);

nCX = numel(Contact);
for iCX = 1:nCX
    idxiCX = strcmp(dbsElec, Contact(iCX));
    
    % Perform paired ttest on CO rows delta H
    Hdbs = R_swp.H_DBS(idxiCX);
    Hpre = R_swp.H_PRE(idxiCX);
    [~, pValue(iCX)] = ttest(Hdbs, Hpre); % paired-sample t-test

    
    
    
end

% Create output table
pValTable = table(Contact, pValue);



end

function pValTable = deltaH_ttest_Fswp(R_swp)


% FIND unique dbs Frequency labels
dbsFreq = R_swp.dbsFrequency(:);
sortedDbsFreq = sort(dbsFreq);
Frequency = unique(sortedDbsFreq);

pValue = zeros(numel(Frequency), 1);

nFX = numel(Frequency);
for iFX = 1:nFX
    idxiFX = (dbsFreq == Frequency(iFX));
    
    % Perform paired ttest on CO rows delta H
    Hdbs = R_swp.H_DBS(idxiFX);
    Hpre = R_swp.H_PRE(idxiFX);
    [~, pValue(iFX)] = ttest(Hdbs, Hpre); % paired-sample t-test

    
    
    
end

% Create output table
pValTable = table(Frequency, pValue);



end
