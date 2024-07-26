% test script for genIntData_spkRate_XXsec_XsecBins



pipeParams.tablepn = 'C:\Users\bello043\Documents\GitHub\Thalamic-DBS-Cx-Record\tables';
pipeParams.subjID = 'Kramer';

% define time of each bin, and standard total time to divide pre/dbs/pos
% into:
pipeParams.ratebinWidth = 1; % seconds
pipeParams.totTime = 60; % seconds

% directory to store all intermediate data files to:


pipeParams.intDataPn = ['L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\intermediateData\spkRate\', pipeParams.subjID];






% load NEXprocfiles_XXX table, where metadata for individual NEX files is
% stored:
load([pipeParams.tablepn, '\', 'NEXprocfiles_', pipeParams.subjID, '.mat']);
NEX = NEXprocfiles; 


% generate Intermediate data table:
genTable_intDataRateBins_XXsec_XsecBins_subjID(NEX, pipeParams);