% script for converting H_DirectISItable into a strictly single-value table
% element table, that is, no fields in the table hold any arrays within
% them, just one value

% script assumes that you first called "script_pipelineParams"
% save all intermediate data within them as individual .mat files


% re-do this table into two:

%% Initial calls

% load tableRoot:
load([ppar.tablePath, '\', 'tableRoot']);

% load Kramer's EntropyDirectISI_Kramer_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots
load([ppar.tablePath, '\', 'EntropyDirectISI_Kramer_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots']);
Hkr = H_DirectISIResults;

% load Uva's EntropyDirectISI_Kramer_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots
load([ppar.tablePath, '\', 'EntropyDirectISI_Uva_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots']);
Huv = H_DirectISIResults;

H = [Hkr; Huv];



%% 

% Somne constants:
intDataPath = 'DataProcessing\intermediateData\entropyDirISI';
savePath = [ppar.projectRootPath, '\', intDataPath];
metaLabel = '_SU_2HzThresh_60pre_60dbs_ordH2_15binsPD_10000boots';


% prefill columns of future analysis table
        nRecs = size(H, 1);
     objectID = H.objectID;
     
      matFile = cell(nRecs, 1);
matFileFolder = cell(nRecs, 1);

     H_PREemp = zeros(nRecs, 1);
     H_DBSemp = zeros(nRecs, 1);
  H_PREbootAv = zeros(nRecs, 1);
H_PREbootStdv = zeros(nRecs, 1);
         pVal = zeros(nRecs, 1);


% 1) Save data to individual matfiles and 2) track data for table columns
for iRec = 1:nRecs
    % 1)
    % build struct with relevant data and SAVE
     entropyDirISI.preEmp = H.H_PREemp{iRec,1};
    entropyDirISI.preBoot = H.H_PREboot{iRec,1};
     entropyDirISI.dbsEmp = H.H_DBS{iRec,1};

    recLabel = H.objectID{iRec,1};
    fileName = [recLabel, metaLabel, '.mat'];

    save([savePath, '\', fileName], 'entropyDirISI');

    % 2)
    % Add in details to table columns:
    e = entropyDirISI;

          matFile{iRec,1} = fileName;
    matFileFolder{iRec,1} = intDataPath;
         H_PREemp(iRec,1) = e.preEmp(1);
         H_DBSemp(iRec,1) = e.dbsEmp(1);
      H_PREbootAv(iRec,1) = mean(e.preBoot(:,1));
    H_PREbootStdv(iRec,1) = std(e.preBoot(:,1));
             pVal(iRec,1) = H.pVal(iRec,1);
      
end

% build table
% EntropyDirectISI_SU_2HzThresh_60pre_60dbs_ordH1_15binsPD_10000boots
tabName = 'EntropyDirectISI_Analysis_SU_2HzThresh_60pre_60dbs_ordH1_15binsPD_10000boots';
tab = [table(objectID), ... 
    table(matFile), ...
    table(matFileFolder), ...
    table(H_PREemp), ...
    table(H_DBSemp), ...
    table(H_PREbootAv), ...
    table(H_PREbootStdv), ...
    table(pVal)];
    

EntropyDirectISI_analysis = tab;

save([ppar.tablePath, '\', tabName], 'EntropyDirectISI_analysis')


