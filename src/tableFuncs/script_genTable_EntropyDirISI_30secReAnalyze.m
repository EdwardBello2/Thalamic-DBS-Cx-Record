% script for calculating Direct Entropy estimate of spike activity,
% bootstrapping the pre-dbs condition to see if the dbs condition is
% significantly different. After getting reviews back for VPLo manuscript,
% I determined to re-do the Entropy analysis based only on 30-sec preDBS,
% and from seconds 30-60 during DBS (in order to omit any initial dynamic
% changes immediately after DBS onset; we want to look at something like a
% steady-state Entropy during DBS...). This script is based on my previous
% script "script_genTable_EntropyDirISI.m", with some changes to account
% for the above mentioned concerns.
%
% function [Hdir_ISI_Results] = calcDirectEntropyISIbootstrap_batch(N, pipeParams)

clear

%% Pipeline Parameters

ppar.projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record'; % string

% full path on local PC where tables to be loaded are kept
% ppar.tablePath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string
ppar.tablePath = 'C:\Users\bello043\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string


% CONSTANTS
% in seconds. Time to include spikes for analysis before DBS onset
PREDBS_TIME = 30; % seconds, maximum time before DBS to get spikes from 
DBS_TIME = 30; % seconds, total time to get spikes from after first DBS pulse
SPARSE_HZ = 2; % Any cells that have av spkRate below this are excluded

MATFILESAVEPATH = 'DataProcessing\intermediateData\entropyDirISI_30s';

% PARAMETERS for Direct Entropy estimate with log binned ISI:
ORD_H = 1; % how many orders of entropy to use for linear extrapolation
BINS_PER_DECADE = 15; % of log-spaced bins for ISI histogram
NBOOTS = 10000; % integer, number of bootstrapped resampled pre-DBS entropy estimates


ppar.neuTypeFilter = 'SU';
ppar.hzThresh = SPARSE_HZ;
ppar.predbsTime = PREDBS_TIME;
ppar.dbsTime = DBS_TIME;
ppar.ordH = ORD_H;
ppar.binsPD = BINS_PER_DECADE;
ppar.nBoot = NBOOTS;
ppar.subjID = 'Analysis';

[codeDirectoryFullPath, codeName] = fileparts(mfilename('fullpath'));

%% 

% load the master table to work on each row
load([ppar.tablePath, '\', 'tableRoot'], 'tableRoot');

% merge the table of rate analysis so that we can exclude any rows with
% too-sparse data (or no spikes in PRE and/or DBS condition)
% load RateBins_60sec_1secBins table
load([ppar.tablePath, '\', 'RateBins_60sec_1secBins'], 'RateBins');

% rename matFile and matFileFolder variables:
varNames = RateBins.Properties.VariableNames;
varNames{strcmp(varNames, 'matFile')} = 'matFile_RateBins';
varNames{strcmp(varNames, 'matFileFolder')} = 'matFileFolder_RateBins';
RateBins.Properties.VariableNames = varNames;

N = join(tableRoot, RateBins, 'LeftKeys', 'objectID', 'RightKeys', 'objectID');


%%  Append two new columns to RateBins table of meanRates (PRE and DBS) based
% on 30-second portions only, as described at the top of this script. 

% Initialize new columns:
nRows = size(N, 1);
meanRatePRE_30s = zeros(nRows, 1);
meanRateDBS_30s = zeros(nRows, 1);


% Fill new columns with mean rate values: load each row's RateBins matfile
% to get the bin values ("bins" struct)
disp('calculating spkRate means...')
tic
for iRow = 1:nRows
    fullFilePath = [ppar.projRootPath, '\', N.matFileFolder_RateBins{iRow}, ...
            '\', N.matFile_RateBins{iRow}, '.mat'];
    load(fullFilePath);

    % get rate-bins from last 30 seconds worth of data for pre and dbs
    % conditions: (in case the collection of rate bins is less than 30 
    % seconds-worth, output all of them)
    ratesPRE_iRow = getElementsSubselect(bins.pre, 30, 'direction', 'RtoL');
    meanRatePRE_30s(iRow) = mean(ratesPRE_iRow);
    ratesDBS_iRow = getElementsSubselect(bins.dbs, 30, 'direction', 'RtoL');
    meanRateDBS_30s(iRow) = mean(ratesDBS_iRow);
    
end
toc
disp('means calculated')

% Append new columns to table
N = [N, table(meanRatePRE_30s), table(meanRateDBS_30s)];



%% 

% filter out any rows that have too-sparse cells, based on 30-second
% portions
Nnew = N;
isPREbelow = Nnew{:, 'meanRatePRE_30s'} < SPARSE_HZ;
isDBSbelow = Nnew{:, 'meanRateDBS_30s'} < SPARSE_HZ;
isRemove = isPREbelow | isDBSbelow;
Nnew(isRemove, :) = [];


%% Perform Entropy (Direct-ISI) estimation on each file's spike train data
% 

nRows = size(Nnew, 1);

% Initialize the variables to be filled for results table:
     objectID = cell(nRows, 1);
      matFile = cell(nRows, 1);
matFileFolder = cell(nRows, 1);
     H_PREemp = zeros(nRows, 1);
     H_DBSemp = zeros(nRows, 1);
  H_PREbootAv = zeros(nRows, 1);
H_PREbootStdv = zeros(nRows, 1);
         pVal = zeros(nRows, 1);

% Calc entropy values for each file:
disp('Calculating Direct Entropy estimates...')
tic
for iNex = 1:nRows
    % load i row of nexFile
    % iNex = 100;
    objectID{iNex} = Nnew.objectID{iNex};
    nexfn = Nnew.nexFile{iNex,1};
    nexpn = Nnew.nexFileFolder{iNex,1};
    nexFile = readNexFile([ppar.projRootPath, '\', nexpn, '\', nexfn]);


    % build entropyDirISI struct
    [spkTimes, StimTs] = parseNexFile(nexFile);

    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = StimTs.DBS;
    stimPeriod = median(diff(dbsTimes));


    % get pre-DBS spikes
    isPreDBS = (spkTimes < dbsTimes(1)) & ...
               (spkTimes >= (dbsTimes(1) - PREDBS_TIME));
    spksPRE = spkTimes(isPreDBS);


    % get DBS-on spikes
    isDBSon = (spkTimes >= (dbsTimes(1) + 30)) & ...
              (spkTimes < (dbsTimes(1) + 30 + DBS_TIME));
    spksDBS = spkTimes(isDBSon);

    
    
    %% CALCULATE DBS-ON Entropies for current Nexfile data
    
    spksDBS = sort(spksDBS); % make sure the seconds are monotonically increasing...
    % Define ISIs, and remove any ISIs == 0
    isiDBS = diff(spksDBS);
    isiDBS(isiDBS == 0) = [];
    
    
    iH_DBS = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);
    entropyDirISI.dbsEmp = iH_DBS;
    H_DBSemp(iNex) = iH_DBS;
    
    
    
    %% CALCULATE PRE-DBS Entropies for current Nexfile data
    % Note that this estimate is not scaled to have the same number of
    % isi's as the DBS case.
    
    spksPRE = sort(spksPRE); % make sure the seconds are monotonically increasing...
    % Define ISIs, and remove any ISIs == 0
    isiPRE = diff(spksPRE);
    isiPRE(isiPRE == 0) = [];
    
    iH_PREemp = entropyISIdirect_multOrder(isiPRE, BINS_PER_DECADE, ORD_H);
    entropyDirISI.preEmp = iH_PREemp;
    H_PREemp(iNex) = iH_PREemp;
    
    
    
    %% CALCULATE Bootstrapped distribution of Pre-DBS Entropies for current Nexfile data
    % Note that these bootstrap estimates all have the same number of isi's
    % as the DBS case.
    
% tic
%     nBoots = pipeParams.nBoot;
    iH_PREboot = zeros(NBOOTS, ORD_H);
    parfor iBoot = 1:NBOOTS
        % Resample PRE ISIs (with replacement) so that number of ISIs
        % matches the DBS case
        [isiPREresamp] = datasample(isiPRE, numel(isiDBS));
        
%         % Get Log-binned version of ISI
%         [isiLogHist, logBinEdges, binIdx] = isiLogBinned(isiPREresamp, ...
%                                                          BINS_PER_DECADE);
        
        % Calculate multi-order Direct-Entropy
        iH_PREboot(iBoot,:) = entropyISIdirect_multOrder(isiPREresamp, BINS_PER_DECADE, ORD_H);
                                                     
    end
% toc   
    entropyDirISI.preBoot = iH_PREboot;
    H_PREbootAv(iNex) = mean(iH_PREboot);
    H_PREbootStdv(iNex) = std(iH_PREboot);


    
    %% Calculate the p-value difference between First-order Entropies
    
%     iH_PREbootAv = mean(iH_PREboot(:,1));
        
    if iH_DBS(1) <= H_PREbootAv(iNex) % DBS-entropy is lower than PRE
        isLess = iH_PREboot(:,1) < iH_DBS(1);
        iPval = sum(isLess) / numel(isLess); % p-value
        
    else % if DBS-entropy is higher than PRE
        isMore = iH_PREboot(:,1) > iH_DBS(1);
        iPval = sum(isMore) / numel(isMore);
        
    end
    
    pVal(iNex) = iPval;


    %% Save it with a proper name and in right location
    formatSpec = '%s_%s_%dhzThresh_%dpre_%ddbs_ordH%d_%dbinsPD_%dboots.mat';
    matfnStr = sprintf(formatSpec, Nnew.objectID{iNex}, ...
                              ppar.neuTypeFilter, ...
                              ppar.hzThresh, ...
                              ppar.predbsTime, ...
                              ppar.dbsTime, ...
                              ppar.ordH, ...
                              ppar.binsPD, ...
                              ppar.nBoot);

    fullPathFn = [ppar.projRootPath, '\', MATFILESAVEPATH, '\', matfnStr];
    matFile{iNex} = matfnStr;
    matFileFolder{iNex} = MATFILESAVEPATH;

    save(fullPathFn, 'entropyDirISI');

% update a row-values in entropy table with the relevant information 


% repeat
end
toc
disp('Done Calculating!')

% halelujah!
load handel;
player = audioplayer(y, Fs);
play(player);



%% Create RESULTS table and SAVE

t = now;
d = datetime(t, 'ConvertFrom', 'datenum');

% combine columns into final table:
E = [table(objectID), ...
     table(matFile), ...
     table(matFileFolder), ...
     table(H_PREemp), ...
     table(H_DBSemp), ...
     table(H_PREbootAv), ...
     table(H_PREbootStdv), ...
     table(pVal)];
                         
E.Properties.UserData.generatingCode.codeName = [codeName, '.m'];
E.Properties.UserData.generatingCode.codeDirectoryFullPath = codeDirectoryFullPath;
E.Properties.UserData.generatingCode.dateExecuted = d;
EntropyDirectISI_analysis = E;     


% build the name of the table correctly:              
formatSpec = 'EntropyDirectISI_%s_%s_%dhzThresh_%dpre_%ddbs_ordH%d_%dbinsPD_%dboots.mat';
    tableStr = sprintf(formatSpec, ppar.subjID, ...
                              ppar.neuTypeFilter, ...
                              ppar.hzThresh, ...
                              ppar.predbsTime, ...
                              ppar.dbsTime, ...
                              ppar.ordH, ...
                              ppar.binsPD, ...
                              ppar.nBoot);     
                          
fullTablePath = [ppar.tablePath, '\', tableStr];
save(fullTablePath, 'EntropyDirectISI_analysis');
                             
disp('SUCCESS!')

%%

% % in seconds. Time to include spikes for analysis before DBS onset
% PREDBS_TIME = pipeParams.predbsTime; 
% DBS_TIME = pipeParams.dbsTime;  % also the time that DBS was on for the experiment 
% 
% % PARAMETERS for Direct Entropy estimate with log binned ISI:
% ORD_H = pipeParams.ordH; % how many orders of entropy to use for linear extrapolation
% BINS_PER_DECADE = pipeParams.binsPD; % of log-spaced bins for ISI histogram
% 
% 
% 
% %%  MAIN FOR-LOOP
% 
% nNex = size(N, 1);
% 
% % Initialize the variables to be filled for each Nexfile
%     H_DBS = cell(nNex, 1);
%  H_PREemp = cell(nNex, 1);
% H_PREboot = cell(nNex, 1);
%      pVal = zeros(nNex, 1);
% 
% % H_delta = zeros(nNex, 1);
% 
% disp('Calculating Direct Entropy estimates...')
% tic
% for iNex = 1:nNex
%     %% LOAD each Nexfile and get spike times
%     
%     % Get each NEX file from table
%     nexfn = N.Filename{iNex,1};
%     nexpn = N.Pathname{iNex,1};
% 
%     nexFile = readNexFile([nexpn, '\', nexfn]);
%     [spkTimes, StimTs] = parseNexFile(nexFile);
% 
%     % separate the spike times into pre-DBS and DBS-on
%     dbsTimes = StimTs.DBS;
%     stimPeriod = median(diff(dbsTimes));
% 
% 
%     % get pre-DBS spikes
%     isPreDBS = (spkTimes < dbsTimes(1)) & ...
%                (spkTimes >= (dbsTimes(1) - PREDBS_TIME));
%     spksPRE = spkTimes(isPreDBS);
% 
% 
%     % get DBS-on spikes
%     isDBSon = (spkTimes >= dbsTimes(1)) & ...
%               (spkTimes < (dbsTimes(end) + stimPeriod));
%     spksDBS = spkTimes(isDBSon);
% 
%     
%     
%     %% CALCULATE DBS-ON Entropies for current Nexfile data
%     
%     % Define ISIs, and remove any ISIs == 0
%     isiDBS = diff(spksDBS);
%     isiDBS(isiDBS == 0) = [];
%     
%     iH_DBS = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);
%     
%     
%     
%     %% CALCULATE PRE-DBS Entropies for current Nexfile data
%     % Note that this estimate is not scaled to have the same number of
%     % isi's as the DBS case.
%     
%     % Define ISIs, and remove any ISIs == 0
%     isiPRE = diff(spksPRE);
%     isiPRE(isiPRE == 0) = [];
%     
%     iH_PREemp = entropyISIdirect_multOrder(isiPRE, BINS_PER_DECADE, ORD_H);
%     
%     
%     
%     %% CALCULATE Bootstrapped distribution of Pre-DBS Entropies for current Nexfile data
%     % Note that these bootstrap estimates all have the same number of isi's
%     % as the DBS case.
%     
% 
%     nBoots = pipeParams.nBoot;
%     iH_PREboot = zeros(nBoots, ORD_H);
%     parfor iBoot = 1:nBoots
%         % Resample PRE ISIs (with replacement) so that number of ISIs
%         % matches the DBS case
%         [isiPREresamp] = datasample(isiPRE, numel(isiDBS));
%         
% %         % Get Log-binned version of ISI
% %         [isiLogHist, logBinEdges, binIdx] = isiLogBinned(isiPREresamp, ...
% %                                                          BINS_PER_DECADE);
%         
%         % Calculate multi-order Direct-Entropy
%         iH_PREboot(iBoot,:) = entropyISIdirect_multOrder(isiPREresamp, BINS_PER_DECADE, ORD_H);
%                                                      
%     end
%         
% 
%     
%     %% Calculate the p-value difference between First-order Entropies
%     
%     iH_PREbootAv = mean(iH_PREboot(:,1));
%         
%     if iH_DBS(1) <= iH_PREbootAv % DBS-entropy is lower than PRE
%         isLess = iH_PREboot(:,1) < iH_DBS(1);
%         iPval = sum(isLess) / numel(isLess); % p-value
%         
%     else % if DBS-entropy is higher than PRE
%         isMore = iH_PREboot(:,1) > iH_DBS(1);
%         iPval = sum(isMore) / numel(isMore);
%         
%     end
%     
%     
% 
%     %% UPDATE Entropy result variables
%     
%     % Initialize the variables to be filled 
%         H_DBS{iNex,1} = iH_DBS;
%      H_PREemp{iNex,1} = iH_PREemp;
%     H_PREboot{iNex,1} = iH_PREboot;
%          pVal(iNex,1) = iPval;
% 
%     
%     
% end % END FOR-LOOP
% toc
% disp('Done Calculating!')
% 
% % halelujah!
% load handel;
% player = audioplayer(y, Fs);
% play(player);
% 
% 
% 
% %% Create RESULTS table and SAVE
% 
% PE = table(H_PREemp);
% PB = table(H_PREboot);
% 
% D = table(H_DBS);
% PVAL = table(pVal);
% 
% Hdir_ISI_Results = [N, PE, PB, D, PVAL];
% 
% disp('SUCCESS!')
% 
% 
% 
% 
% 
% 
% % end
% 
% 
% 
% 
% 






