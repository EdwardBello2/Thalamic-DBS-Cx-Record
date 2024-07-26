% Function for calculating Direct Entropy estimate of spike activity,
% bootstrapping the pre-dbs condition to see if the dbs condition is
% significantly different
%
% INPUTS:
% N - matlab table of values to be analyzed
% pipeParams - struct with user-spec analysis parameters

% TO-DO:
%

function [Hdir_ISI_Results] = calcDirectEntropyISIbootstrap_batch(N, pipeParams)


% in seconds. Time to include spikes for analysis before DBS onset
PREDBS_TIME = pipeParams.predbsTime; 
DBS_TIME = pipeParams.dbsTime;  % also the time that DBS was on for the experiment 

% PARAMETERS for Direct Entropy estimate with log binned ISI:
ORD_H = pipeParams.ordH; % how many orders of entropy to use for linear extrapolation
BINS_PER_DECADE = pipeParams.binsPD; % of log-spaced bins for ISI histogram



%%  MAIN FOR-LOOP

nNex = size(N, 1);

% Initialize the variables to be filled for each Nexfile
    H_DBS = cell(nNex, 1);
 H_PREemp = cell(nNex, 1);
H_PREboot = cell(nNex, 1);
     pVal = zeros(nNex, 1);

% H_delta = zeros(nNex, 1);

disp('Calculating Direct Entropy estimates...')
tic
for iNex = 1:nNex
    %% LOAD each Nexfile and get spike times
    
    % Get each NEX file from table
    nexfn = N.Filename{iNex,1};
    nexpn = N.Pathname{iNex,1};

    nexFile = readNexFile([nexpn, '\', nexfn]);
    [spkTimes, StimTs] = parseNexFile(nexFile);

    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = StimTs.DBS;
    stimPeriod = median(diff(dbsTimes));


    % get pre-DBS spikes
    isPreDBS = (spkTimes < dbsTimes(1)) & ...
               (spkTimes >= (dbsTimes(1) - PREDBS_TIME));
    spksPRE = spkTimes(isPreDBS);


    % get DBS-on spikes
    isDBSon = (spkTimes >= dbsTimes(1)) & ...
              (spkTimes < (dbsTimes(end) + stimPeriod));
    spksDBS = spkTimes(isDBSon);

    
    
    %% CALCULATE DBS-ON Entropies for current Nexfile data
    
    % Define ISIs, and remove any ISIs == 0
    isiDBS = diff(spksDBS);
    isiDBS(isiDBS == 0) = [];
    
    iH_DBS = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);
    
    
    
    %% CALCULATE PRE-DBS Entropies for current Nexfile data
    % Note that this estimate is not scaled to have the same number of
    % isi's as the DBS case.
    
    % Define ISIs, and remove any ISIs == 0
    isiPRE = diff(spksPRE);
    isiPRE(isiPRE == 0) = [];
    
    iH_PREemp = entropyISIdirect_multOrder(isiPRE, BINS_PER_DECADE, ORD_H);
    
    
    
    %% CALCULATE Bootstrapped distribution of Pre-DBS Entropies for current Nexfile data
    % Note that these bootstrap estimates all have the same number of isi's
    % as the DBS case.
    

    nBoots = pipeParams.nBoot;
    iH_PREboot = zeros(nBoots, ORD_H);
    parfor iBoot = 1:nBoots
        % Resample PRE ISIs (with replacement) so that number of ISIs
        % matches the DBS case
        [isiPREresamp] = datasample(isiPRE, numel(isiDBS));
        
%         % Get Log-binned version of ISI
%         [isiLogHist, logBinEdges, binIdx] = isiLogBinned(isiPREresamp, ...
%                                                          BINS_PER_DECADE);
        
        % Calculate multi-order Direct-Entropy
        iH_PREboot(iBoot,:) = entropyISIdirect_multOrder(isiPREresamp, BINS_PER_DECADE, ORD_H);
                                                     
    end
        

    
    %% Calculate the p-value difference between First-order Entropies
    
    iH_PREbootAv = mean(iH_PREboot(:,1));
        
    if iH_DBS(1) <= iH_PREbootAv % DBS-entropy is lower than PRE
        isLess = iH_PREboot(:,1) < iH_DBS(1);
        iPval = sum(isLess) / numel(isLess); % p-value
        
    else % if DBS-entropy is higher than PRE
        isMore = iH_PREboot(:,1) > iH_DBS(1);
        iPval = sum(isMore) / numel(isMore);
        
    end
    
    

    %% UPDATE Entropy result variables
    
    % Initialize the variables to be filled 
        H_DBS{iNex,1} = iH_DBS;
     H_PREemp{iNex,1} = iH_PREemp;
    H_PREboot{iNex,1} = iH_PREboot;
         pVal(iNex,1) = iPval;

    
    
end % END FOR-LOOP
toc
disp('Done Calculating!')

% halelujah!
load handel;
player = audioplayer(y, Fs);
play(player);



%% Create RESULTS table and SAVE

PE = table(H_PREemp);
PB = table(H_PREboot);

D = table(H_DBS);
PVAL = table(pVal);

Hdir_ISI_Results = [N, PE, PB, D, PVAL];

disp('SUCCESS!')






end











