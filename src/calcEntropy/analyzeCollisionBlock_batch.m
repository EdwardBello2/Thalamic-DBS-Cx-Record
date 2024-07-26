% Function for Collision-block analysis for each row of the input table N
% 
%
% INPUTS:
% N - matlab table of rows to be analyzed
% pipeParams - struct with user-spec analysis parameters

% TO-DO:
%

function CollisionBlockResults = analyzeCollisionBlock_batch(N, pipeParams); %<---------


% in seconds. Time to include spikes for analysis before DBS onset
PREDBS_TIME = pipeParams.predbsTime; 
DBS_TIME = pipeParams.dbsTime;  % also the time that DBS was on for the experiment 


bBeg = pipeParams.binEdgeLeft;
bw = pipeParams.binWidth;
nBins = pipeParams.binNum;
binEdges = bBeg:bw:(bw * nBins);
% % PARAMETERS for Direct Entropy estimate with log binned ISI:
% ORD_H = pipeParams.ordH; % how many orders of entropy to use for linear extrapolation
% BINS_PER_DECADE = pipeParams.binsPD; % of log-spaced bins for ISI histogram
% 


%%  MAIN FOR-LOOP

nNex = size(N,1);

Hpre = zeros(nNex,1);
Hdbs = zeros(nNex,1);
% H_delta = zeros(nNex,1);

NexpMin_av = 0.0358;
minEst = NexpMin_av * nNex;
disp(['Calculating PSTH-Letter Entropy... (est. ', num2str(minEst), ' mins)'])


% Initialize variables to fill for each Nexfile
H_PRE = zeros(nNex,1);
H_DBS = zeros(nNex,1);

H_PREbootdistr = cell(nNex,1);
H_DBSbootdistr = cell(nNex,1);

pVal = zeros(nNex,1);

% Columns to update for each Nex file
           nStims = zeros(nNex, 1);% # of DBS pulses in trial
         nPhsLock = zeros(nNex, 1);% # of pulses that had a phase-lock spike

         nTcMinus = zeros(nNex, 1);
nMissDueToTcMinus = zeros(nNex, 1);
      nMissSimple = zeros(nNex, 1);

tic
for iNex = 1:nNex
% iNex
    % Get each NEX file from table
    nexfn = N.Filename{iNex,1};
    nexpn = N.Pathname{iNex,1};

    nexFile = readNexFile([nexpn, '\', nexfn]);
    [spkTimes, stims] = parseNexFile(nexFile);

    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = stims.DBS;
%     stimPeriod = median(diff(dbsTimes));


%     % get pre-DBS spikes
%     isPreDBS = (spkTimes < dbsTimes(1)) & ...
%                (spkTimes >= (dbsTimes(1) - PREDBS_TIME));
%     spksPRE = spkTimes(isPreDBS);


    % get DBS-on spikes
    isDBSon = (spkTimes >= dbsTimes(1) - binEdges(end)) & ...
              (spkTimes < (dbsTimes(end) + binEdges(end)));
    spksDBS = spkTimes(isDBSon);



    %% Count phase-locked phenomena
    
    nPulses = numel(dbsTimes);
    % number of standard deviations to include in the borders of the 
    % phase-lock bin
    nStd = 3;
    
    % Define antidromic conduction time tc based on mean spike latency
    tc = N.latency_mean(iNex); % time for anitdromic conduction
    
%     % Create "phase-locked" spike bin, based on mean and stdv spike latency
%     phsLockBin = [(tc - (nStd * N.latency_std(iNex))), ...
%                   (tc + (nStd * N.latency_std(iNex)))];
%     if phsLockBin(1) < binEdges(1), phsLockBin(1) = binEdges(1); end
%     if phsLockBin(2) > binEdges(end), phsLockBin(2) = binEdges(end); end

    % Create "phase-locked" spike bin, based on pipeParams.latencyWin
    phsLockBin = pipeParams.latencyWin;
    
    

    % Count # of "hits" in the phase-locked bin
    spkLatencies = findEvokedSpikeLatency(spksDBS, dbsTimes, [phsLockBin(1), phsLockBin(2)]);
    isMiss = isnan(spkLatencies);
    isHit = ~isMiss;
        
    
    % Count # of misses due to tcMinus interruption
    missIdx = find(isMiss);
    nMiss = numel(missIdx);
    
    % see if any pulse has a spike preceding it within tc-
    isInTcMinus = false(nPulses, 1); 
    for iPulse = 1:nPulses
        dSpks = spksDBS - dbsTimes(iPulse);
        dSpksPre = dSpks(dSpks < 0);
        
        if max(dSpksPre) > -tc
            isInTcMinus(iPulse,1) = true;
            
        end
        
    end
    
    isTcInterrMiss = (isMiss & isInTcMinus);
    isSimpleMiss = (isMiss & ~isInTcMinus); 
    
    
    % Columns to update for each Nex file
               nStims(iNex,1) = nPulses;
             nPhsLock(iNex,1) = sum(isHit);
             nTcMinus(iNex,1) = sum(isInTcMinus);
    nMissDueToTcMinus(iNex,1) = sum(isTcInterrMiss);
          nMissSimple(iNex,1) = sum(isSimpleMiss);

    
    allMiss = sum(isTcInterrMiss) + sum(isSimpleMiss);
    if (sum(isHit) + allMiss) ~= nPulses, error('Unaccounted pulses'); end
        
end % END iNex for loop
    
    
   
toc
disp('Done Calculating!')

% % halelujah!
% load handel;
% player = audioplayer(y, Fs);
% play(player);



%% Create RESULTS table and SAVE

NSTIM = table(nStims);
PL = table(nPhsLock);
TCM = table(nTcMinus);
TCMiss = table(nMissDueToTcMinus);
SM = table(nMissSimple);


CollisionBlockResults = [N, NSTIM, PL, TCM, TCMiss, SM];






disp('SUCCESS!')






end

% % Three reasons why a phase-locked spike might have missed:
%     % 1) a spike occurred just before the DBSpulse in the tc- period,
%     % 2) a spike occurred just after the DBSpulse in the tc+ period,
%     % 3) unrelated reason, perhaps "simply" just misses at random...
%     isInterrTcMinus = false(nPulses, 1);
%     isInterrTcPlus = false(nPulses, 1);
%     isSimpleMiss = false(nPulses, 1);
%     
%     % For each stim that had a phase-lock MISS, find out why
%     for iMiss = 1:nMiss
%          pulseTime = dbsTimes(missIdx(iMiss));
%          dSpks = spkTimes - pulseTime;
%          dSpksPre = dSpks(dSpks < 0);
%          dSpksPos = dSpks(dSpks >= 0);
%          
%          % is the closest spike before this pulse within tc-?
%          if max(dSpksPre) > -tc
%              isInterrTcMinus(missIdx(iMiss),1) = true;
%              
%          elseif min(dSpksPos) < tc
%              isInterrTcPlus(missIdx(iMiss),1) = true;
%              
%          else
%              isSimpleMiss(missIdx(iMiss),1) = true;
%          
%          end
%     end
% 
%     
%     % Columns to update for each Nex file
%                    nStims(iNex,1) = nPulses;
%                  nPhsLock(iNex,1) = sum(isHit);
%     nPhsLockInterrTcMinus(iNex,1) = sum(isInterrTcMinus);
%      nPhsLockInterrTcPlus(iNex,1) = sum(isInterrTcPlus);
%        nPhsLockSimpleMiss(iNex,1) = sum(isSimpleMiss);














