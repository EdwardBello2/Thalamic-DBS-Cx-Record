% Function for calculating PSTH-based "Letter" Entropy estimate of
% stim-locked neural activity
%
% INPUTS:
% N - matlab table of values to be analyzed
% pipeParams - struct with user-spec analysis parameters

% TO-DO:
%

function [EntropyResults] = calcEntropyPsthLetterVarBins_batch(N, pipeParams)


% in seconds. Time to include spikes for analysis before DBS onset
PREDBS_TIME = pipeParams.predbsTime; 
DBS_TIME = pipeParams.dbsTime;  % also the time that DBS was on for the experiment 





%%  MAIN FOR-LOOP

nNex = size(N,1);

Hpre = zeros(nNex,1);
Hdbs = zeros(nNex,1);
% H_delta = zeros(nNex,1);

% NexpMin_av = 0.0358;
% minEst = NexpMin_av * nNex;
% disp(['Calculating PSTH-Letter Entropy... (est. ', num2str(minEst), ' mins)'])


% Initialize variables to fill for each Nexfile
H_PRE = zeros(nNex,1);
H_DBS = zeros(nNex,1);

H_PREbootdistr = cell(nNex,1);
H_DBSbootdistr = cell(nNex,1);

pVal = zeros(nNex,1);


tic
for iNex = 1:nNex

    %% Calculate binEdges depending on DBS period length
    
    bBeg = pipeParams.binEdgeLeft;
    bw = pipeParams.binWidth;
    DBSint = 1 / N.dbsFrequency(iNex,1);
    nBins = floor(DBSint / bw);
    binEdges = bBeg:bw:(bw * nBins);
    
    
    
    %% Get spike times for both PRE- and DBS- times
    
    % Get each NEX file from table
    nexfn = N.Filename{iNex,1};
    nexpn = N.Pathname{iNex,1};

    nexFile = readNexFile([nexpn, '\', nexfn]);
    [spkTimes, stims] = parseNexFile(nexFile);

    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = stims.DBS;
    stimPeriod = median(diff(dbsTimes));


    % get pre-DBS spikes
    isPreDBS = (spkTimes < dbsTimes(1)) & ...
               (spkTimes >= (dbsTimes(1) - PREDBS_TIME));
    spksPRE = spkTimes(isPreDBS);


    % get DBS-on spikes
    isDBSon = (spkTimes >= dbsTimes(1)) & ...
              (spkTimes < (dbsTimes(end) + stimPeriod));
    spksDBS = spkTimes(isDBSon);



    %% Calculate Entropies

    nBoot = pipeParams.nBoot;
    
    
    % Calculate empirical Entropy for PRE and DBS
    psthPRE = psth(spksPRE, stims.VirtPre, binEdges);
    % If user specified to remove first bin, perform that now:
    if pipeParams.trimPSTH
        psthPRE(1) = [];
        
    end
    
    Hpre = entropyLetter_bitpSpike(psthPRE);
    
    psthDBSct = psth(spksDBS, stims.DBS, binEdges, 'count');
    % If user specified to remove first and last bins, perform that now:
    if pipeParams.trimPSTH
        psthDBSct(1) = [];
        
    end
    
    nSpksDBS = sum(psthDBSct);
    
    psthDBS = psthDBSct / nSpksDBS;
    Hdbs = entropyLetter_bitpSpike(psthDBS);
    
    
    % Calculate bootstrapped-resampled Entropy for PRE
    psthRows = psthBootstrap(spksPRE, stims.VirtPre, binEdges, nBoot, nSpksDBS);
    % If user specified to remove first and last bins, perform that now:
    if pipeParams.trimPSTH
        psthRows(:,1) = [];
        
    end
    nRow = size(psthRows, 1);
    Hprebootdistr = zeros(nRow, 1);
    for iRow = 1:nRow
        Hprebootdistr(iRow) = entropyLetter_bitpSpike(psthRows(iRow,:));

    end
    
    
    % Calculate bootstrapped-resampled Entropy for DBS
    psthRows = psthBootstrap(spksDBS, stims.DBS, binEdges, nBoot, nSpksDBS);
    if pipeParams.trimPSTH
        psthRows(:,1) = [];
        
    end
    nRow = size(psthRows, 1);
    Hdbsbootdistr = zeros(nRow, 1);
    for iRow = 1:nRow
        Hdbsbootdistr(iRow) = entropyLetter_bitpSpike(psthRows(iRow,:));

    end
    
      
    % Calculate the p-value. Usually the modulated cells have lower H in 
    % the on-stim period compared to the pr-stim period, so look for the 
    % left tail of the p-value.
    isLess = Hprebootdistr < Hdbs;
    p = sum(isLess) / nBoot; % p-value
    
    
    % Store results to be made into table later
    H_PRE(iNex,1) = Hpre;
    H_DBS(iNex,1) = Hdbs;
    
    H_PREbootdistr{iNex,1} = Hprebootdistr;
    H_DBSbootdistr{iNex,1} = Hdbsbootdistr;
    
    pVal(iNex,1) = p;
    
    
end % END for loop
    
    
   
toc
disp('Done Calculating!')

% halelujah!
load handel;
player = audioplayer(y, Fs);
play(player);



%% Create RESULTS table and SAVE

HPRE = table(H_PRE);
HDBS = table(H_DBS);

HPREb = table(H_PREbootdistr);
HDBSb = table(H_DBSbootdistr);

P = table(pVal);


EntropyResults = [N, HPRE, HDBS, HPREb, HDBSb, P];






disp('SUCCESS!')






end












