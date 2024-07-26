function HpsthBoots = func_calcHPSTHnex_bootstrap(nexFile, objectID, SPKINTERVselect, NBOOTS, ppar)
% output firing rate calculations; if already cached, simply load the data, or else perform the calculation... 

%% Assert that intermediate data cache folder exists
[~, scriptName] = fileparts(mfilename('fullpath'));

% First make sure that this script has a folder within the project folder's
% intermediate data section
scriptIntermediateDataFolder = [ppar.projRootPath, '\', ppar.intDataRootFolder, '\', scriptName];
if ~exist(scriptIntermediateDataFolder, 'dir')
    mkdir(scriptIntermediateDataFolder)

end


%% Set up data cache name

if ppar.trimPSTH
    trimStr = 'Trim';

else
    trimStr = 'noTrim';

end

formatSpecFn = 'h%s_%ss-%ss_fromDbsOnset_%smsBW_%sBins_%sboots_%s';
matfnStr = sprintf(formatSpecFn, objectID, ...
                          num2str(SPKINTERVselect(1)), ...
                          num2str(SPKINTERVselect(2)), ...
                          num2str(ppar.psthBinWidth*1000), ...
                          num2str(ppar.psthNumBins), ...
                          num2str(NBOOTS), ...
                          trimStr);
fullPathFn = [scriptIntermediateDataFolder '\' matfnStr '.mat'];



%% Calculate entropy, or load it if already cached previously 

% load/calculate entropy for this nexfile
if exist(fullPathFn, 'file')
    load(fullPathFn, 'HpsthBoots')

else
    % extract spike times and DBS stim times from nexFile struct
    [spkTimes, StimTs] = parseNexFile(nexFile);

    % make sure spkTimes are all of time values that are monotonically
    % increasing:
    spkTimes = sort(spkTimes);


    % Set binEdges based on input parameters:
    bBeg = ppar.psthTimeBeg;
    bw = ppar.psthBinWidth;
    nBins = ppar.psthNumBins;
    binEdgesPsth = bBeg:bw:(bw * nBins);


    % If SPKINTERV specifies times that fall outside of the recording,
    % adjust SPKINTERV to have boundaries within recorded times
    spkIntervalSelect = SPKINTERVselect; 
    tbegRef = nexFile.tbeg - StimTs.DBS(1);
    tendRef = nexFile.tend - StimTs.DBS(1);
    spkIntervalSelect(1) = max(spkIntervalSelect(1), tbegRef);
    spkIntervalSelect(2) = min(spkIntervalSelect(2), tendRef);


    % Get spike times occurring within DBS interval 30-60
    refT = StimTs.DBS(1);

    % get spike times in relevent interval
    [spkTimesInterv] = getIntervalEvents(spkTimes, refT, spkIntervalSelect);
    spkTimesInterv = sort(spkTimesInterv); % make sure the seconds are monotonically increasing...

    % get DBS times in relevent interval
    dbsTimes = getIntervalEvents(StimTs.DBS, refT, spkIntervalSelect);

    % Calculate empirical psth-Entropy for DBS
    if isempty(spkTimesInterv)
        psth_Intervct = zeros(1, length(binEdgesPsth)-1);

    else
        psth_Intervct = psth(spkTimesInterv, dbsTimes, binEdgesPsth, 'count');

    end
    % If user specified to remove first bin, perform that now:
    if ppar.trimPSTH % if trimPSTH is '[]', this step is skipped
        psth_Intervct(ppar.trimPSTH) = 0;

    end
    nSpksInterv = sum(psth_Intervct);


    % Calculate bootstrapped distribution of psth-Entropies for DBS
    HpsthBoots = zeros(NBOOTS, 1);
    unifRange = binEdgesPsth;
    unifRange(ppar.trimPSTH) = [];
    A = unifRange(1);
    B = unifRange(end);
    if (nSpksInterv < 1) || isempty(nSpksInterv)
        HpsthBoots(:) = NaN;

    else
        for iBoot = 1:NBOOTS
            % for each bootstrap PSTH: based on current "allowed" bins (trimming?), 
            % random-uniformly sample spike counts for the bins
            R = unifrnd(A, B, nSpksInterv, 1);
            [psth_iBOOTct,~] = histcounts(R, binEdgesPsth);

            psth_iBOOTprob = psth_iBOOTct / nSpksInterv;

            % gather bootstrapped entropy distribution
            HpsthBoots(iBoot) = entropyLetter_bitpSpike(psth_iBOOTprob);


        end

    end


    % Save results in cache for later retreival, safe a crap ton of
    % time!
    save(fullPathFn, 'HpsthBoots');

end

   
    
end