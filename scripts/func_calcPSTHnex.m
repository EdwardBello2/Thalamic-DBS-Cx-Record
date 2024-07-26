function psthnex = func_calcPSTHnex(nexFile, objectID, SPKINTERVselect, psthType, ppar)

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

formatSpecFn = 'h%s_%ss-%ss_fromDbsOnset_%smsBW_%sBins_%s';
matfnStr = sprintf(formatSpecFn, objectID, ...
                          num2str(SPKINTERVselect(1)), ...
                          num2str(SPKINTERVselect(2)), ...
                          num2str(ppar.psthBinWidth*1000), ...
                          num2str(ppar.psthNumBins), ...
                          trimStr);
fullPathFn = [scriptIntermediateDataFolder '\' matfnStr '.mat'];



%% Get the PSTH

% load/calculate entropy for this nexfile
if exist(fullPathFn, 'file')
    load(fullPathFn, 'psthnex')

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
    [spkTimesDbs] = getIntervalEvents(spkTimes, refT, spkIntervalSelect);
    spkTimesDbs = sort(spkTimesDbs); % make sure the seconds are monotonically increasing...

    % get DBS times in relevent interval
    dbsTimes = getIntervalEvents(StimTs.DBS, refT, spkIntervalSelect);

    % Calculate empirical psth for DBS
    if isempty(spkTimesDbs)
        psth_DBSct = zeros(1, length(binEdgesPsth)-1);

    else
        psth_DBSct = psth(spkTimesDbs, dbsTimes, binEdgesPsth, 'count');

    end
    % If user specified to remove first bin, perform that now:
    if ppar.trimPSTH % if trimPSTH is '[]', this step is skipped
        psth_DBSct(ppar.trimPSTH) = 0;

    end

    switch psthType
        case 'count'
            psthnex = psth_DBSct;

        case 'fr'
            psthnex = psth_DBSct / (numel(dbsTimes) * bw); % norm by spikes per second

        case 'prob'
            psthnex = psth_DBSct / sum(psth_DBSct);

        otherwise
            error('Incorrect string for psthType')

    end     

    
    save(fullPathFn, 'psthnex');
    
end




end
