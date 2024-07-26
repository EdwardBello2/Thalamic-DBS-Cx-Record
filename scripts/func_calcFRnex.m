function FR = func_calcFRnex(nexFile, objectID, SPKINTERV, ppar)
% output firing rate calculations; if already cached, simply load the data, or else perform the calculation... 

    % Get the name of the currently running script:
    [scriptDirectoryFullPath, scriptName] = fileparts(mfilename('fullpath'));

    % First make sure that this script has a folder within the project folder's
    % intermediate data section
    scriptIntermediateDataFolder = [ppar.projRootPath, '\', ppar.intDataRootFolder, '\', scriptName];
    if ~exist(scriptIntermediateDataFolder, 'dir')
        mkdir(scriptIntermediateDataFolder)

    end


    % extract spike times and DBS stim times from nexFile struct
    [spkTimes, StimTs] = parseNexFile(nexFile);

    % make sure spkTimes are all of time values that are monotonically
    % increasing:
    spkTimes = sort(spkTimes);



    % Get PRE-DBS period spike rate, both observed and corrected for artifact
    % blanking

    % First check if this calculation has already been done and stored in
    % intermediate data
%     subFolder = 'FRpre';
%     if ~exist([scriptIntermediateDataFolder, '\', subFolder], 'dir')
%         mkdir([scriptIntermediateDataFolder, '\', subFolder]);
%         
%     end
    formatSpecFn = 'fr%s_%ss-%ss_fromDbsOnset';
    matfnStr = sprintf(formatSpecFn, objectID, ...
                              num2str(SPKINTERV(1)), ...
                              num2str(SPKINTERV(2)));
    fullPathFn = [scriptIntermediateDataFolder '\' matfnStr '.mat'];
    
    % load/calculate firing rates for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn)
        
    else
        % If SPKINTERV specifies times that fall outside of the recording,
        % adjust SPKINTERV to have boundaries within recorded times
        spkInterval = SPKINTERV; 
        tbegRef = nexFile.tbeg - StimTs.DBS(1);
        tendRef = nexFile.tend - StimTs.DBS(1);
        spkInterval(1) = max(spkInterval(1), tbegRef);
        spkInterval(2) = min(spkInterval(2), tendRef);
        
        % Get various data for firing rate, outputted in struct
        FR = calcFiringRates(spkTimes, StimTs, spkInterval);
        FR.intervalLimits = spkInterval;
        save(fullPathFn, 'FR');
    
    end
    
%     FRpreCorr(iNex,1) = FR.rateCorrected;

end