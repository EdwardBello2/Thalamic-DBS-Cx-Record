function H = func_calcHnex(nexFile, objectID, SPKINTERV, ppar)
% output firing rate calculations; if already cached, simply load the data, or else perform the calculation... 

%% Assert that intermediate data cache folder exists
[scriptDirectoryFullPath, scriptName] = fileparts(mfilename('fullpath'));

% First make sure that this script has a folder within the project folder's
% intermediate data section
scriptIntermediateDataFolder = [ppar.projRootPath, '\', ppar.intDataRootFolder, '\', scriptName];
if ~exist(scriptIntermediateDataFolder, 'dir')
    mkdir(scriptIntermediateDataFolder)

end


%% Get spike data out of nexfile and set up data cache name

% extract spike times and DBS stim times from nexFile struct
[spkTimes, StimTs] = parseNexFile(nexFile);

% make sure spkTimes are all of time values that are monotonically
% increasing:
spkTimes = sort(spkTimes);


ORD_H = ppar.ORD_H;
BINS_PER_DECADE = ppar.BINS_PER_DECADE;

formatSpecFn = 'h%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD';
matfnStr = sprintf(formatSpecFn, objectID, ...
                          num2str(SPKINTERV(1)), ...
                          num2str(SPKINTERV(2)), ...
                          num2str(ORD_H), ...
                          num2str(BINS_PER_DECADE));
fullPathFn = [scriptIntermediateDataFolder '\' matfnStr '.mat'];


%% Calculate entropy, or load it if already cached previously 

% load/calculate entropy for this nexfile
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
        
        
        % Get spike times occurring within interval 
        refT = StimTs.DBS(1);
          
        [spkTimesDbs] = getIntervalEvents(spkTimes, refT, spkInterval);

        % Define ISIs, and remove any ISIs == 0
        isiDBS = diff(spkTimesDbs);
        isiDBS(isiDBS == 0) = [];
        nISI = numel(isiDBS);

        if nISI < 2
            H = NaN;

        else
            H = entropyISIdirect_multOrder(isiDBS, BINS_PER_DECADE, ORD_H);

        end
    %     entropyDirISI.dbsEmp = H_DBS;


        
        save(fullPathFn, 'H');
    
    end

   
    
end