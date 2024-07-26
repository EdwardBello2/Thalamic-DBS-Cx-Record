function Hboots = func_calcHnex_bootstrap(nexFile, objectID, SPKINTERVselect, NBOOTS, SPKINTERVmatch, ppar)
% output firing rate calculations; if already cached, simply load the data, or else perform the calculation... 

%% Assert that intermediate data cache folder exists
[~, scriptName] = fileparts(mfilename('fullpath'));

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

formatSpecFn = 'h%s_%ss-%ss_fromDbsOnset_ordH%s_%sbinsPD_matchToInterv_%ss-%ss';
matfnStr = sprintf(formatSpecFn, objectID, ...
                          num2str(SPKINTERVselect(1)), ...
                          num2str(SPKINTERVselect(2)), ...
                          num2str(ORD_H), ...
                          num2str(BINS_PER_DECADE), ...
                          num2str(SPKINTERVmatch(1)), ...
                          num2str(SPKINTERVmatch(2)));
fullPathFn = [scriptIntermediateDataFolder '\' matfnStr '.mat'];


%% Calculate entropy, or load it if already cached previously 

% load/calculate entropy for this nexfile
    if exist(fullPathFn, 'file')
        load(fullPathFn, 'Hboots')
        
    else
        % If SPKINTERV specifies times that fall outside of the recording,
        % adjust SPKINTERV to have boundaries within recorded times
        spkIntervalSelect = SPKINTERVselect; 
        tbegRef = nexFile.tbeg - StimTs.DBS(1);
        tendRef = nexFile.tend - StimTs.DBS(1);
        spkIntervalSelect(1) = max(spkIntervalSelect(1), tbegRef);
        spkIntervalSelect(2) = min(spkIntervalSelect(2), tendRef);
        
        % also perform for reference interval for bootstrap
        spkIntervalMatch = SPKINTERVmatch;
        spkIntervalMatch(1) = max(spkIntervalMatch(1), tbegRef);
        spkIntervalMatch(2) = min(spkIntervalMatch(2), tendRef);
        
        
        % Get spike times occurring within interval 
        refT = StimTs.DBS(1);
          
        spkTimesSelect = getIntervalEvents(spkTimes, refT, spkIntervalSelect);
        spkTimesMatch = getIntervalEvents(spkTimes, refT, spkIntervalMatch);
        

        % Define ISIs, and remove any ISIs == 0
        isiSelect = diff(spkTimesSelect);
        isiSelect(isiSelect == 0) = [];
        nISIselect = numel(isiSelect);
        
        isiMatch = diff(spkTimesMatch);
        isiMatch(isiMatch == 0) = [];
        nISImatch = numel(isiMatch);
                

         
        % Fill get bootstrapped values of entropy for the select interval
        Hboots = zeros(NBOOTS, 1);
        
        if (nISImatch < 2) || isempty(isiMatch) % if the match interval has too few spikes...
            Hboots(:) = NaN;
            
        elseif (nISIselect < 2) || isempty(isiSelect)% if the select interval itself has too few spikes...
            Hboots(:) = NaN;

        else % generate all NBOOTS bootstrapped entropy values
            for iBoot = 1:NBOOTS
                % Resample PRE ISIs (with replacement) so that number of ISIs
                % matches the DBS case
                [isiMatchResamp] = datasample(isiSelect, numel(isiMatch));

                % Calculate multi-order Direct-Entropy
                Hboots(iBoot,:) = entropyISIdirect_multOrder(isiMatchResamp, BINS_PER_DECADE, ORD_H);

            end

        end
    
    
        % Save results in cache for later retreival, safe a crap ton of
        % time!
        save(fullPathFn, 'Hboots');
    
    end

   
    
end