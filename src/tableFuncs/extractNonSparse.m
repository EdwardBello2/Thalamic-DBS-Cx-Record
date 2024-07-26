function Extracted = extractNonSparse(NexTable, timePRE, timeDBS, HzThresh)
% Looking at NEXprocfiles_"name", extract only those rows that have
% spk-rates above HzThresh in both preDBS and DBS-on condition.


nRows = size(NexTable,1);
extractRows = true(nRows, 1);

for iRow = 1:nRows
    nexfn = NexTable.Filename{iRow};
    nexpn = NexTable.Pathname{iRow};

    nexFile = readNexFile([nexpn, '\', nexfn]);
    [spkTimes, StimTs] = parseNexFile(nexFile);


    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = StimTs.DBS;
    stimPeriod = median(diff(dbsTimes));


    % get pre-DBS spikes
    isPreDBS = (spkTimes < dbsTimes(1)) & ...
               (spkTimes >= dbsTimes(1) - timePRE);
    spksPRE = spkTimes(isPreDBS);


    % get DBS-on spikes
    isDBSon = (spkTimes >= dbsTimes(1)) & ...
              (spkTimes < (dbsTimes(end) + stimPeriod));
    spksDBS = spkTimes(isDBSon);


    % Check if any section has spikes under 2Hz rate
    spksPRErate = numel(spksPRE) / timePRE;
    spksDBSrate = numel(spksDBS) / timeDBS;

    if (spksPRErate < HzThresh) || (spksDBSrate < HzThresh)
        extractRows(iRow) = false;
    end

end

Extracted = NexTable(extractRows,:);

end