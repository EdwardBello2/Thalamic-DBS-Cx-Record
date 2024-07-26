function bins = genIntData_spkRate_XXsec_XsecBins(spkTimes, dbsTimes, pipeParams)
% gets known spike times and DBS times, and counts all spikes as they fall
% into the specified time bins, divides by total bin time to get rate
% estimate within each bin.


binWidth = pipeParams.ratebinWidth; % second 
totTime = pipeParams.totTime; % seconds

tbeg = pipeParams.tbeg;
tend = pipeParams.tend;

stimPeriod = median(diff(dbsTimes));

% calc dbs-on binned spk rates
dbsOnset = dbsTimes(1);
binEdgesDBS = dbsOnset:binWidth:(dbsOnset + totTime);
binCtsDBS = binTimeEvents(spkTimes, binEdgesDBS);
binRatesDBS = binCtsDBS / binWidth;


% calc pre-dbs binned spk rates
tStartPRE = dbsOnset - totTime;
tStartPRE = max(tStartPRE, tbeg); % in case totTime exceeds amount of time in PRE condition
binEdgesPRE = dbsOnset:-binWidth:tStartPRE;
binEdgesPRE = fliplr(binEdgesPRE); % make sure the bin edges are in ascending order
binCtsPRE = binTimeEvents(spkTimes, binEdgesPRE);
binRatesPRE = binCtsPRE / binWidth;


% calc post-dbs binned spk rates
tBegPOS = dbsTimes(end) + stimPeriod;
tEndPOS = tBegPOS + totTime;
tEndPOS = min(tEndPOS, tend); % in case totTime exceeds amount of time in POS condition
binEdgesPOS = tBegPOS:binWidth:tEndPOS; % similar to dbs bins, but shifted 60 seconds further in time
binCtsPOS = binTimeEvents(spkTimes, binEdgesPOS);
binRatesPOS = binCtsPOS / binWidth;



% Save intermediate mat file to intermediate processing location, and track
% this saved location in an "Intermediate" table
bins.dbs = binRatesDBS;
bins.pre = binRatesPRE;
bins.pos = binRatesPOS;




end