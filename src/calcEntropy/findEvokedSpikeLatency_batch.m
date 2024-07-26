function EvokedSpikeLatencyResults = findEvokedSpikeLatency_batch(NEX_2analyze, pipeParams)
% Perform the "findEvokedSpikeLatency.m" function every NEX file specified
% in NEX-2analyze table, and append the resulting latency data to an
% updated table "EvokedSpikeLatencyResults". If the user specifed 
% "pipeParams.trimPSTH = true", then the first bin in the PSTH will be
% ignored. 


N = NEX_2analyze;
nNex = size(N, 1);

latency_mean = zeros(nNex, 1);
latency_std = zeros(nNex, 1);
latency_nSpks = zeros(nNex, 1);
nStims = zeros(nNex, 1);

tic
for iNex = 1:nNex

    % Get each NEX file from table
    nexfn = N.Filename{iNex,1};
    nexpn = N.Pathname{iNex,1};

    nexFile = readNexFile([nexpn, '\', nexfn]);
    [spkTimes, StimTs] = parseNexFile(nexFile);

    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = StimTs.DBS;
    stimPeriod = median(diff(dbsTimes));


    % get DBS-on spikes
    isDBSon = (spkTimes >= dbsTimes(1)) & ...
              (spkTimes < (dbsTimes(end) + stimPeriod));
    spksDBS = spkTimes(isDBSon);





%     % Find the first spike after each Stimulation 
%     % (estimates spike-latency after each pulse)
% 
%     if pipeParams.trimPSTH
%         edgeLeft = pipeParams.binEdgeLeft + pipeParams.binWidth;
% 
%     else
%         edgeLeft = pipeParams.binEdgeLeft;
% 
%     end
%     
%     switch pipeParams.latencyWin
%         case 'binEdgeRight'
%             edgeRight = pipeParams.binWidth * pipeParams.binNum;
% 
%         otherwise
%             edgeRight = stimPeriod;
% 
%     end

    
    win = pipeParams.latencyWin;

    spkLatencies = findEvokedSpikeLatency(spksDBS, dbsTimes, win);
    spkLatencies(isnan(spkLatencies)) = []; % Remove any NaN values

    % Fill the latency values for each NEX file
       latency_mean(iNex,1) = mean(spkLatencies);
        latency_std(iNex,1) = std(spkLatencies);
    latency_nSpks(iNex,1) = numel(spkLatencies); 
%     nStims(iNex,1) = numel(dbsTimes);

end % END for-loop
toc





    
    
   
% toc
disp('Done Calculating!')

% % halelujah!
% load handel;
% player = audioplayer(y, Fs);
% play(player);
% 


%% Create RESULTS table and SAVE

Lmean = table(latency_mean);
 Lstd = table(latency_std);
LnSpk = table(latency_nSpks);
% NSTIM = table(nStims);

% E = [N, Lmean, Lstd, LnSpk, NSTIM];
E = [N, Lmean, Lstd, LnSpk];


% specify some metadata for the resulting table:
E.Properties.Description = 'Table of individual DBS trial spike-recordings, with stim-evoked spike latency data appended.';
% E.Properties.VariableUnits = {'', '', '', '', '', 'DBS-electrode label', ...
%                               'Hz', '', 'bits/spk', 'bits/spk', '', '', ...
%                               '', 'seconds', 'seconds',}; 

EvokedSpikeLatencyResults = E;


end