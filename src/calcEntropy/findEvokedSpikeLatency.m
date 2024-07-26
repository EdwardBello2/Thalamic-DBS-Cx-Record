function spkTevoked = findEvokedSpikeLatency(spkTimes, stimTimes, win )
% Find all spike times of spikes evoked by a stim event. Function gets the
% FIRST spike detected after each stim, and falling within the time-bin
% specifed by the two values in "win". Outputs an N-length colmun, where N
% is the number of stimTimes. If no spikes were found in "win", then output
% for that row is NaN. 



nStims = numel(stimTimes);
spkTevoked = zeros(nStims, 1);
spkTevoked(:) = NaN;
for iStim = 1:nStims
    % Start by referencing all spike times to the current stim
    dTimes = spkTimes - stimTimes(iStim);
    spkLatency = [];
    
    
    isEvoked = (dTimes >= win(1)) & (dTimes < win(2));
    spkLatency = dTimes(isEvoked);

    
    if isempty(spkLatency)
        spkTevoked(iStim,1) = NaN;
        
    else
        spkTevoked(iStim,1) = min(spkLatency);
        
    end
        
    
end 

% spkTevoked(isnan(spkTevoked)) = [];




end