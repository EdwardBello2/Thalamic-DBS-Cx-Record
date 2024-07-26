function [Psth] = psth(evTimes,refTimes,binEdges)
% Peri-stimulus (Peri-event) time histogram. THIS FUNCTION PRESUMES THAT
% ALL REFTIME EVENTS ARE IN ORDER OF EARLIEST TO LATEST.
%
% SYNTAX:
% [Psth] = psth(evTimes,reftTimes,binEdges)
%
% DESCRIPTION:
% Creates a peth (psth) of one event's timestamps (evTimes) around
% another's (refTimes). Time bins are specified with the last input, where
% binEdges(1) can be a negative number to the left of the reference event. 
%
% Author: Ed Bello
% Created: 9/18/2018




if isempty(evTimes) || isempty(refTimes)
    Psth = NaN;
    return
end
    


% simply remove any event times that are too early or too late to be in the
% psth analysis
evTimes(evTimes < (refTimes(1) - binEdges(1))) = [];
evTimes(evTimes > (refTimes(end) + binEdges(end))) = [];



xMin = binEdges(1);
nBins = numel(binEdges) - 1;
nRefs = length(refTimes);
binCounts = zeros(nRefs, nBins);

for iRefs = 1:nRefs
    
    dTimes = evTimes - refTimes(iRefs);
    dTimes(dTimes < xMin) = [];
    
    for jBins = 1:nBins
        
        isInjBin = (dTimes >= binEdges(jBins)) & ...
                   (dTimes < binEdges(jBins + 1)) ;
        count = numel(dTimes(isInjBin));
        binCounts(iRefs, jBins) = count;
        
    end
    
end 



Psth = sum(binCounts,1);

end
