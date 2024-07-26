function [binCounts] = binTimeEvents(spkTimes, binEdgeTimes)
% Gather bin counts of spike event times as they fall within specified
% bin-time edges.
%
% Syntax:
% binCounts = getRateBinsDurDBS(spkTimes, binEdgeTimes)
%
%
% Description:
 

% Author: Ed Bello
% Created: 2019/4/29

%% CODE

% countsDBS = zeros(nBins, 1);

nBins = length(binEdgeTimes) - 1;
binCounts = zeros(nBins, 1);
for i = 1:nBins
    isIn_iBin = spkTimes >= binEdgeTimes(i) & ...
                spkTimes < binEdgeTimes(i+1);
                   
    binCounts(i) = sum(isIn_iBin);
                   
end
% ratesDBS = countsDBS / rateWin;
                   


end




