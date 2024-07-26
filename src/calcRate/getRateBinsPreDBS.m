function [ratesPRE] = getRateBinsPreDBS(spkTimes, dbsStart, rateWin, nBins)
% Gather bins of rate estimates for the spiking activity right BEFORE DBS
% starts.
%
% Syntax:
% ratesPRE = getRateBinsPreDBS(spkTimes, dbsStart, rateWin, nBins)
%
%
% Description:
% Get a collection of rate estimates for spkTimes (in seconds) before DBS,
% in order of earliest to latest. Specify time bin "rateWin" in seconds.
% nBins indicates the number of "rateWin"-second bins to collect
% spike-times for, going BACKWARDS in time from "dbsStart" (in seconds). 

% Author: Ed Bello
% Created: 2019/4/11


%% CODE 

binEdgePRE = dbsStart:-rateWin:(dbsStart - rateWin*nBins);
binEdgePRE = fliplr(binEdgePRE);

countsPRE = zeros(nBins, 1);

for i = 1:nBins
    isIn_iBin = spkTimes >= binEdgePRE(i) & ...
                spkTimes < binEdgePRE(i+1);
                   
    countsPRE(i) = sum(isIn_iBin);
                   
end
ratesPRE = countsPRE / rateWin;



end




