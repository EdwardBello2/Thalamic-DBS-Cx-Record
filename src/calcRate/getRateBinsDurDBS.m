function [ratesDBS] = getRateBinsDurDBS(spkTimes, dbsStart, rateWin, nBins)
% Gather bins of rate estimates for the spiking activity DURING DBS.
%
% Syntax:
% ratesDBS = getRateBinsDurDBS(spkTimes, dbsStart, rateWin, nBins)
%
%
% Description:
% Get a collection of rate estimates for spkTimes (in seconds) after DBS
% begins, in order of earliest to latest. Specify time bin "rateWin" in 
% seconds. "nBins" indicates the number of "rateWin"-second bins to collect
% spike-times for, going FORWARDS in time from "dbsStart" (in seconds). 

% Author: Ed Bello
% Created: 2019/4/11

%% CODE

countsDBS = zeros(nBins, 1);

binEdgeDBS = dbsStart:rateWin:(dbsStart + rateWin*nBins);

for i = 1:nBins
    isIn_iBin = spkTimes >= binEdgeDBS(i) & ...
                spkTimes < binEdgeDBS(i+1);
                   
    countsDBS(i) = sum(isIn_iBin);
                   
end
ratesDBS = countsDBS / rateWin;
                   


end




