% Get the logorithmically-binned histogram of inter-spike intervals for
% spike-time data.
%
% Syntax:
% [isiLogHist, binEdges, binIdx] = isiLogHistogram(X, binPD)
% [isiLogHist, binEdges, binIdx] = isiLogHistogram(X, [binEdgeLeft, binPD])
%
% 
% Description:
% [isiLogHist, binEdges, binIdx] = isiLogHistogram(X, binPD) computes the ISI-histogram of
% X (n x 1 array) using binPD bins per decade. The leftmost and rightmost
% bin-edges are estimated automatically.
%
% [isiLogHist, binEdges, binIdx] = isiLogHistogram(X, [binEdgeLeft, binPD]) computes the
% ISI-histogram of X (n x 1 array) using binPD bins per decade, with the
% leftmost bin edge specified by binEdgeLeft. 
%
%
% INPUTS:
%           X - vector of spike times (in seconds)
%       binPD - scalar indicating logarithmic bins per decade
% binEdgeLeft - scalar indicating the leftmost bin edge (shortest time in
%                seconds to include in binning
%
% OUTPUTS:
%  isiLogHist - Counts of ISIs found within each bin
%    binEdges - edges of logarithmic bins, units in seconds (num edges = (num bins + 1)
%      binIdx - index array same size as X whose elements are the bin 
%                indices for the corresponding elemsnts in X.

% Author: Ed Bello
% Created: 1/13/2019
% 
% TO-DO
% - need to update help stuff above

function [isiLogHist, binEdges, binIdx] = isiLogHistogram(spikeTimes, binParam)
%%

if ~isvector(spikeTimes)
    error('Input is not vector array')
end


% If defining the left-most bin edge as a value less than the smallest ISI,
% use this scale factor to get that smaller value
SCALE_ISI0 = 0.9; % 



isiTimes = diff(spikeTimes); 
isiMin = min(isiTimes);

while isiMin == 0 % Remove this ISI because it should not happen
    zeroInds = find(isiTimes == isiMin);
    spikeTimes(zeroInds+1) = []; 
    isiTimes = diff(spikeTimes); 
    isiMin = min(isiTimes);
    
end

isiMax = max(isiTimes);



%%  Decide on the leftmost bin edge based on user-input

if isscalar(binParam) % automatically choose binEdges based on bins/decade
    binsPerDecade = binParam;
    
    % Set to be just below the shortest ISI
    binEdgeLeft = SCALE_ISI0 * isiMin; 
 
elseif isvector(binParam) && numel(binParam) == 2 % left edge user-defined
    binsPerDecade = binParam(2);
    
    % left end of the log ISI specified by user
    binEdgeLeft = binParam(1); 
   
else
    error('Unknown input for "binParam"');
    
end
   


%% Set up logarithmically-spaced bin-edges

nBins = 1;
binEdgeRight = binEdgeLeft * 10^(nBins/binsPerDecade);

while binEdgeRight < isiMax
    
    nBins = nBins + 1; 
    binEdgeRight = binEdgeLeft * 10^(nBins/binsPerDecade);
    
end

% This contains all the bin edges, from the left-most to the right most 
binEdges = binEdgeLeft*10.^((0:nBins)/binsPerDecade); 



%% Create ISI histogram with specified bin edges

[isiLogHist, ~, binIdx]= histcounts(isiTimes, binEdges); 



end