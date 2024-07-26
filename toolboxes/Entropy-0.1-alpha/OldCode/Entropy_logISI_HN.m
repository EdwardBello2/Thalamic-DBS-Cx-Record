function [H, ISI, bins] = Entropy_HN(spikeTimes, binPD, order)
% This function calculates the Nth order direct estimate of the entropy
% from a spike train (spike times are in miliseconds)
%
% Input:
%
% spikeTimes - This is the spike times in (ms)
% binPD      - This is the number of bins per ISI decade
% order      - Order number
%
% Output:
%
% H    - This is the Nth order direct entropy estimate from the spike train,
%        unit is in (bits/ISI).
% ISI  - Interspike intervals from the input spike train
% bins - The bins in log time, given the binPD

ISI = diff(spikeTimes); 

% Perform log binning of ISIs
ISImin = min(ISI);
if ISImin == 0
    % Remove this ISI because it should not happen
    zeroInds = find(ISI == ISImin);
    spikeTimes(zeroInds+1) = []; 
    ISI = diff(spikeTimes); 
    ISImin = min(ISI);
end

ISImax = max(ISI);
% Need to round down to find the left end of the log ISI time axis
ISIleft = roundDown(ISImin);
ISIright = roundUp(ISImax); % This is the right end of the log ISI axis. This is set to be just above the longest ISI
decadeNum = log10(ISIright/ISIleft); % Number of decades btw ISImax and ISIleft
% Need to round up the decade number 
binNum = (decadeNum)*binPD; % This is the total number of bins
bins = ISIleft*10.^(linspace(0,decadeNum, binNum+1)); % This contains all the bin edges, from the left-most to the right most
[ISIcounts,ind]= histc(ISI, bins); % ind, an array the same size as x indicating the bin number that each entry in x sorts into. 
p_joint = blockCount(ind, order); 
H = -1/order*sum((p_joint.*log2(p_joint))); % bits/ISI
end