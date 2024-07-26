function H = Entropy_H2(spikeTimes, binPD)
% This function calculates the second order direct estimate of the entropy
% from a spike train (spike times are in miliseconds)
%
% Input:
%
% spikeTimes - This is the spike times in (ms)
% binPD     - This is the number of bins per ISI decade
%
% Output:
%
% H - This is the 1st order direct entropy estimate from the spike train,
%     unit is in (bits/ISI).  

ISI = diff(spikeTimes); 
% Perform log binning of ISIs
ISImin = min(ISI);
ISImax = max(ISI);
% Need to round down to find the left end of the log ISI time axis
ISIleft = roundDown(ISImin); 
ratio = ceil(log10(ISImax/ISImin)); % Round up the log ratio 
ISIright = ISIleft*(10^(ratio)); % This is the right end of the log ISI axis. This is set to be just above the longest ISI
binNum = (ratio)*binPD; % This is the total number of bins
bins = ISIleft*10.^([0:binNum]/binPD); % This contains all the bin edges, from the left-most to the right most 
[ISIcounts,ind]= histc(ISI, bins); 
p_joint = blockCount(ind, 2);
H = -1/2*sum((p_joint.*log2(p_joint))); % bits/ISI
end