function H = Entropy_H1(spikeTimes, binPD)
% This function calculates the first order direct estimate of the entropy
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
[ISIcounts,~]= histc(ISI, bins); 
% The last bin contains the number of ISIs that have exactly the same value as ISIright: see 'binranges' in http://www.mathworks.com/help/matlab/ref/histc.html#inputarg_binranges
ISIcounts(end-1) = ISIcounts(end-1)+ ISIcounts(end); % Basically combine the last two bins
ISIcounts = ISIcounts(1:end-1); % Rid the last bin 
% Calculate the entropy 
% First find the bins that have zero count, in this case, P*logP = 0, and thus doesn't add anything to the entropy calculation 
ISIcounts = ISIcounts(find(ISIcounts~=0)); 
H = -sum(((ISIcounts/numel(ISI)).*(log2(ISIcounts/numel(ISI))))); % bits/ISI
end