function HbitpSec = entropyLetter_bitpSec(distrib, nSpks, nRef, binTime)
% Get the Entropy of a probability distribution using the "letter" method, 
% in bits/second.
%
% SYNTAX:
% H = entropyLetter_bitpSec(distrib, nSpks, nRef, binTime)
%


% Remove any parts of the discrete distribution that are zero-probability
distrib(distrib == 0) = [];

% Entropy in bits/spike
HbitpSpk = -sum((distrib.*log2(distrib)));

% Estimate of spikes/sec based only on PSTH
spkpSec = nSpks / (nRef * binTime);

% Get into units of bits/second
HbitpSec = HbitpSpk * spkpSec;



end