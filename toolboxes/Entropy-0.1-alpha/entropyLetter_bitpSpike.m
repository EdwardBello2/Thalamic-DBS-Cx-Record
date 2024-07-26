function HbitpSpk = entropyLetter_bitpSpike(distrib)
% Get the Entropy of a probability distribution using the "letter" method, 
% in bits/spike.
%
% SYNTAX:
% H = entropyLetter_bitpSec(distrib)
%


% Remove any parts of the discrete distribution that are zero-probability
distrib(distrib == 0) = [];

% Entropy in bits/spike
HbitpSpk = -sum((distrib.*log2(distrib)));




end