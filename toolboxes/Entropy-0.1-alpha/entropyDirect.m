function H = entropyDirect(spikeTimes, binPD, order)
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
% H    - This is the Nth order direct entropy estimate from the spike-
%        train, unit is in (bits/ISI).
% ISI  - Interspike intervals from the input spike train
% bins - The bins in log time, given the binPD


[~, ~, binIdx] = isiLogHistogram(spikeTimes, binPD);

binBlocks = groupSubsequent(binIdx, order);

pJoint = isiJointProb(binBlocks); 

H = (-1 / order) * sum((pJoint .* log2(pJoint))); % bits/spike
end

function p = isiJointProb(isi)
% Find the frequency of occurrence of each unique "word" of 
% subsequent ISIs
% 
% Input: 
% 
% isi - isi block
% n        - number of adjacent vector elements to group together 
% 
% Output: 
%
% p - frequency of occurrence of each unique block, this is a vector 


% find unique rows and their unique id
[~,~,u_id] = unique(isi, 'rows'); 

% count occurrences of unique ids
occurrences = histc(u_id, unique(u_id)); 

p = occurrences/sum(occurrences);

end 
