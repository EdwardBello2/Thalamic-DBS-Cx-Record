function H = entropyISIdirect_multOrder(isi, binPD, nOrders)
% This function calculates the Nth order direct estimates of the entropy
% from a spike train (spike times are in miliseconds)
%
% Input:
%
% isi        - This is the inter-spike interval times in (ms)
% binPD      - This is the number of bins per ISI decade
% nOrders    - Order numbers such that 1 to nOrders will be computed
%
% Output:
%
% H    - Vector of  Nth order direct entropy estimates from the spike-
%        train, unit is in (bits/spike). Nth order is H(N).
% 


[~, binEdges, binIdx] = isiLogBinned(isi, binPD);

H = zeros(1, nOrders);

for iOrd = 1:nOrders

binBlocks = groupSubsequent(binIdx, iOrd);

pJoint = isiJointProb(binBlocks); 

H(1,iOrd) = (-1 / iOrd) * sum((pJoint .* log2(pJoint))); % bits/spike
end

end

function p = isiJointProb(isiWords)
% Find the frequency of occurrence of each unique "word" of 
% subsequent ISIs
% 
% Input: 
% 
% isiWords - matrix containing rows of isi "words"
% n        - number of adjacent vector elements to group together 
% 
% Output: 
%
% p - frequency of occurrence of each unique block, this is a vector 


% find unique row-words and their unique id
[~,~,u_id] = unique(isiWords, 'rows'); 

% a vector of "word-labels" for word-1, word-2, etc
uniqueWords = unique(u_id);

nWords = numel(uniqueWords);
occurrences = zeros(1, nWords);
for iWord = 1:nWords
    occurrences(1,iWord) = sum(u_id == uniqueWords(iWord));
   
end

% count occurrences of unique row-words
occurrences = histc(u_id, unique(u_id)); 


p = occurrences / sum(occurrences);

end 
