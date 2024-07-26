function p = blockCount(inputVec, n)
% This function will find segment a vector into overlapping blocks of
% length n, consisted of adjacent vector elements. Then find the frequency
% of occurrence of each unique block 
% 
% Input: 
% 
% inputVec - input vector (1xn)
% n        - number of adjacent vector elements to group together 
% 
% Output: 
%
% p - frequency of occurrence of each unique block, this is a vector 

blockMat = zeros(n, length(inputVec) - n + 1);
for i = 1:n
    blockMat(i, :) = inputVec(i:end-n+i);
end 
blockMat = transpose(blockMat); % transpose to be able to find unique rows (not cols)
% solution: http://stackoverflow.com/questions/26860149/count-the-duplicate-columns-of-a-matrix-in-matlab
[~,~,u_id] = unique(blockMat, 'rows'); % // find unique rows and their unique id
occurrences = histc(u_id, unique(u_id)); % // count occurrences of unique ids
p = occurrences/sum(occurrences);

end 