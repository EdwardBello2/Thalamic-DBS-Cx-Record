% Convert output of histcounts to x and y vectors that show up as a stair
% function when plotted using "plot" and "plot3" (and possibly others).
% Inputs "N" and "edges" correspond exactly to the outputs of histcounts,
% but user can put their desired binned data in the same format.
%
%
% Syntax:
% [xS, yS] = histcounts2stairplot(N, edges)
%
% 
% Description:
% [xS, yS] = histcounts2stairplot(N, edges) converts histogram
% count/percentage data (1 x n) with accompanying bin edges (1 x n+1) into
% stair-wise data by repeating data values and interleaving them in such a 
% way that they will plot as a stair-wise function on many Matlab plotting 
% functions.
%
%
% INPUTS:
%           N - 1 x n vector of counts/percentages for each bin
%       edges - 1 x n+1 vector of bin edges corresponding to data in N
%
% OUTPUTS:
%          xS - stair-wise representation of the bin edges data
%          yS - stair-wise representation of the counts/percentages data


% Author: Ed Bello
% Created: 6/20/2019
% 
% TO-DO
%
function [xS, yS] = histcounts2stairplot(N, edges)

if size(edges, 1) > 1 % check that edges is row vector
    edges = edges';
    
end
x = [edges(1:end-1); edges(2:end)];
xS = reshape(x, numel(x), 1);

if size(N, 1) > 1 % check that N is row vector
    N = N';
    
end

if (numel(edges) - 1) ~= numel(N)
    error('an input of "n" bins must be accompanied by "n+1" edges');
    
end

y = [N; N];
yS = reshape(y, numel(y), 1);



end