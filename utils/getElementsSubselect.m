% Gets a subselection of values from a linear array, and handles the cases
% where the desired number of sub-selected elements may exceed the actual
% number of elements by simply outputting itself, and marking the situation
% with "outputState" variable as well. for "outputState": 0 = desired
% subselection exactly equals the number of elements in "x"; 1 = desired
% subselection is indeed smaller than the total number of available
% elements in "x"; -1 = desired number of subselected elements exceeded number of
% elements in "x". Using the extra parameter "direction", User may also
% specify if the subselection is extracted from left-to-right ('LtoR') or
% right-to-left ('RtoL'), left pertaining to position 1 in the index and
% right pertaining to position "end" in the index. 



function [select, outputState] = getElementsSubselect(x, nElementsDesired, varargin)
%% INPUTS
p = inputParser;

addRequired(p, 'bins');
addRequired(p, 'nDesiredBins');

defaultDirection = 'LtoR';
addParameter(p, 'direction', defaultDirection, @ischar)

parse(p, x, nElementsDesired, varargin{:});

direction = p.Results.direction;



%% CODE

% nBinsDesired = 30;
if nElementsDesired > numel(x)
    nElementsDesired = numel(x);
    outputState = -1;
    
elseif nElementsDesired < numel(x)
    outputState = 1;
    
else
    outputState = 0;
    
end

switch direction
    case 'RtoL'
        select = x((end-(nElementsDesired-1)):end);
        
    case 'LtoR'
        select = x(1:nElementsDesired);

end