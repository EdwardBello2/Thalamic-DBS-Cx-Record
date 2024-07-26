% Find out if the number of spikes within a certain interval cross a
% Hz-threshold defined by user
%
%
% Syntax:
%
% ThrResult = detectIntervalSpkRateThresh(spkTimes, hzThresh)
% ThrResult = detectIntervalSpkRateThresh(spkTimes, hzThresh, intervals)
% ThrResult = detectIntervalSpkRateThresh(__, Name, Value)
%
%
% 
% Description:
%
% ThrResult = detectIntervalSpkRateThresh(spkTimes, hzThresh) Outputs a
% true/false value for whether the input events, spkTimes, have an average
% rate of occurrence above hzThresh, for the default interval of the first
% event to the last event in spkTimes. 
%
% ThrResult = detectIntervalSpkRateThresh(spkTimes, hzThresh, intervals)
% finds out the above and outputs Nx1 true/false vector for user-specified
% Nx2 time intervals within spkTimes. 
% 
% ThrResult = detectIntervalSpkRateThresh(__, Name, Value) specify options
% for the kind of rtheshold crossing that gets detected as true/false. 
%
%
%
% Input Arguments:
%
%  spkTimes - vector of event times in seconds
%  hzThresh - scalar value of threshold in Hz
% intervals - Nx2 matrix of values in seconds for N time intervals to
%             perform the test on
%
% 
%
% Name-Value Pair Arguments:
%
% 'ThreshCross' - Type of hzThresh crossing to check for
% 'above'(default)|'below'|'equalabove'|'equalbelow'
%   Type of threshold crossing test corresponds to a matlab boolean operator
%   test according to the following:
%        'above': spkRate > hzThresh
%        'below': spkRate < hzThresh
%   'equalabove': spkRate >=hzThresh
%   'equalbelow': spkRate <=hzThresh


% Author: Ed Bello
% Created: 4/24/2019
% 
% TO-DO

function ThrResult = detectIntervalSpkRateThresh(spkTimes, hzThresh, varargin)
%% INPUTS

p = inputParser;

% Check that "intervals" was entered is an Nx2 matrix and that all
% intervals show increasing time
checkIntervals = @(x) (size(x, 2) == 2) && ...
                      all(diff(x, [], 2) > 0);


% Specify inputs to parse
addRequired(p, 'spkTimes', @isnumeric);
addRequired(p, 'hzThresh', @isnumeric);

defaultIntervals = [spkTimes(1), spkTimes(end)];
addOptional(p, 'intervals', defaultIntervals, checkIntervals)

defaultThreshCross = 'above';
addParameter(p, 'ThreshCross', defaultThreshCross, @checkThreshCross)


% PARSE and ASSIGN values for code:
parse(p, spkTimes, hzThresh, varargin{:})

  intervals = p.Results.intervals;
ThreshCross = p.Results.ThreshCross;



%% CODE

% For each interval specified by user, find the average spike rate within,
% and see how it crosses hzThresh:

nIntervals = size(intervals, 1);
ThrResult = false(nIntervals, 1);

for iInt = 1:nIntervals
    
    isWithinInt = (spkTimes >=intervals(iInt,1)) & ...
                  (spkTimes < intervals(iInt,2));
    spksWithinInt = spkTimes(isWithinInt);

    timeInt = intervals(iInt,2) - intervals(iInt,1);

    intervalSpkRate = numel(spksWithinInt) / timeInt;

    switch ThreshCross
        case 'above'
            ThrResult(iInt,1) = (intervalSpkRate > hzThresh);

        case 'below'
            ThrResult(iInt,1) = (intervalSpkRate < hzThresh);

        case 'equalabove'
            ThrResult(iInt,1) = (intervalSpkRate >= hzThresh);

        case 'equalbelow'
            ThrResult(iInt,1) = (intervalSpkRate <= hzThresh);

        otherwise
            error('ThreshCross option not valid!')

    end
        
    
end
        
    

end

function TF = checkThreshCross(x)
validThreshCross = {'above', 'below', 'equalabove', 'equalbelow'};
TF = false;

if ~any(strcmp(x, validThreshCross))
    error('Valid string inputs: above | below | equalabove | equalbelow')
    
else
    TF = true;

end

end