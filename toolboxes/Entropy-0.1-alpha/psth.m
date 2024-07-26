function [PSTH] = psth(evTimes,refTimes,binEdges, varargin)
% Peri-stimulus (Peri-event) time histogram. 
%
% SYNTAX:
% [Psth] = psth(evTimes, refTimes, binEdges)
% [Psth] = psth(evTimes, refTimes, binEdges, psthMode)
%
% DESCRIPTION:
% Creates a peth (psth) of one event's timestamps (evTimes) around
% another's (refTimes). Time bins are specified with the last input, where
% binEdges(1) can be a negative number to the left of the reference event.
% specify psthMode as either 'probability' or 'count', probability is
% default.
%
% THIS FUNCTION PRESUMES THAT ALL REFTIME EVENTS ARE IN ORDER OF EARLIEST 
% TO LATEST. 

% Author: Ed Bello
% Created: 1/20/2018

%% DEFAULT VALUES

if numel(varargin) == 1
    psthOutput = varargin{1};
else
    psthOutput = 'probability';
end



%% Create re-referenced versions of evTimes

% make sure evTimes is column-vector
if size(evTimes, 2) > size(evTimes, 1)
    evTimes = evTimes';
end

% make sure refTimes is row-vector
if size(refTimes, 2) < size(refTimes, 1)
    refTimes = refTimes';
end 


nEv = numel(evTimes);
nRef = numel(refTimes);

% Create N = nRef columns of evTimes
ev_nRefCols = evTimes * ones(1, nRef);

% Create M = nEv rows of refTimes
ref_nEvRows = ones(nEv, 1) * refTimes;

% Now get refTime-referenced evTimes columns, one for each refTime 
evTimesRef = ev_nRefCols - ref_nEvRows;



%% Remove any evTimes outside the range of bins, and align evTimes to nearest ref event
evTimesRefNaN = evTimesRef;
evTimesRefNaN(evTimesRefNaN < binEdges(1)) = NaN;

psthTimes = min(evTimesRefNaN, [], 2);

% remove any spike time events that lie outside of the binranges around
% each reference stim event.
psthTimes(isnan(psthTimes)) = [];
psthTimes(psthTimes >= binEdges(end)) = [];



%% Final psth 

[psthCount,~] = histcounts(psthTimes, binEdges);

switch psthOutput
    
    case 'probability'
        PSTH = psthCount / (sum(psthCount));
        
    case 'count'
        PSTH = psthCount;
        
    otherwise
        error('unknown psth output type specified')
        
end



end % END function
