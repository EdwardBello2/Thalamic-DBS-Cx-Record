function PSTH = psthBootstrap(evTimes,refTimes,binEdges, nBoot, nSpks, varargin)
% Peri-stimulus (Peri-event) time histogram, resampling peristimulus spike
% times (with replacement). 
%
% SYNTAX:
% [Psth] = psth(evTimes, refTimes, binEdges, resampNum,)
% [Psth] = psth(evTimes, refTimes, binEdges, resampNum, psthMode)
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
%

% TO-DO
% - add in inputParser so that we have a name-value option to trim the PSTH
% bins

%% DEFAULT VALUES

% tic
if numel(varargin) == 1
    psthOutput = varargin{1};
else
    psthOutput = 'probability';
end

% toc



%% Create re-referenced versions of evTimes

% tic
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

% toc



%% Remove any evTimes outside the range of bins, and align evTimes to nearest ref event

% tic
evTimesRefNaN = evTimesRef;
evTimesRefNaN(evTimesRefNaN < binEdges(1)) = NaN;

psthTimes = min(evTimesRefNaN, [], 2);

% remove any spike time events that lie outside of the binranges around
% each reference stim event.
psthTimes(isnan(psthTimes)) = [];

isGreater = psthTimes > binEdges(end);
psthTimes(isGreater) = [];

% toc



%% GENERATE rows of bootstrapped psths

% pre-allocate psth rows matrix
% nBoot = bootstrap.number;
nBins = numel(binEdges) - 1;

psthCount = zeros(nBoot, nBins);

% tic
% nSpks = bootstrap.psthResampN;
for iBoot = 1:nBoot
    % get resampled version of psthTimes
    psthResamped = datasample(psthTimes, nSpks, 'Replace', true);


    % generate psth for resampled psthTimes
    [psthCount(iBoot,:),~] = histcounts(psthResamped, binEdges);

end
% toc



%% Final psth 


switch psthOutput
    
    case 'probability'
        PSTH = psthCount / nSpks;
        
    case 'count'
        PSTH = psthCount;
        
    otherwise
        error('unknown psth output type specified')
        
end



end % END function
