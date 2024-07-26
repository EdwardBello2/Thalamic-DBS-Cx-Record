function [intervEvents] = getIntervalEvents(evTimes, refT, refInterval)
% Gets a subselection of point event times within the reference interval.
% refInterval is 1x2 vector and is relative to refT. For example, this
% function will take arbitrary spike times and grab a subselection of them
% that are within 0 to 30 seconds from refT. All units in Seconds. 
rerefTimes = evTimes - refT;

% get observed spike rate "FRobs" based on count and time duration
idxInInterv = (rerefTimes >= refInterval(1)) & (rerefTimes < refInterval(2));
rerefSubselect = rerefTimes(idxInInterv);
intervEvents = rerefSubselect + refT;



end