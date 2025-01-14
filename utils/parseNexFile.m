function [spkTimes, StimTimes] = parseNexFile(nexFile)
% Get spike-time and DBS-time data out of the nexfile
%
% Syntax:
% 
% [spkTimes, StimTimes] = parseNexFile(nexFile)
%
% 
% Description:
%
% [spkTimes, StimTimes] = parseNexFile(nexFile) extracts the data from the
% nexFile object generated by the output of "readNexFile.m". spkTimes is a
% vector of spike-times, StimTimes is a struct containing vectors for DBS
% times as well as virtual stim times:
% StimTimes
%          .DBS      -- DBS pulse timestamps
%          .VirtPre  -- Virtual stim pulse timestamps in trial before DBS
%          .VirtPost -- Virtual stim pulse timestamps in trial after DBS 
%

nEv = size(nexFile.events, 1);

for iEv = 1:nEv
    eventLabels{iEv,1} = nexFile.events{iEv,1}.name;
end

% gather all DBS stim times (including fake virtual stims)
dbsEv = find(strcmp(eventLabels,'DBS_stims'));
StimTimes.DBS = nexFile.events{dbsEv}.timestamps;

preEv = find(strcmp(eventLabels,'VirtPre_stims'));
StimTimes.VirtPre = nexFile.events{preEv}.timestamps;

posEv = find(strcmp(eventLabels,'VirtPost_stims'));
StimTimes.VirtPost = nexFile.events{posEv}.timestamps;


% there should only be one neuron in these NEX files...      
spkTimes = nexFile.neurons{1,1}.timestamps;


end