function event = nexGetEvent(nexFile, eventLabel)
%% INPUTS

p = inputParser;

% Specify inputs to parse
addRequired(p, 'nexFile', @isstruct);
addRequired(p, 'eventLabel', @ischar);



%% CODE

% eventStr = convertCharsToStrings(eventLabel);

nEv = size(nexFile.events, 1);

for iEv = 1:nEv
    evLabels{iEv,1} = nexFile.events{iEv,1}.name;
end


% extract only the specified event
evIdx = strcmp(eventLabel, evLabels);

if ~any(evIdx)
    evLabels
    error(['Your label: ', eventLabel, ' does not match any of the above labels in the nexFile struct.'])
    
end

event = nexFile.events{evIdx};

% % gather all DBS stim times (including fake virtual stims)
% dbsEv = find(strcmp(evLabels,'DBS_stims'));
% StimTimes.DBS = nexFile.events{dbsEv}.timestamps;
% 
% preEv = find(strcmp(evLabels,'VirtPre_stims'));
% StimTimes.VirtPre = nexFile.events{preEv}.timestamps;
% 
% posEv = find(strcmp(evLabels,'VirtPost_stims'));
% StimTimes.VirtPost = nexFile.events{posEv}.timestamps;
% 
% 
% % there should only be one neuron in these NEX files...      
% spkTimes = nexFile.neurons{1,1}.timestamps;





end