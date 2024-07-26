% test script for "psth.m"


stimFreq = 100; % Hz
dt = 1/stimFreq;


refs = [60:dt:120]'; % seconds

nEvents = 10000;
eventLims = [0, 150];

% randomly distributed "spike" events; a psth of this should look like a
% Uniform distribution.
events = eventLims(2) * rand(nEvents, 1);
events = sort(events);


binEdges = 0:0.001:dt;

probDistrib = psth(events, refs, binEdges);

figure; bar(probDistrib);


