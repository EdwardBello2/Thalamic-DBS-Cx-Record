% filter out unwanted data from the final analysis for ISI-based Entropy
%
% 1) Selected data only for Single Unit neurons (SU)
% 2) Remove any trials from analysis with spk-rate < 2Hz
%
% INPUT: only change subjID needed
% other inputs like pn may need to change depending on the location of the
% project folder on your PC


function Selected = createSelectedDataTable(NEXprocfiles, SortedUnits, pipeParams)

%% 1)

% Use the SortedUnits table to choose only the NEXprocfiles trials that
% belong to Neurons of type specified by neuTypeFilter
NEX_byType = selectDataRows_byNeuType(NEXprocfiles, SortedUnits, pipeParams.neuTypeFilter);

% disp([pipeParams.neuTypeFilter, 's selected from NEXprocfiles table'])



%% 2)

% script assumes that a) DBS took exactly 60 seconds each time and b) that
% you want to look at only the preceding 60 seconds of spikes before DBS


PREDBS_TIME = pipeParams.predbsTime; % seconds
DBS_TIME = pipeParams.dbsTime; % seconds
HZ_THRESH = pipeParams.hzThresh; % Hz

% whittle down the rows further by choosing only trials with spk rates
% above 2 hz.
Selected = extractNonSparse(NEX_byType, PREDBS_TIME, DBS_TIME, HZ_THRESH);





end



