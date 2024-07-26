% script for calculating firing-rate differences in Cortex before and
% during DBS:

% Include code for manipulating nexFiles:
addpath(genpath('C:\Users\bello043\OneDrive\Tools\NEXtools'));

% Load one NEX file, with has ONE DBS trial of data for ONE unit
pn = 'D:\PROJECTS\Thalamic DBS Cx Record\DataProcessing\NEX_proc';
nexfn = '17080401block01-ch01a.nex';
nexFile = readNexFile([pn, '\', nexfn]);



% Get the spike times and DBS times
[spkTimes, stimTimes] = parseNexFile(nexFile);
dbsTimes = stimTimes.DBS;



%% Collect 1-sec bins of spike counts for each condition:

nBins = 60;
rateWin = 1; % seconds

ratesPRE = getRateBinsPreDBS(spkTimes, dbsTimes(1), 1, 60);

ratesDBS = getRateBinsDurDBS(spkTimes, dbsTimes(1), 1, 60);



figure; histogram(ratesPRE); hold on
histogram(ratesDBS);
legend('PRE', 'DBS'); title('Histograms of rates');

figure; plot(ratesPRE); hold on;
plot(ratesDBS);
legend('PRE', 'DBS'); title('Rates vs. time')


[p, h, stats] = ranksum(ratesPRE, ratesDBS)




