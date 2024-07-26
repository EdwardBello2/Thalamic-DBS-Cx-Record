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


% divide spkTimes into pre and during DBS groups
timeWin = 60; %seconds
spksPRE = spkTimes(spkTimes < dbsTimes(1) & ...
                   spkTimes >= (dbsTimes(1) - timeWin));
               
spksDBS = spkTimes(spkTimes >= dbsTimes(1) & ...
                   spkTimes < (dbsTimes(1) + timeWin));



%% Display distribution of ISI-based Instantaneous frequencies 1) Before dBS
% and 2) During DBS
% isiPRE = diff(spksPRE);
% isiDBS = diff(spksDBS);
% 
% instFrPRE = 1 ./ isiPRE;
% instFrDBS = 1 ./ isiDBS;
% 
% figure; histogram(instFrPRE); hold on
% histogram(instFrDBS);
% legend ('PRE', 'DBS')


%% Collect 1-sec bins of spike counts for each condition:



nBins = 60;
rateWin = 1; % seconds

countsPRE = zeros(nBins, 1);
countsDBS = zeros(nBins, 1);

% First for DBS on times
binEdgeDBS = dbsTimes(1):rateWin:(dbsTimes(1) + rateWin*nBins);

for i = 1:nBins
    isIn_iBin = spkTimes >= binEdgeDBS(i) & ...
                spkTimes < binEdgeDBS(i+1);
                   
    countsDBS(i) = sum(isIn_iBin);
                   
end
ratesDBS = countsDBS / rateWin;
                   
 
% Next for PRE-DBS times
binEdgePRE = dbsTimes(1):-rateWin:(dbsTimes(1) - rateWin*nBins);
binEdgePRE = fliplr(binEdgePRE);

for i = 1:nBins
    isIn_iBin = spkTimes >= binEdgePRE(i) & ...
                spkTimes < binEdgePRE(i+1);
                   
    countsPRE(i) = sum(isIn_iBin);
                   
end
ratesPRE = countsPRE / rateWin;



figure; histogram(ratesPRE); hold on
histogram(ratesDBS);
legend('PRE', 'DBS'); title('Histograms of rates');

figure; plot(ratesPRE); hold on;
plot(ratesDBS);
legend('PRE', 'DBS'); title('Rates vs. time')


[p, h, stats] = ranksum(ratesPRE, ratesDBS)




