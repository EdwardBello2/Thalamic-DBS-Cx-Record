function  [eEPF, pfs, pfsb] = calcEEPF(nexFile)
% test script for classifying PSTH's into anti, ortho, and other
% based on 1) first peak latency, 2) peak jitter



%% CONSTANTS:

% number of standard deviations above mean firing rate to set a threshold
% for PSTH peak to cross
numSTD = 6; 

% gather PSTHs:
bw = 0.1 / 1000; % seconds
bMax = 7.5 / 1000; % seconds
binEdges = 0:bw:bMax;
t = binEdges(1:end-1) * 1000; % ms



%%

% % for one NEX file:
% pn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\Manuscripts\WorkingTitle\CURRENT DRAFT\Figures_Versions\Fig4\ExampleNexFiles';
% nexFolder = 'PhaseLock-H_decr';
% 
% fn = '17080401block12-ch12a';
% 
% nexFile = readNexFile([pn, '\', nexFolder, '\', fn, '.nex']);




[spkTimes, stims] = parseNexFile(nexFile); % outputs timestamps for spiketimes and DBS times (seconds)


% for PRE-dbs
PSTHcount_pre = psth(spkTimes, stims.VirtPre, binEdges, 'count');

PSTHratenorm_pre = (PSTHcount_pre / bw) / numel(stims.VirtPre); % rate-normalized PSTH
FRmean_pre = mean(PSTHratenorm_pre);
FRmean_preSTD = std(PSTHratenorm_pre);


% for DBS-on
PSTHcount_dbs = psth(spkTimes, stims.DBS, binEdges, 'count');
PSTHratenorm_dbs = (PSTHcount_dbs / bw) / numel(stims.DBS); % rate-normalized PSTH

% FRmean = mean(PSTHratenorm_dbs);
% FRthreshold = FRmean + (0.66 * (FRpeak - FRmean));
FRthreshold = FRmean_pre + (numSTD * FRmean_preSTD);

figure; plot(t, PSTHratenorm_dbs);
hold on; plot([t(1), t(end)], [FRthreshold, FRthreshold], '--r')



%% Detect FRpeak and peakJitter

[FRpeak, peakIdx] = max(PSTHratenorm_dbs);

tPeakLatency = t(peakIdx)

% find left and right borders above FRtrehshold
isAboveThresh = PSTHratenorm_dbs >= FRthreshold;
tPhaseLock = t(isAboveThresh);
tPeakLeft = min(tPhaseLock);
tPeakRight = max(tPhaseLock);
peakJitter = tPeakRight - tPeakLeft



%% Detect Excitatory Effective Pulse Fraction

% based on Filippo's 2015 paper

if isempty(tPeakLeft) & isempty(tPeakRight)
    pfs = 0;
    pfsb = 0;
    
else
    pfs = psth(spkTimes, stims.DBS, [tPeakLeft, tPeakRight], 'count');
    pfsb = psth(spkTimes, stims.VirtPre, [tPeakLeft, tPeakRight], 'count');

end

eEPF = (pfs - pfsb) / (numel(stims.DBS) - pfsb);






























end