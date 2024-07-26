% test script for isiLogHistogram



load('16_03_25u_0001_proc.mat');

% get one trial's worth of spike data
trial = 3;

fs = F.s_rate;
spkTimeInterval = (F.stim(trial).interval) / fs;
dbsTimeInterval = [F.stim(trial).on(1), F.stim(trial).on(end)];

isWithinTrial = (F.CC1_n1 >= spkTimeInterval(1)) & ...
              (F.CC1_n1 < spkTimeInterval(2));          
spksInTrial = F.CC1_n1(isWithinTrial);


% get pre-DBS spikes
isPreDBS = spksInTrial < dbsTimeInterval(1);
preDBSspks = spksInTrial(isPreDBS);


% get DBS-on spikes
isDBSon = (spksInTrial >= dbsTimeInterval(1)) & ...
          (spksInTrial < dbsTimeInterval(2));
DBSspks = spksInTrial(isDBSon);



%%

binPD = 20;
scale = 1000; % to rescale bins from seconds to ms

[isiLogPre, binEdgesPre] = isiLogHistogram(preDBSspks, binPD);
figure;
h1 = histogram('BinEdges', (binEdgesPre * scale), 'BinCounts', isiLogPre)

hold on

[isiLogDBS, binEdgesDBS] = isiLogHistogram(DBSspks, binPD);
h2 = histogram('BinEdges', (binEdgesDBS * scale), 'BinCounts', isiLogDBS)

% h1.Normalization = 'probability';
% h2.Normalization = 'probability';

set(gca,'xscale','log');
legend('Pre-DBS', 'DBS-on');
xlabel('Interspike Interval (ms)');
ylabel('Count');
title(['Log-binned Interspike Interval of Trial ', num2str(trial)])






