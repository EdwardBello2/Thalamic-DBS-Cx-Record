% Entropy distributions

% load NEXprocfiles_Uva
load('K:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables\NEXprocfiles_Uva.mat')
% row 439 has a good example orthodromic: C0, 100Hz, Freq Sweep

row = 439;

nexfn = NEXprocfiles.Filename{row};
nexpn = NEXprocfiles.Pathname{row};

%% Get pre and DBS spike times

PREDBS_TIME = 60; % seconds

nexFile = readNexFile([nexpn, '\', nexfn]);
[spkTimes, StimTs] = parseNexFile(nexFile);

% separate the spike times into pre-DBS and DBS-on
dbsTimes = StimTs.DBS;
stimPeriod = median(diff(dbsTimes));


% get pre-DBS spikes
isPreDBS = (spkTimes < dbsTimes(1)) & ...
           (spkTimes >= (dbsTimes(1) - PREDBS_TIME));
spksPRE = spkTimes(isPreDBS);


% get DBS-on spikes
isDBSon = (spkTimes >= dbsTimes(1)) & ...
          (spkTimes < (dbsTimes(end) + (7.5/1000)));
spksDBS = spkTimes(isDBSon);



%% See what bootstrap resampling looks like for DBS-on spikes


% Use initial PSTH to get H value and true psth spike count:
binW = 0.5/1000; % seconds
binMax = 7.5/1000; % seconds
binEdges = 0:binW:binMax;
resampNum = numel(spksDBS);

psthCt = psth(spksDBS, StimTs.DBS, binEdges, 'count');
nSpks = sum(psthCt);

psthDBS = psthCt/sum(psthCt);
H_DBS = entropyLetter_bitpSpike(psthDBS);


% Now generate bootstrapped resampled psths
bootstrap.number = 10000;
bootstrap.psthResampN = 3000;


psthRows = psthBootstrap(spksDBS, StimTs.DBS, binEdges, bootstrap);

nRow = size(psthRows, 1);
H_bootDBS = zeros(nRow, 1);
for iRow = 1:nRow
    
    H_bootDBS(iRow) = entropyLetter_bitpSpike(psthRows(iRow,:));

end

figure;
histogram(H_bootDBS)

spks = bootstrap.psthResampN
title(['spks: ', num2str(spks), ', mean: ', num2str(mean(H_bootDBS)), ', std: ', num2str(std(H_bootDBS))])


%% See what bootstrap resampling looks like for PRE-dbs spikes


% Use initial PSTH to get H value and true psth spike count:
binW = 0.5/1000; % seconds
binMax = 7.5/1000; % seconds
binEdges = 0:binW:binMax;

psthCt = psth(spksPRE, StimTs.VirtPre, binEdges, 'count');
nSpks = sum(psthCt);

psthPRE = psthCt/sum(psthCt);
H_PRE = entropyLetter_bitpSpike(psthPRE);


% Now generate bootstrapped resampled psths
bootstrap.number = 10000;
bootstrap.psthResampN = 1200;


psthRows = psthBootstrap(spksPRE, StimTs.VirtPre, binEdges, bootstrap);

nRow = size(psthRows, 1);
H_bootDBS = zeros(nRow, 1);
for iRow = 1:nRow
    
    H_bootPRE(iRow) = entropyLetter_bitpSpike(psthRows(iRow,:));

end


figure;
histogram(H_bootPRE)

spks = bootstrap.psthResampN
title(['spks: ', num2str(spks), ', mean: ', num2str(mean(H_bootPRE)), ', std: ', num2str(std(H_bootPRE))])


%%

H_bootDBS = zeros(bootNum, 1);
for iBoot = 1:bootNum
    psthBoot = psthResamp(spksDBS, StimTs.DBS, binEdges, resampNum);
    
    H_bootDBS(iBoot) = entropyLetter_bitpSpike(psthBoot);
    
end













