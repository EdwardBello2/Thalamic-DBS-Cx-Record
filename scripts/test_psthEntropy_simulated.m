% test script for comparing entropy calculations normalized by psth events


% % bin parameters
% bw = 0.5; 
% nBins = 15;
% binEdges = bw * 0:nBins; % this array will be nBins + 1 in length
% 
% % assumed FR of test neuron:
% fr = 10; % spikes/sec
% 
% 
% nBoots = 40;
% 
% % 
% lambda = 0.1;
% R = poissrnd(lambda, 1000000, 1);
% figure; histogram(R);
% 
% n = 10000;
% p = 0.9;
% isSpk = zeros(n, 1);
% for i = 1:n
%     isSpk(i) = binornd(1, p);
% 
%    
% end
% spkyes = 100 * sum(isSpk == 1) / n;
% disp(['% trials that spiked: ', num2str(spkyes)]);
% 
% 
% y = poisscdf(41, lambda)


%% simulate 60 seconds of spiking during DBS
% 1) simulated spike times are poisson-process, uncorrelated to DBS
% 

% sampling frequency
% generate a poisson-distrib spike-train:
fs = 24414.0625;
dt = 1/fs;
fr = 10; % hz
tSim = 60; % seconds
nTrials = 1;

% psth bin parameters
bw = 0.5 / 1000; % seconds
nBins = 15;
binEdges = bw * (0:nBins); % this array will be nBins + 1 in length
    

nSims = 1000;
freqs = [10, 20, 30, 50, 100, 130];
nFreqs = numel(freqs);

H_bitSpk = zeros(nSims, nFreqs);
H_bitSec = zeros(nSims, nFreqs);
H_bitPulse = zeros(nSims, nFreqs);
FR   = zeros(nSims, nFreqs);
NUMSPKS = zeros(nSims, nFreqs);
tic
for iSim = 1:nSims

    [spkMat, tVec] = poissonSpikeGen(fr, tSim, nTrials, dt);

    spkTimes = tVec(spkMat);
    % nSpks = numel(spkTimes);
    % simFR = nSpks / tSim
    % spkISI = diff(spkTimes);
    % figure; histogram(spkISI);


    for iFr = 1:nFreqs
        % generate DBS stim times
        DBSrate = freqs(iFr); % hz
        DBSdt = 1 / DBSrate;
        dbsTimes = 0:DBSdt:tSim;

        % get a psth
             
            
        % for each bootstrap PSTH: based on current "allowed" bins (trimming?), 
        % random-uniformly sample spike counts for the bins
%                 psth_iBOOTct = zeros(length(binEdgesPsth)-1, 1)


        psthCt = psth(spkTimes, dbsTimes, binEdges, 'count');
%         figure; plot(psthCt)


        % calculate bit/spk Entropy of PSTH
        psthProb = psthCt / sum(psthCt);
        NUMSPKS(iSim,iFr) = sum(psthCt);
        Hpsth = entropyLetter_bitpSpike(psthProb);
        H_bitSpk(iSim,iFr) = Hpsth;
        
        H_bit = Hpsth * NUMSPKS;

        psthSpksec = (psthCt / bw) / numel(dbsTimes);
    %     figure; plot(psthSpksec);

        % get bit/second Entropy by normalized PSTH-firing-rate
        psthTime = binEdges(end) - binEdges(1);
        psthFR = mean(psthSpksec);
        
        
        FR(iSim,iFr) = psthFR;
        
        % calc bit/sec Entorpy of PSTH
        H_bitSec(iSim,iFr) = Hpsth * psthFR;
        
        % calc bit/pulse Entropy of PSTH
        H_bitPulse(iSim,iFr) = Hpsth * psthFR * (1/freqs(iFr));

        % % did we observe the number of spikes we'd expect to within our psth?
        % sum(psthCt) / (psthTime * numel(dbsTimes))
    
    end

end
toc

freqStrs = {'10', '20', '30', '50', '100', '130'};

% show Entropy bits
figure; ax = axes;
boxplot(exp(H_bitSpk));
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/spike)')
title('simulated poisson spikes H bit/spike');

% show Entropy bit/spike
figure; ax = axes;
boxplot(H_bitSpk);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/spike)')
title('simulated poisson spikes H bit/spike');

% show Entropy bit/second
figure; ax = axes;
boxplot(H_bitSec);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/second)')
title('simulated poisson spikes H bit/sec');

% show psthFR
figure; ax = axes;
boxplot(FR);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('firing rate (spike/sec)')
title('simulated poisson spikes psthFR');

% show Entropy bit/pulse
figure; ax = axes;
boxplot(log(H_bitPulse));
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/pulse)')
title('simulated poisson spikes H bit/sec');


disp('done!')


% look at relationship between pairs of features
rangeFR = [min(min(FR)), max(max(FR))];
rangeBitSpk = [min(min(H_bitSpk)), max(max(H_bitSpk))];
rangeBitSec = [min(min(H_bitSec)), max(max(H_bitSec))];

% Bit/sec vs. Bit/spk
figure;
for iFr = 1:nFreqs
    ax = subplot(2, 3, iFr);
    scatter(H_bitSpk(:,iFr), H_bitSec(:,iFr));
    ax.XLim = rangeBitSpk;
    ax.YLim = rangeBitSec;
    
end

% Bit/sec vs. FR
figure;
for iFr = 1:nFreqs
    ax = subplot(2, 3, iFr);
    scatter(FR(:,iFr), H_bitSec(:,iFr));
    ax.XLim = rangeFR;
    ax.YLim = rangeBitSec;
    
end

% Bit/spk vs. FR
figure;
for iFr = 1:nFreqs
    ax = subplot(2, 3, iFr);
    scatter(FR(:,iFr), H_bitSpk(:,iFr));
    ax.XLim = rangeFR;
    ax.YLim = rangeBitSpk;
    
end


% Look at correlation between psth Entropy and num of spikes in psth
figure;
for iFr = 1:nFreqs
    scatter(NUMSPKS(:,iFr), H_bitSpk(:,iFr));
    hold on
%     ax.XLim = rangeFR;
%     ax.YLim = rangeBitSpk;
    
end
legend(freqStrs, 'Location', 'northeastoutside')


% for each column, rescale all Entropy values by zscore
H_bitSpk_Z = zscore(H_bitSpk);
H_bitSec_Z = zscore(H_bitSec);
FR_Z = zscore(FR);


% show Entropy bit/spike
figure; ax = axes;
boxplot(H_bitSpk);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/spike)')
title('simulated poisson spikes H bit/spike');


% show Entropy bit/Sec


% show Entropy bit/second
figure; ax = axes;
boxplot(H_bitSec_Z);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/second)')
title('simulated poisson spikes H bit/sec');

% show psthFR
figure; ax = axes;
boxplot(FR_Z);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('firing rate (spike/sec)')
title('simulated poisson spikes psthFR');




%%

% Version of simulation where I control for number of spikes in each psth,
% assuming a 10Hz firing cell with no preference in psth bin (uncorrelated
% to DBS)

fr = 10; % hz
tSim = 60; % seconds
timePsthWin = binEdges(end) - binEdges(1);
nSims = 1000;

for iSim = 1:nSims

%     [spkMat, tVec] = poissonSpikeGen(fr, tSim, nTrials, dt);
% 
%     spkTimes = tVec(spkMat);
    % nSpks = numel(spkTimes);
    % simFR = nSpks / tSim
    % spkISI = diff(spkTimes);
    % figure; histogram(spkISI);


    for iFr = 1:nFreqs
        % generate DBS stim times
%         DBSrate = freqs(iFr); % hz
%         DBSdt = 1 / DBSrate;
%         dbsTimes = 0:DBSdt:tSim;

        % get a psth
        unifRange = binEdges;
%         unifRange(trimIdx) = [];
        A = unifRange(1);
        B = unifRange(end);
        
        nPulses = freqs(iFr) * tSim;
        timePsthTot = timePsthWin * nPulses;
        nSpksTot = fr * tSim;
        nSpksPsth = nSpksTot * (timePsthTot / tSim);
        
            
        % for each bootstrap PSTH: based on current "allowed" bins (trimming?), 
        % random-uniformly sample spike counts for the bins
%                 psth_iBOOTct = zeros(length(binEdgesPsth)-1, 1)
        R = unifrnd(binEdges(1), binEdges(end), nSpksPsth, 1);
        [psthCt,~] = histcounts(R, binEdges);

%         psthCt = psth(spkTimes, dbsTimes, binEdges, 'count');
%         figure; plot(psthCt)


        % calculate bit/spk Entropy of PSTH
        psthProb = psthCt / sum(psthCt);
        NUMSPKS(iSim,iFr) = sum(psthCt);
        Hpsth = entropyLetter_bitpSpike(psthProb);
        H_bitSpk(iSim,iFr) = Hpsth;
        
        H_bit = Hpsth * NUMSPKS;

        psthSpksec = (psthCt / bw) / nPulses;
    %     figure; plot(psthSpksec);

        % get bit/second Entropy by normalized PSTH-firing-rate
%         psthTime = binEdges(end) - binEdges(1);
        psthFR = mean(psthSpksec);
        
        
        FR(iSim,iFr) = psthFR;
        
        % calc bit/sec Entorpy of PSTH
        H_bitSec(iSim,iFr) = Hpsth * psthFR;
        
        % calc bit/pulse Entropy of PSTH
        H_bitPulse(iSim,iFr) = Hpsth * psthFR * (1/freqs(iFr));

        % % did we observe the number of spikes we'd expect to within our psth?
        % sum(psthCt) / (psthTime * numel(dbsTimes))
    
    end

end
toc


freqStrs = {'10', '20', '30', '50', '100', '130'};

% show Entropy bits
figure; ax = axes;
boxplot(H_bitSpk);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/spike)')
title('simulated poisson spikes H bit/spike');

% % show Entropy bit/spike
% figure; ax = axes;
% boxplot(H_bitSpk);
% ax.XTickLabels = freqStrs;
% xlabel('DBS frequency')
% ylabel('Entropy (bit/spike)')
% title('simulated poisson spikes H bit/spike');

% show Entropy bit/second
figure; ax = axes;
boxplot(H_bitSec);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/second)')
title('simulated poisson spikes H bit/sec');

% % show psthFR
% figure; ax = axes;
% boxplot(FR);
% ax.XTickLabels = freqStrs;
% xlabel('DBS frequency')
% ylabel('firing rate (spike/sec)')
% title('simulated poisson spikes psthFR');

% show Entropy bit/pulse
figure; ax = axes;
boxplot(H_bitPulse);
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('Entropy (bit/pulse)')
title('simulated poisson spikes H bit/sec');

% show Entropy bit/pulse
figure; ax = axes;
boxplot(log(H_bitPulse));
ax.XTickLabels = freqStrs;
xlabel('DBS frequency')
ylabel('log-Entropy (bit/pulse)')
title('simulated poisson spikes H bit/sec');

