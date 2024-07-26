% Entropy thought experiments

%% Case: a steadily increasing peak in a distribution
% as more spikes added, spike-rate increases;
% how do bits/spike and bits/second scale with each other? 


spkNum = 10:100;

nRef = 100;
binTime = 10/1000; % seconds
spkRate = spkNum /(nRef * binTime);


nValues = numel(spkNum);

for i = 1:nValues    
    count = [1, 1, (spkNum(i)-9), 1, 1, 1, 1, 1, 1, 1];


    distrib = count/sum(count);
    HbitpSec(i) = entropyLetter_bitpSec(distrib, spkNum(i), nRef, binTime);
    HbitpSpk(i) = entropyLetter_bitpSpike(distrib);

end

figure;
% plot(spkNum, spkRate); hold on

plot(spkNum, HbitpSpk); hold on


plot(spkNum, HbitpSec);
legend('bits/spike', 'bits/second')
title('Entropy of a steadily increasing single sharp peak distrib')
ylabel('Bits per "X"'); xlabel('Spike-rate (Hz)')

%% Case: number of stims change -- almost perfect "antidromic"
% This is to test the case of comparing, say, 10Hz stim to 100Hz stim, to
% see how Entropy would change as number of stim events change. 
% Suppose that in 60 seconds of stim, a spike only spikes 600 times
% regardless whether 10Hz or 100Hz stim is used. It's PSTH is also the same
% in terms of probability. With PSTHs spike rate is estimated based just on
% the PSTH window length and the number of windows (i.e. number of stims).
% The experiment below shows what would happen if the 10Hz stim gave us
% almost perfect antidromic following -- at the same time the 100Hz case
% does not affect the total number of spikes seen, 
% 

spkNum = 600;

nRef = 600:6000;
binTime = 10/1000; % seconds
% spkRate = spkNum / (nRef * binTime);

nValues = numel(nRef);

for i = 1:nValues
    count = [1, 1, 10, 1, 1, 1, 1, 1, 1, 1];
    distrib = count / sum(count);
    
    HbitpSec(i) = entropyLetter_bitpSec(distrib, spkNum, nRef(i), binTime);
    HbitpSpk(i) = entropyLetter_bitpSpike(distrib);
    
    
    
end

figure;
plot(nRef, HbitpSpk); hold on
plot(nRef, HbitpSec);
legend('bits/spike', 'bits/second')
title('Entropy of a steadily increasing number of stims')
ylabel('Bits per "X"'); xlabel('Stim-Events')




%% see how bits per spike does when number of spikes is held constant across 10Hz and 100Hz stim
% Now we just scale bits/spike with the number of spikes found
spkNum = 600;

nRef = 600:6000;
binTime = 10/1000; % seconds

nValues = numel(nRef);

for i = 1:nValues
    count = [1, 1, 10, 1, 1, 1, 1, 1, 1, 1];
    distrib = count / sum(count);
    
    HbitpSec(i) = entropyLetter_bitpSec(distrib, spkNum, nRef(i), binTime);
    HbitpSpk(i) = entropyLetter_bitpSpike(distrib);
    
    
    
end

figure;
plot(nRef, HbitpSpk); hold on
plot(nRef, HbitpSec);
plot(nRef, HbitpSpk * spkNum);
legend('bits/spike', 'bits/second', 'bits')
title('Entropy of a steadily increasing number of stims')
ylabel('Bits per "X"'); xlabel('Stim-Events')







