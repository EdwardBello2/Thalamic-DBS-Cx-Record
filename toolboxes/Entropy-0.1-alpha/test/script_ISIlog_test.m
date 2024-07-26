% test script for my ISI Entropy code!


% load the NEX file of interest


% clear; close all

% locPath = 'D:\PROJECTS\Thalamic DBS Cx Record\DataProcessing';
% folder = 'NEX_proc';
% filename = [folder, '_Frequency-Sweep_C0-30Hz.nex'];
nexFile = readNexFile(['D:\PROJECTS\Thalamic DBS Cx Record\DataProcessing\NEX_proc\17080401block08-ch02a.nex']);

T = 1/100;

spkTimes = nexFile.neurons{1,1}.timestamps;
dbsTimes = nexFile.events{5,1}.timestamps;


predbsSpkTimes = spkTimes(spkTimes<dbsTimes(1));
durdbsSpkTimes = spkTimes(spkTimes>dbsTimes(1) & spkTimes<dbsTimes(end));

predbsISIs = diff(predbsSpkTimes);
durdbsISIs = diff(durdbsSpkTimes);

binwidth = 1/1000; % seconds
binedges = 0:binwidth:(300/1000);

predbsISIhist = histcounts(predbsISIs, binedges);
durdbsISIhist = histcounts(durdbsISIs, binedges);



H1pre = entropyH1(predbsISIhist);

%%
dim = 1:6;
for d = dim
    Hpre(d) = Entropy_logISI_HN(predbsSpkTimes/1000, 20, d);
end
Hpre = fliplr(Hpre);
dimrecip = 1./(fliplr(dim));

figure; plot(dimrecip, Hpre, '-o'); hold on

coeff = polyfit(dimrecip(end-1:end), Hpre(end-1:end), 1);
fitLine = polyval(coeff, [0, dimrecip]);
plot([0,dimrecip], fitLine)
title(['Direct Estimate of Entropy: ', num2str(fitLine(1)), ' bits/spk' ])
ylabel('H: bits/spk');
xlabel('Reciprocal dimension');


%%
dim = 1:6;
for d = dim
    Hdur(d) = Entropy_logISI_HN(durdbsSpkTimes/1000, 20, d);
end
Hdur = fliplr(Hdur);
dimrecip = 1./(fliplr(dim));

figure; plot(dimrecip, Hdur, '-o'); hold on

coeff = polyfit(dimrecip(end-1:end), Hdur(end-1:end), 1);
fitLine = polyval(coeff, [0, dimrecip]);
plot([0,dimrecip], fitLine)
title(['Direct Estimate of Entropy: ', num2str(fitLine(1)), ' bits/spk' ])
ylabel('H: bits/spk');
xlabel('Reciprocal dimension');



%%

H2pre = Entropy_logISI_H2(predbsSpkTimes/1000, 20);
H3pre = Entropy_logISI_HN(predbsSpkTimes/1000, 20, 3);



H1dur = entropyH1(durdbsISIhist);











% [~, psthSpkTimes,~] = peth(spkTimes, dbsTimes, [0:0.5:7.5]/1000);
% ISI = median(diff(dbsTimes));
% % specify output PSTH times for next section:
% PSTH = psthSpkTimes;
% % PSTH = psthSpkTimes(psthSpkTimes>=0 & psthSpkTimes<(7.5/1000));
% 
% % display results:
% % histogram(PSTH);
% % 
% % histogram(psthSpkTimes);
