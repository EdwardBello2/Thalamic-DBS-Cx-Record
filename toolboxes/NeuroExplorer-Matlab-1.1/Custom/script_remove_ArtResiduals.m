% script for easily excluding waveforms that are obviously the dbs
% artifact, by virtue of crossing the threshold line right at the time of
% DBS pulse, or VERY near it


addpath('C:\Users\bello043\OneDrive\Tools\NEXtools')

% read in nexfile
nexFile = readNexFile('D:\DataProcessing\170804session01\170804session01_SPKch12_Block07_inProgress.nex');




% plot a pile of each units waveforms that occurr during the DBS pulse 
ldbs = nexFile.events{end,1}.name;
tdbs = nexFile.events{end,1}.timestamps;

unit = 1;

spkName = nexFile.waves{unit,1}.name;
spkWavs = nexFile.waves{unit,1}.waveforms;
spkTimes = nexFile.waves{unit,1}.timestamps;
% figure; plot(spkWavs)

window = 3/1000;
winSamps = round(window*fs);

lWav = size(spkWavs,1);
midpt = winSamps+1;
t = (-winSamps)/fs:1/fs:winSamps/fs;
t=t*1000; % to put the timescale in ms

f = figure; ax = axes; 
ax.XLim = [-winSamps winSamps];


for j = 1:numel(tdbs)
% for the nth ref event (i.e. a DBS pulse) see if any spikes occur near it
ref = tdbs(j);
refSamp = round(ref*fs);

suspectMask = false(length(spkTimes),1);


for i = 1:numel(spkTimes)
    if spkTimes(i) > ref-window & spkTimes(i) < ref+window
        suspectMask(i) = true;
    end
    
end

if any(suspectMask)
    suspectWavs = spkWavs(:,suspectMask);
    suspectTimes = spkTimes(suspectMask);
    suspectSamps = round(suspectTimes*fs);
    
    adjWavs = zeros((2*winSamps+1),size(suspectWavs,2));
    for k = 1:size(adjWavs,2)
        spkpt = suspectSamps-refSamp;
        wavpt = midpt+spkpt-13;
        % put current spkWav into adjWavs
        if wavpt<1
            offset = 1-wavpt;
            adjWavs(1:(1+lWav-1-offset),k) = suspectWavs(1:end-offset,k);

        elseif wavpt+lWav>2*winSamps+1
            offset = (wavpt+lWav)-(2*winSamps+1);
            adjWavs(wavpt:(wavpt+lWav-1-offset),k) = suspectWavs(1:end-offset,k);

        else
            adjWavs(wavpt:(wavpt+lWav-1),k) = suspectWavs(:,k);

        end


    plot(t,adjWavs); hold on
    end
    
    
end



end
plot([0 0],ax.YLim)
grid on; grid minor
xlabel('Time (milliseconds)')
title([spkName ' waves aligned to ' ldbs],'Interpreter','none')


% for that unit, decide if waveforms can be excluded based on what time
% window the valley falls within


% organize the offending waveforms it into a new "DBSartifact" unit, but
% keep innocent waveforms in their original unit

