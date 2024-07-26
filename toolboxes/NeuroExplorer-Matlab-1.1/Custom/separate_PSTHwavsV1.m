function separate_PSTHwavs(nexFile)
% go thru each sorted unit and separate out any spike waveforms that do not
% belong. Waveforms are shown aligned to a reference event.


% Default Settings
suspectPSTHwindow = 3; %milliseconds




window = suspectPSTHwindow/1000; % seconds

%% Code


dbsLabel = nexFile.events{end,1}.name;
dbsTimes = nexFile.events{end,1}.timestamps;

numUnits = size(nexFile.waves,1);

% plot a pile of each units waveforms that occurr during the DBS pulse 

for u = 1:numUnits

fs = nexFile.waves{u,1}.WFrequency;
spkName = nexFile.waves{u,1}.name;
spkWavs = nexFile.waves{u,1}.waveforms;
spkTimes = nexFile.waves{u,1}.timestamps;
% figure; plot(spkWavs)

winSamps = round(window*fs);
lWav = size(spkWavs,1);
midpt = winSamps+1;
t = (-winSamps)/fs:1/fs:winSamps/fs;
t=t*1000; % to put the timescale in ms

figure; ax = axes; 
ax.XLim = [-winSamps winSamps];


for j = 1:numel(dbsTimes)
% for the nth ref event (i.e. a DBS pulse) see if any spikes occur near it
ref = dbsTimes(j);
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
    
%     try
    [adjWavs] = adjustWaveforms(refSamp,midpt,lWav,winSamps,suspectWavs,suspectSamps);
    plot(t,adjWavs); hold on
%     catch
%         error(['problem in loop ' num2str(j)])
%     end


end %any(suspectMask)



end %j = 1:numel(dbsTimes)
plot([0 0],ax.YLim)
grid on; grid minor
xlabel('Time (milliseconds)')
title([spkName ' waves aligned to ' dbsLabel],'Interpreter','none')

[ticks,~] = ginput;
tickSamps = round((ticks/1000)*fs);

end %unit
end

%% Subfunctions

function [adjWavs] = adjustWaveforms(refSamp,midpt,lWav,winSamps,suspectWavs,suspectSamps)


adjWavs = zeros((2*winSamps+1),size(suspectWavs,2));
    for k = 1:size(adjWavs,2)
        spkpt = suspectSamps(k)-refSamp;
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


    end
end



% for that unit, decide if waveforms can be excluded based on what time
% window the valley falls within


% organize the offending waveforms it into a new "DBSartifact" unit, but
% keep innocent waveforms in their original unit

