function displayPSTHwaveforms(nexFile,edgeLeft,winWidth,prethreshTime)
% Plot a PSTH of of all spike waveforms for each unit in the nexFile.
% Waveforms are shown aligned to event times of interest (default event is
% the "last" event in the list of events in your nexFile).
%
% edgeLeft: time in seconds to the left of your event to start the PSTH
% winWidth: total time window (seconds) of your PSTH
% prethreshSamps: number of samples before spike detection 
% 
% read in nexfile
% nexFile = readNexFile('D:\DataProcessing\170804session01\170804session01_SPKch12_Block07_inProgress.nex');

fs = nexFile.waves{1,1}.WFrequency;


prethreshSamps = round(prethreshTime*fs);
lSamp = round(edgeLeft*fs);
winSamp = round(winWidth*fs);


% DBS related variables
dbsTimes = nexFile.events{end,1}.timestamps;
dbsTSamps = round(dbsTimes*fs);
ISI = median(diff(dbsTSamps));

%% Plot all units PSTHs in this nexFile


numUnits = size(nexFile.waves,1);
for unit = 1:numUnits
    
    % unit-specific variables
    spkWaveforms = nexFile.waves{unit,1}.waveforms';
    spkTimes = nexFile.waves{unit,1}.timestamps;
    spkTSamps = round(spkTimes*fs);

    wavLen = size(spkWaveforms,2);
    spkIdx = [1:length(spkTSamps)]';



    % gather all SUSPECTS -- the waveforms that fall within the user-specified PSTH window 
    [suspectIdx,suspectTSamps,suspectWaveforms] = gatherSuspects(spkIdx,spkTSamps,spkWaveforms,dbsTSamps,lSamp,winSamp);


    % Arrange each SUSPECT waveform into place within the PSTH 
    [PSTHwaveforms] = gatherPSTH(suspectTSamps,suspectWaveforms,prethreshSamps,winSamp);

    eventLoc = lSamp + 1;
    t = (-lSamp)/fs:1/fs:(winSamp-lSamp-1)/fs;
    t=t*1000; % to put the timescale in ms

    figure; ax = axes;
    plot(t,PSTHwaveforms'); hold on;
    
    % Add in a line showing the time of the EVENT in the PSTH (i.e. DBS
    % pulse)
    plot([0 0],(ax.YLim)*1.2,'r','LineWidth',1);
    grid on; grid minor
    xlabel('Time (milliseconds)')

    spkName = nexFile.waves{unit,1}.name;
    ldbs = nexFile.events{end,1}.name;

    title([spkName ' waves aligned to ' ldbs],'Interpreter','none')


end






        
        
        

end
%% SUBFUNCTIONS

function [suspectIdx,suspectTSamps,suspectWaveforms] = gatherSuspects(spkIdx,spkTSamps,spkWaveforms,dbsTSamps,lSamp,winSamp)
    
done = false;
totspks = numel(spkTSamps);
n = 1;
sus = 1;
for st = 1:numel(dbsTSamps)
  

    while spkTSamps(n) < (dbsTSamps(st)-lSamp)
        n = n+1;
        if n>totspks
            done = true; 
        end
        if done == true, break; end

    end
    
    if done == true, break; end
    
    while spkTSamps(n) >= (dbsTSamps(st)-lSamp) && spkTSamps(n) < (dbsTSamps(st)-lSamp+winSamp)

            suspectTSamps(sus) = spkTSamps(n) - (dbsTSamps(st)-lSamp);
            suspectWaveforms(sus,:) = spkWaveforms(n,:);
            suspectIdx(sus) = spkIdx(n);

            sus = sus+1;
            n = n+1;
            if n>totspks, done = true; end
            if done == true, break; end
    end  

    if done == true, break; end

end

end

function [PSTHwaveforms] = gatherPSTH(suspectTSamps,suspectWaveforms,prethreshSamps,winSamp)


% prethreshSamps = 13;
totPSTH = size(suspectWaveforms,1);
wavSampL = size(suspectWaveforms,2);

PSTHwaveforms = zeros(totPSTH,winSamp);
susPSTHIdx = zeros(totPSTH,wavSampL);


for p = 1:totPSTH
    
   
    % fill the susPSTHIdx with it's marker of where it fits in the PSTH
    for i = 1:wavSampL
        susPSTHIdx(p,i) = suspectTSamps(p)-prethreshSamps+i;
    end
   
   
    % fill each PSTH window with suspected waveforms point by point
    s = 1;
    for w = 1:winSamp
        if w == susPSTHIdx(p,s)
            
            PSTHwaveforms(p,w) = suspectWaveforms(p,s);
            s = s+1;
            if s>wavSampL, break; end
            
        end
    end
    
end


end
        
        
    