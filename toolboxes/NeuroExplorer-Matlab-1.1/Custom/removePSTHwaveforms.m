function removePSTHwaveforms(nexFile,edgeLeft,winWidth,prethreshSamps,unit)

% read in nexfile
% nexFile = readNexFile('D:\DataProcessing\170804session01\170804session01_SPKch12_Block07_inProgress.nex');

fs = nexFile.waves{1,1}.WFrequency;
% edgeLeft = 0/1000; % seconds
% winWidth = 6/1000; % seconds

% % test case
% edgeLeft = 0; % seconds
% winWidth = 2; % seconds


lSamp = round(edgeLeft*fs);
winSamp = round(winWidth*fs);


% row = 1;
spkWaveforms = nexFile.waves{unit,1}.waveforms';

spkTimes = nexFile.waves{unit,1}.timestamps;
dbsTimes = nexFile.events{end,1}.timestamps;

% % test case
% dbsTimes = [10 20 30 40];
% spkTimes = [50 60];

spkTSamps = round(spkTimes*fs);

wavLen = size(spkWaveforms,2);
spkIdx = [1:length(spkTSamps)]';

dbsTSamps = round(dbsTimes*fs);

ISI = median(diff(dbsTSamps));
% clear suspectSamps suspectWaveforms suspectIdx suspectTimes
[suspectIdx,suspectTSamps,suspectWaveforms] = gatherSuspects(spkIdx,spkTSamps,spkWaveforms,dbsTSamps,lSamp,winSamp);




[PSTHwaveforms] = gatherPSTH(suspectTSamps,suspectWaveforms,prethreshSamps,winSamp);

figure; plot(PSTHwaveforms'); 
spkName = nexFile.waves{unit,1}.name;
ldbs = nexFile.events{end,1}.name;

title([spkName ' waves aligned to ' ldbs],'Interpreter','none')
%--------------------
% 
% prethreshSamps = 13;
% totPSTH = size(suspectWaveforms,1);
% wavSampL = size(suspectWaveforms,2);
% 
% PSTHwaveforms = zeros(totPSTH,winSamp);
% susPSTHIdx = zeros(totPSTH,wavSampL);
% 
% 
% for p = 1:totPSTH
%     
%    
%     % fill the susPSTHIdx with it's marker of where it fits in the PSTH
%     for i = 1:wavSampL
%         susPSTHIdx(p,i) = suspectTSamps(p)-prethreshSamps+i;
%     end
%    
%    
%     % fill each PSTH window with suspected waveforms point by point
%     s = 1;
%     for w = 1:winSamp
%         if w == susPSTHIdx(p,s)
%             
%             PSTHwaveforms(p,w) = suspectWaveforms(p,s);
%             s = s+1;
%             if s>wavSampL, break; end
%             
%         end
%     end
%     
% end



%-----------------------





        
        
        

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
        
        
    