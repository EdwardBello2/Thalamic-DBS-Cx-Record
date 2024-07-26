


% read in nexfile
nexFile = readNexFile('D:\DataProcessing\170804session01\170804session01_SPKch12_Block07_inProgress.nex');


edgeLeft = 0/1000; % seconds
winWidth = 6/1000; % seconds

% test case
edgeLeft = 0; % seconds
winWidth = 2; % seconds


lSamp = round(edgeLeft*fs);
winSamp = round(winWidth*fs);


row = 8;
spkWaveforms = nexFile.waves{row,1}.waveforms';

spkTimes = nexFile.waves{row,1}.timestamps;
dbsTimes = nexFile.events{end,1}.timestamps;

% test case
dbsTimes = [10 20 30 40];
spkTimes = [50 60];

spkSamps = round(spkTimes*fs);

wavLen = size(spkWaveforms,2);
spkIdx = [1:length(spkSamps)]';

dbsSamps = round(dbsTimes*fs);

stim = 1;
sus = 1;
ISI = median(diff(dbsSamps));
maxStim = 4;
clear suspectSamps suspectWaveforms suspectIdx suspectTimes
done = false;
totspks = numel(spkTimes);
n = 1;
for st = 1:numel(dbsTimes)
  

    while spkTimes(n) < dbsTimes(st)
        n = n+1;
        if n>totspks
            done = true; 
        end
        if done == true, break; end

    end
    
    if done == true, break; end

    while spkTimes(n) > dbsTimes(st) && spkTimes(n) < dbsTimes(st)+2

            suspectTimes(sus) = spkTimes(n);
%             suspectWaveforms(sus,:) = spkWaveforms(n,:);
%             suspectIdx(sus) = spkIdx(n);

            sus = sus+1;
            n = n+1;
            if n>totspks, done = true; end
            if done == true, break; end
    end  

    if done == true, break; end

end
        
        
    