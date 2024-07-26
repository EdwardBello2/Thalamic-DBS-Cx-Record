


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
spkTimes = [2 9 21 22 31 41 60];

spkSamps = round(spkTimes*fs);

wavLen = size(spkWaveforms,2);
spkIdx = [1:length(spkSamps)]';

dbsSamps = round(dbsTimes*fs);

stim = 1;
sus = 1;
ISI = median(diff(dbsSamps));
maxStim = 4;
clear suspectSamps suspectWaveforms suspectIdx
for n = 1:numel(spkSamps)
    n
%     while spkSamps(n) < (dbsSamps(stim)-lSamp+winSamp)
%         stim = stim+1;
%     end
    if spkSamps(n) > (dbsSamps(stim)-lSamp)
        while spkSamps(n) > (dbsSamps(stim)+ISI-lSamp)
            stim = stim+1;
            if stim>maxStim
                break; 
            end
        end
        if spkSamps(n) > (dbsSamps(stim)-lSamp) && spkSamps(n) < (dbsSamps(stim)-lSamp+winSamp)

            suspectSamps(sus) = spkSamps(n)-dbsSamps(stim)-lSamp;
%             suspectWaveforms(sus,:) = spkWaveforms(n,:);
%             suspectIdx(sus) = spkIdx(n);

            sus = sus+1;

        else     

        end
        
        
    end
    
    
end
