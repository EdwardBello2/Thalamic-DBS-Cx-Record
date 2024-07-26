function nexFileNEW = separateArtifactWaveforms(nexFile,edgeLeft,winWidth,prethreshTime)
% Plot a PSTH of of all spike waveforms for each unit in the nexFile.
% Waveforms are shown aligned to event times of interest (default event is
% the "last" event in the list of events in your nexFile).
%
% edgeLeft: time in seconds to the left of your event to start the PSTH
% winWidth: total time window (seconds) of your PSTH
% prethreshTime: the stretch of time in the spike waveform before spike detection point  
% 
% read in nexfile
% nexFile = readNexFile('D:\DataProcessing\170804session01\170804session01_SPKch12_Block07_inProgress.nex');
nexFileNEW = nexFile;

global PSTHs

if ~isfield(nexFile,'waves')

    error('function requires the .NEX file to have the "waves" field!')
    
end


fs = nexFile.waves{1,1}.WFrequency;


prethreshSamps = round(prethreshTime*fs);
lSamp = round(edgeLeft*fs);
winSamp = round(winWidth*fs);


% DBS related variables
dbsTimes = nexFile.events{end,1}.timestamps;
dbsTSamps = round(dbsTimes*fs);
ISI = median(diff(dbsTSamps));

%% Sifting process for each individual unit


numUnits = size(nexFile.waves,1);
for unit = 1:numUnits
%     unit = 1;
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

    f = figure; ax = axes;
    plot(t,PSTHwaveforms'); hold on;
    
    % Add in a line showing the time of the EVENT in the PSTH (i.e. DBS
    % pulse)
    plot([0 0],(ax.YLim)*1.2,'r','LineWidth',1);
    grid on; grid minor
    xlabel('Time (milliseconds)')

    spkName = nexFile.waves{unit,1}.name;
    ldbs = nexFile.events{end,1}.name;

    title([spkName ' waves aligned to ' ldbs],'Interpreter','none')
    
    
    % Gather all SUSPECT waveforms whose timestamp (in samples) falls
    % within the tickmarks input by the user -- declare them GUILTY
    [ticks,~] = ginput;
    
    % If bounds were specified, 1) create new unit with guilty spks and 2)
    % remove guilty spikes from the original unit of the nexfile
    if ~isempty(ticks)
        [guiltyIdx,guiltyTSamps,guiltyWaveforms] = gatherGuilty(suspectIdx,suspectTSamps,suspectWaveforms,prethreshSamps,lSamp,winSamp,fs,ticks);   
    
        
               
        
        guiltyTimes = spkTimes(guiltyIdx);
        guiltyWaveforms = spkWaveforms(guiltyIdx,:)';
        
        % create new unit in nexFile for guilty spikes
        nexFileNEW = addguilty2newUnit(nexFileNEW,unit,guiltyTimes,guiltyWaveforms)
        
        % Remove all guilty spikes from original unit in nexFile
        spkTimes(guiltyIdx) = [];
        spkWaveforms(guiltyIdx,:) = [];
        
        nexFileNEW.waves{unit,1}.timestamps = spkTimes;
        nexFileNEW.waves{unit,1}.waveforms = spkWaveforms';
        
        if isfield(nexFile,'neurons')  
            
            nexFileNEW.neurons{unit,1}.timestamps = spkTimes;
            
        end
        
    
    end

    close(f)

end

%% Cleanup of empty units

% check if there are any units that ended up empty
nospikes = false(size(nexFileNEW.waves,1),1);
for j = 1:size(nexFileNEW.waves,1)
    times = nexFileNEW.waves{j, 1}.timestamps  ;
    if isempty(times)
        nospikes(j) = true
    end
    
    
end

% remove any such empty units
nexFileNEW.waves(nospikes) = [];
nexFileNEW.neurons(nospikes) = [];

% rename units to account for any removals
nexFileNEW = redoLabels(nexFileNEW)


% display all units PSTHs at the end for reference while spike sorting
numUnits = size(nexFileNEW.waves,1)

for k = 1:numUnits
    figure;
    
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
        
function [guiltyIdx,guiltyTSamps,guiltyWaveforms] = gatherGuilty(suspectIdx,suspectTSamps,suspectWaveforms,prethreshSamps,lSamp,winSamp,fs,ticks)
% from "suspect spikes", label those that feel within the ginput boundaries
% as "guilty" (of being DBS artifacts); the rest are "innocent", but are
% not output in this function. 

t = (-lSamp)/fs:1/fs:(winSamp-lSamp-1)/fs;
t=t*1000; % to put the timescale in ms

ticks = ticks(1:2); % only the first two points are used.
shiftedTicks = ticks - t(1);

bounds = round((shiftedTicks/1000)*fs);

guil = find(suspectTSamps>bounds(1) & suspectTSamps<bounds(2));
inno = find(suspectTSamps<=bounds(1) | suspectTSamps>=bounds(2));


% define guilty spikes and waveforms + plot
guiltyIdx = suspectIdx(guil);
guiltyTSamps = suspectTSamps(guil);
guiltyWaveforms = suspectWaveforms(guil,:);

[PSTHguilty] = gatherPSTH(guiltyTSamps,guiltyWaveforms,prethreshSamps,winSamp);

% t = (-lSamp)/fs:1/fs:(winSamp-lSamp-1)/fs;
% t=t*1000; % to put the timescale in ms


f1 = figure; 
f1.Position = [386 359 1251 556];
ax1 = subplot(1,2,1);


if isempty(guil)
    
    PSTHguilty = zeros(1,winSamp);

    plot(t,PSTHguilty'); hold on;
    plot([0 0],(ax1.YLim)*1.2,'r','LineWidth',1);
    grid on; grid minor
    xlabel('Time (milliseconds)')
    title('GUILTY SPIKES')
    
else
%     f1 = figure; ax = axes;

    plot(t,PSTHguilty'); hold on;
    plot([0 0],(ax1.YLim)*1.2,'r','LineWidth',1);
    grid on; grid minor
    xlabel('Time (milliseconds)')
    title('GUILTY SPIKES')
end

ax2 = subplot(1,2,2);

% define guilty spikes and waveforms + plot
if isempty(inno)
    PSTHinnocent = zeros(1,winSamp);
    
%     f2 = figure; ax = axes;
    plot(t,PSTHinnocent'); hold on;
    plot([0 0],(ax2.YLim)*1.2,'r','LineWidth',1);
    grid on; grid minor
    xlabel('Time (milliseconds)')
    title('INNOCENT SPIKES')
else
    
    innocentTSamps = suspectTSamps(inno);
    innocentWaveforms = suspectWaveforms(inno,:);

    [PSTHinnocent] = gatherPSTH(innocentTSamps,innocentWaveforms,prethreshSamps,winSamp);

%     f2 = figure; ax = axes;
    plot(t,PSTHinnocent'); hold on;
    plot([0 0],(ax2.YLim)*1.2,'r','LineWidth',1);
    grid on; grid minor
    xlabel('Time (milliseconds)')
    title('INNOCENT SPIKES')

end

input('Press ENTER when you are done with this unit')

try
close(f1)
% close(f2)
catch
    
end

end
    
function nexFileUPDATED = addguilty2newUnit(nexFile,currentUnit,guiltyTimes,guiltyWaveforms)
% create new unit in nexFile for guilty spikes

nexFileUPDATED = nexFile;

exampleWave = nexFile.waves{currentUnit,1}

exampleName = exampleWave.name;

exampleWave.name = [];
exampleWave.unitNumber = [];
exampleWave.timestamps = [];
exampleWave.waveforms = [];

newUnit = size(nexFile.waves,1)+1;

nexFileUPDATED.waves{newUnit,1} = exampleWave;

nexFileUPDATED.waves{newUnit,1}.name = generateUnitNameWF(exampleName,newUnit);
nexFileUPDATED.waves{newUnit,1}.unitNumber = newUnit;
nexFileUPDATED.waves{newUnit,1}.timestamps = guiltyTimes;
nexFileUPDATED.waves{newUnit,1}.waveforms = guiltyWaveforms;

% if the "neurons" field exists, do the same as the above for the new unit
if isfield(nexFile,'neurons')
    
    exampleNeuron = nexFile.neurons{currentUnit,1};
    exampleNeuron.name = [];
    exampleNeuron.unitNumber = [];
    exampleNeuron.timestamps = [];
    
    nexFileUPDATED.neurons{newUnit,1} = exampleNeuron;

    
    nexFileUPDATED.neurons{newUnit,1}.name = generateUnitNameNeu(exampleName,newUnit);
    nexFileUPDATED.neurons{newUnit,1}.unitNumber = newUnit;
    nexFileUPDATED.neurons{newUnit,1}.timestamps = guiltyTimes;
end


end

function unitName = generateUnitNameWF(exampleName,newUnit)

base = exampleName;
base(end-3:end) = [];

albet = 'abcdefghijklmnopqrstuvwxyz';

unitName = [base albet(newUnit) '_wf' ];



end

function unitName = generateUnitNameNeu(exampleName,newUnit)

base = exampleName;
base(end) = [];

albet = 'abcdefghijklmnopqrstuvwxyz';

unitName = [base albet(newUnit) ];



end

function nexFileNEW = redoLabels(nexFile)

nexFileNEW = nexFile;

numUnits = size(nexFile.waves,1);
neuronName = nexFile.neurons{1,1}.name;
wavName = nexFile.waves{1,1}.name;
start = 1;
if strcmp(neuronName(end),'U')
    start = 2;
end
alphbet = 1;
for j = start:numUnits
    
    nexFileNEW.neurons{j,1}.name = generateUnitNameNeu(neuronName,alphbet)
    nexFileNEW.waves{j,1}.name = generateUnitNameWF(wavName,alphbet)
    alphbet = alphbet+1;
end


end