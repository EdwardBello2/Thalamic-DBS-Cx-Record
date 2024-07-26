%%
clear all, clc
% find(diff(F.CC1_n1)<0)
%% Script to calculate the entropy of a neuron from its PSTH (Letters method)
% First have the user select all the data to load in (_proc files)
[FileName,PathName] = uigetfile('*.mat','Select the MATLAB code file', 'MultiSelect', 'on');


Settings.FileName = FileName; 
Settings.PathName = PathName; 
Settings.binSize = 0.5; % ms
Settings.cropTime_front = 0.5;% (ms), this is the time to crop out at the beginning of the PSTH (to remove blanking artifacts)
Settings.cropTime_end = 1; % (ms), this is the time to crop out at the end of the PSTH (to remove blanking artifacts)
Settings.psthLength = 10; % (ms), this is the total length of the PSTH, or the time between 2 successive stim pulses
Settings.colNum = 6; % number of columns per stim amplitude
Settings.bootNum = 10000; % number of bootstraps trials
Settings.stimNum = 3; % number of stimulations for each recording to compute. (some recordings might have 5, so will only use the 3 with highest stim amplitudes)


Settings.dispData = zeros(numel(Settings.FileName), Settings.stimNum*Settings.colNum); % the format is: amp1, H1-pr, H1-on, H2-pr,H2-on, H3-pr, H3-on, amp2, ...
Settings.bins = 0:Settings.binSize:Settings.psthLength; % The edges of the bins, from left to right


% Load each file individually and compute all the data
h = waitbar(0,'Please wait...');
for rr = 1: numel(Settings.FileName)
    % Load the _proc.mat file 
    load([Settings.PathName, Settings.FileName{rr}]); % Now 'F' is loaded 
    spks = F.CC1_n1; % get all spike times, change to 'F.CC1_n2' if there's more than 1 cell 
    [~,stimOrder] = sort([F.stim(:).splineMeanArtMax], 'descend');
   
    % Loop through each stim session and calculate the pre-stim H and during-stim H
    for i = 1:Settings.stimNum
        % Convert the spike times into PSTH format, i.e. into time since the last stim(virtual stim) pulse
        pr_spks = spks(spks>F.stim(stimOrder(i)).pr(1) & spks<F.stim(stimOrder(i)).pr(end));
        on_spks = spks(spks>F.stim(stimOrder(i)).on(1) & spks<F.stim(stimOrder(i)).on(end));
        po_spks = spks(spks>F.stim(stimOrder(i)).po(1) & spks<F.stim(stimOrder(i)).po(end));
       
        % This calculates the distance between all spikes and the immediately preceding stimualtion pulse
        pr_st_sp_int = pr_spks*ones(1,length(F.stim(stimOrder(i)).pr))-ones(length(pr_spks),1)*F.stim(stimOrder(i)).pr'; % these 3 lines calculates the distance between every spike and every pulse (spike time - pulse time)
        on_st_sp_int = on_spks*ones(1,length(F.stim(stimOrder(i)).on))-ones(length(on_spks),1)*F.stim(stimOrder(i)).on';
        po_st_sp_int = po_spks*ones(1,length(F.stim(stimOrder(i)).po))-ones(length(po_spks),1)*F.stim(stimOrder(i)).po';
        
        % these 3 lines substitutes negative values (corresponding to pulses following the spike) with very large positive numbers so to select the minimum positive distance for every spike (the distance between the spike and the immediately preceding pulse)
       
        % Here are now the times since the previous stim pulse for all 3 periods
        pr_st_sp_int = min(pr_st_sp_int+(pr_st_sp_int<0)*10^6,[],2); % pr period
        on_st_sp_int = min(on_st_sp_int+(on_st_sp_int<0)*10^6,[],2); % on period
        po_st_sp_int = min(po_st_sp_int+(po_st_sp_int<0)*10^6,[],2); % po period
        
        % Convert these times to ms (from s)
        pr_st_sp_int = pr_st_sp_int*1000;
        on_st_sp_int = on_st_sp_int*1000;
        po_st_sp_int = po_st_sp_int*1000;
        
        % Crop the beginning and end of the PSTH to remove blanking artifacts
        pr_st_sp_int(find(pr_st_sp_int < Settings.cropTime_front | pr_st_sp_int > (Settings.psthLength - Settings.cropTime_end))) = [];
        on_st_sp_int(find(on_st_sp_int < Settings.cropTime_front | on_st_sp_int > (Settings.psthLength - Settings.cropTime_end))) = [];
        po_st_sp_int(find(po_st_sp_int < Settings.cropTime_front | po_st_sp_int > (Settings.psthLength - Settings.cropTime_end))) = [];
        
        % Bin the PSTH data and calculate the H using the relative frequency of the 'letters' in the PSTH
        [H_pr, binCounts_pr] = Entropy_psthLetter_H1(pr_st_sp_int, Settings.bins);
        [H_on, binCounts_on] = Entropy_psthLetter_H1(on_st_sp_int, Settings.bins);
        H_diff = H_pr - H_on; % absolute diff
        H_percentDiff = (H_pr - H_on)/H_pr; % percent diff
        
        % Calculate the boostrap H distribution and p-value
        bootDist = zeros(1, Settings.bootNum); % This is the bootstrapped H distribution from sub-sampled spikes from the pr period
        for j = 1: Settings.bootNum
            Sample = datasample(pr_st_sp_int, numel(on_st_sp_int), 'Replace', true); % Sample from pr-period spikes (w/ replacement), number is equal to number of stim-on spikes
            [H_Sample, ~] = Entropy_psthLetter_H1(Sample, Settings.bins); % Calculate the entropy
            bootDist(j) = H_Sample;
            
        end
        
        % Calculate the p-value. Usually the modulated cells have lower H in the on-stim period compared to the pr-stim period, so look for the left tail of the p-value.
        pValue = sum(bootDist < H_on)/(Settings.bootNum); % p-value
        Settings.dispData(rr,(i-1)*Settings.colNum+1)= [F.stim(stimOrder(i)).splineMeanArtMax]; % Store the artifact amplitudes
        Settings.dispData(rr,(i-1)*Settings.colNum+2: (i-1)*Settings.colNum+Settings.colNum) = [H_pr, H_on, H_diff, H_percentDiff, pValue]; % Store the other info
        
    end % END FOR-LOOP for each stim-trial
    
    
    % Clear variables to save memory
    clearvars -except Settings rr h
    
    waitbar(rr/numel(Settings.FileName))
    
    
    
end % END FOR-LOOP for individual files

close(h)

%% Percent Diff vs amp 

plot(Hdata.amp, Hdata.percentDiff*100, 'k.', 'markersize',10); hold on;
ind = find(Hdata.pValue < 0.05); 
plot(Hdata.amp(ind), Hdata.percentDiff(ind)*100, 'r.', 'markersize',10)
legend('p>0.05', 'p<0.05');
title('% Change in H using PSTH letter method, binsize = 0.5ms');
xlabel('Stim Amp (uV)')
ylabel('Percent Change (%)')
%% pValue vs amp
plot(Hdata.amp, Hdata.pValue, 'k.', 'markersize',10); hold on; 
ind = find(Hdata.pValue < 0.05); 
plot(Hdata.amp(ind), Hdata.pValue(ind), 'r.', 'markersize',10)
hline = refline([0, 0.05]);
hline.Color = 'r';
title('P-values using PSTH letters method, binsize = 0.5ms');
xlabel('Stim Amp (uV)')
ylabel('p-value')

%% percent change vs. p-value 
plot(Hdata.percentDiff*100, Hdata.pValue, 'k.', 'markersize',10); hold on; 
ind = find(Hdata.pValue < 0.05); 
plot(Hdata.percentDiff(ind)*100, Hdata.pValue(ind), 'r.', 'markersize',10)
hline = refline([0, 0.05]);
hline.Color = 'r';
title('% change vs. p-value');
xlabel('Percent Change (%)')
ylabel('p-value')