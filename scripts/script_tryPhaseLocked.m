% script for assessing pulse % spike following for DBS, and estimate
% maximum following frequency, see if groups form



% CONSTANTS

% string input to affect the behavior of mergeTables_Master and
% filterTable_Master:
CURR_FUNC = 'script_displayPopRates_PREandDBS_barplot';

% PSTH parameters:
BW = 0.1 / 1000; % seconds
BMAX = 7.5 / 1000; % seconds
binEdges = 0:BW:BMAX;
timePSTH = binEdges(1:end-1) * 1000; % ms
timePSTHsec = binEdges(1:end-1); % s



%% LOAD all relevant tables and MERGE them

close all; clear


script_pipelineParams

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);

% isFiveNeurons = Tselect.neuCount == 5;
% Tselect = Tselect(isFiveNeurons,:);

% Detect those individual trials RECORDS in which neurons phase-locked 
% phase-locked must have a p-value less than ppar.pAlphaPSTH
% detection gets added as an extra column of information about each RECORD

% Add in phase-locked TF value for each row
nRows = size(Tselect, 1);
PhaseLocked = false(nRows, 1);

for iRec = 1:nRows
    if Tselect.pVal_Hpsth(iRec,1) <= ppar.pAlphaPSTH
        PhaseLocked(iRec,1) = true;
        
    end
    
end

Tselect = [Tselect, table(PhaseLocked)];



% DETECT which neurons had at least one phase-lock reaction to DBS 

% go thru each detected phase-lock, find its Unit_objectID, mark every row
% for that unit at PhaseLockNeuron = 1

idxPhsLck = find(Tselect.PhaseLocked);
nPhsLckRows = numel(idxPhsLck);
isPhsLckNeuron = false(size(Tselect, 1) , 1);

for i = 1:nPhsLckRows 
    iRow = idxPhsLck(i);
        
    % find the Unit_objectID belonging to this row
    currentNeuron = Tselect.Unit_objectID{iRow};
    
    % find all rows that pertain to the current neuron
    idxCurrNeuron = strcmp(currentNeuron, Tselect.Unit_objectID);
    
    % update column tracking all neurons showing at least one phase-lock
    isPhsLckNeuron = isPhsLckNeuron | idxCurrNeuron;
    
end

Tselect = [Tselect, table(isPhsLckNeuron)];


Tselect(~Tselect.isPhsLckNeuron,:) = [];

unique(Tselect.Unit_objectID)



%%

% Look at just one phase-locked cell's trial responses
isN = strcmp(Tselect.Unit_objectID, 'nr17082801-ch05a');
Tfinal = Tselect(isN,:);
Tfinal = sortrows(Tfinal, 'dbsFrequency', 'ascend');



%%





% initialize maximum phase-lock binEdges with something that will for sure
% get updated
maxEdgeRight = -Inf;
maxEdgeLeft = Inf;
maxYLim = 0;

figure; 

nRows = size(Tfinal, 1);
for iRow = 1:nRows

    % % for one NEX file:
    nexFile = readNexFile([ppar.projRootPath, '\', Tfinal.nexFileFolder{iRow}, ...
                          '\', Tfinal.nexFile{iRow}]);
    [spkTimes, stims] = parseNexFile(nexFile); % outputs timestamps for spiketimes and DBS times (seconds)

    % for DBS-on
    PSTHcount_dbs = psth(spkTimes, stims.DBS, binEdges, 'count');
    PSTHratenorm_dbs = (PSTHcount_dbs / BW) / numel(stims.DBS); % rate-normalized PSTH

    % FRmean = mean(PSTHratenorm_dbs);
    % FRthreshold = FRmean + (0.66 * (FRpeak - FRmean));
    % FRthreshold = FRmean_pre + (numSTD * FRmean_preSTD);
    FRthreshold = sum(PSTHcount_dbs) / (numel(stims.DBS) * (binEdges(end) - binEdges(1)));

    
    % Plot data into current axes pane within larger figure:
    ax(iRow) = subplot(2, 3, iRow);
    plot(timePSTH, PSTHratenorm_dbs);
    hold on; plot([timePSTH(1), timePSTH(end)], [FRthreshold, FRthreshold], '--r')
    ylabel('Spikes/sec');
    xlabel('ms');
    t(iRow) = title([num2str(Tfinal.dbsFrequency(iRow)), 'Hz']);

    if ax(iRow).YLim(2) > maxYLim, maxYLim = ax(iRow).YLim(2); end
    
    if Tfinal.pVal_Hpsth(iRow) <= 0.05
        isLocked = true;
        
    else
        isLocked = false;
        
    end
    
    
    % Update "maxEdgeLeft" and "maxEdgeRight" using only Phase-locked
    % trials:
    if isLocked
        % Among those trial responses that sig-phase-lock, choose the widest
        % sig-bin for spike-counting
        isAboveThresh = PSTHratenorm_dbs >= FRthreshold;

        [FRpeak, peakIdx] = max(PSTHratenorm_dbs);

        tPeakLatency = timePSTH(peakIdx);

        % find left and right borders above FRtrehshold
        startIdx = peakIdx;
        plot(timePSTH(startIdx), FRpeak, 'k*');


        % find right border
        currIdx = startIdx;
        maxIdx = numel(PSTHratenorm_dbs) - 1;

        % finding right-border I pad the PSTH with -1 at the end to ENSURE that
        % there is at least one zero-crossing at the very end:
        PSTHratenorm_dbs_Padded = [PSTHratenorm_dbs, -1];

        while ((PSTHratenorm_dbs_Padded(currIdx) - FRthreshold) > 0) && ...
              ((PSTHratenorm_dbs_Padded(currIdx+1) - FRthreshold) > 0)

            % If no zero-cross before we run out of samples, assign last sample
            % as the zero-crossing, otherwise increment current Index till we
            % find a zero-crossing:
            if currIdx >= maxIdx
                currIdx = maxIdx;
                break

            else
                currIdx = currIdx+1;

            end

        end
        binEdgeRight = currIdx + 1;
        plot(timePSTH(binEdgeRight), FRthreshold, 'r*');
        if binEdgeRight > maxEdgeRight, maxEdgeRight = binEdgeRight; end


        % find left border
        currIdx = startIdx;
        minIdx = 2;
        while ((PSTHratenorm_dbs(currIdx) - FRthreshold) > 0) && ...
              ((PSTHratenorm_dbs(currIdx-1) - FRthreshold) > 0)

            % If no zero-cross before we run out of samples, assign first sample
            % as the zero-crossing, otherwise increment current Index till we
            % find a zero-crossing:
            if currIdx <= minIdx
                currIdx = minIdx;
                break

            else
                currIdx = currIdx-1;

            end

        end
        binEdgeLeft = currIdx - 1;
        plot(timePSTH(binEdgeLeft), FRthreshold, 'r*');
        if binEdgeLeft < maxEdgeLeft, maxEdgeLeft = binEdgeLeft; end
        
    end % END "if isLocked"
    
    
    

%     
%     isAboveThresh = PSTHratenorm_dbs >= FRthreshold;
%     tPhaseLock = t(isAboveThresh);
%     tPeakLeft = min(tPhaseLock);
%     tPeakRight = max(tPhaseLock);
%     peakJitter = tPeakRight - tPeakLeft





% disp([num2str(Tfinal.dbsFrequency(iRow)), ' Hz'])
% [eEPF, pfs, pfsb] = calcEEPF(nexFile)
% title([num2str(Tfinal.dbsFrequency(iRow)), ' Hz, eEPF: ', num2str(eEPF)])

end

for i = 1:numel(ax), ax(i).YLim(2) = maxYLim; end


% Get final phase-lock time-bin to apply to each PSTH for this cell
tPhs = [timePSTHsec(maxEdgeLeft), timePSTHsec(maxEdgeRight)];




%% Find: 1) % pulse entrainment, and 2) max following frequency based on pulses

nRows = size(Tfinal, 1);
for iRow = 1:nRows
    % add dashed lines to show updated phase-locked window at each PSTH
    plot(ax(iRow), [timePSTH(maxEdgeLeft), timePSTH(maxEdgeLeft)], ...
         [ax(iRow).YLim(1), ax(iRow).YLim(2)], '--');
    plot(ax(iRow), [timePSTH(maxEdgeRight), timePSTH(maxEdgeRight)], ...
         [ax(iRow).YLim(1), ax(iRow).YLim(2)], '--'); 
     
     
    % for one recording:
    nexFile = readNexFile([ppar.projRootPath, '\', Tfinal.nexFileFolder{iRow}, ...
                          '\', Tfinal.nexFile{iRow}]);
    [spkTimes, stims] = parseNexFile(nexFile); % outputs timestamps for spiketimes and DBS times (seconds)


    % 1) pulse entrainment percent:
    nLocked = psth(spkTimes, stims.DBS, [tPhs(1), tPhs(2)], 'count');
    percLocked = 100 * nLocked / numel(stims.DBS);
    t(iRow).String = [t(iRow).String, ', Locked: ', num2str(percLocked), ' %'];
    
    
    % 2) phase-locked following-frequency:
    followFreq = Tfinal.dbsFrequency(iRow) * (percLocked / 100);
    t(iRow).String = [t(iRow).String, ', follow-Freq: ', num2str(followFreq), ' Hz'];
    
    
    
    
    
    
end




% SUB-FUNCTIONS
% 
% function findPSTHcrossThresh(PSTHratenorm_dbs, FRthreshold)
% % finds the leftmost and rightmost bin-edges in a PSTH that cross the
% % specified threshold line
% 
% % Among those trial responses that sig-phase-lock, choose the widest
% % sig-bin for spike-counting
% % isAboveThresh = PSTHratenorm_dbs >= FRthreshold;
% 
% [FRpeak, peakIdx] = max(PSTHratenorm_dbs);
% 
% tPeakLatency = timePSTH(peakIdx);
% 
% % find left and right borders above FRtrehshold
% startIdx = peakIdx;
% plot(timePSTH(startIdx), FRpeak, 'k*');
% 
% 
% % find right border
% currIdx = startIdx;
% maxIdx = numel(PSTHratenorm_dbs) - 1;
% 
% % finding right-border I pad the PSTH with -1 at the end to ENSURE that
% % there is at least one zero-crossing at the very end:
% PSTHratenorm_dbs_Padded = [PSTHratenorm_dbs, -1];
% 
% while ((PSTHratenorm_dbs_Padded(currIdx) - FRthreshold) > 0) && ...
%       ((PSTHratenorm_dbs_Padded(currIdx+1) - FRthreshold) > 0)
% 
%     % If no zero-cross before we run out of samples, assign last sample
%     % as the zero-crossing, otherwise increment current Index till we
%     % find a zero-crossing:
%     if currIdx >= maxIdx
%         currIdx = maxIdx;
%         break
% 
%     else
%         currIdx = currIdx+1;
% 
%     end
% 
% end
% binEdgeRight = currIdx + 1;
% plot(timePSTH(binEdgeRight), FRthreshold, 'r*');
% if binEdgeRight > maxEdgeRight, maxEdgeRight = binEdgeRight; end
% 
% 
% % find left border
% currIdx = startIdx;
% minIdx = 2;
% while ((PSTHratenorm_dbs(currIdx) - FRthreshold) > 0) && ...
%       ((PSTHratenorm_dbs(currIdx-1) - FRthreshold) > 0)
% 
%     % If no zero-cross before we run out of samples, assign first sample
%     % as the zero-crossing, otherwise increment current Index till we
%     % find a zero-crossing:
%     if currIdx <= minIdx
%         currIdx = minIdx;
%         break
% 
%     else
%         currIdx = currIdx-1;
% 
%     end
% 
% end
% binEdgeLeft = currIdx - 1;
% plot(timePSTH(binEdgeLeft), FRthreshold, 'r*');
% if binEdgeLeft < maxEdgeLeft, maxEdgeLeft = binEdgeLeft; end
% 
% 
% 
% 
% end




% count all phase-lock spikes for each trial and get %, rate