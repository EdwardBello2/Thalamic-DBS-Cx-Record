function f = exploreNeuronRatesFsweep(dataTable,trialOrderStr)
% Show data for one neuron's rate changes across DBS conditions
% input is a matlab table for one neuron's data, assembled by the function
% "buildAnalysisTable_Rates.m"


% Author: Ed Bello
% Created: 2019/06/07



% check trialOrderStr input:
if ~strcmp(trialOrderStr, 'dbsFrequency') & ~strcmp(trialOrderStr, 'blockNum')
        error('acceptable string-inputs for trialOrderStr: dbsFrequency | blockNum');
        
end


% Initial variables
UnitID = dataTable.Unit_objectID{1};
freqs = [10, 20, 30, 50, 100, 130]; % dbs-frequencies for experiment
numFreqs = length(freqs);


binSeconds = 1; % number of seconds for each bin of spike counts


% Generate figure window first
f = figure;
f.Position = [2052 165 957 762];

numAx = 13;
for iAx = 1:numAx
    Ax(iAx) = subplot(4,6,iAx);

end

Ax(numAx) = subplot(4,6,13:24);



% initialize arrays for rate medians, quartiles
   MedsDBS = NaN(1, numFreqs);
quartUpDBS = NaN(1, numFreqs);
quartDnDBS = NaN(1, numFreqs);
   MedsPRE = NaN(1, numFreqs);
quartUpPRE = NaN(1, numFreqs);
quartDnPRE = NaN(1, numFreqs);




%% For one neuron, plot the Rates over time for each record in the upper
% axes

maxRate = 0;


% re-order the table based on trialOrderStr for the upper axes:
dataTable = sortrows(dataTable, trialOrderStr);

for iRow = 1:size(dataTable,1)
    % load the individual file with rate bins
    matfn = dataTable.ratebinFn{iRow};
    matpn = dataTable.ratebinPn{iRow};
    counts = load([matpn, '\', matfn]);

    ratesPRE = counts.bins.pre / binSeconds;
    ratesDBS = counts.bins.dbs / binSeconds;

    plot(Ax(iRow), ratesPRE, 'b');
    title(Ax(iRow), [num2str(dataTable.dbsFrequency(iRow)), 'Hz']);
    plot(Ax(iRow+6), ratesDBS, 'r');

%                 % Update descriptive stats for rate data
%                    MedsDBS(iRow) = prctile(ratesDBS, 50);
%                 quartUpDBS(iRow) = prctile(ratesDBS, 75) - MedsDBS(iRow);
%                 quartDnDBS(iRow) = MedsDBS(iRow) - prctile(ratesDBS, 25);
% 
%                    MedsPRE(iRow) = prctile(ratesPRE, 50);
%                 quartUpPRE(iRow) = prctile(ratesPRE, 75) - MedsPRE(iRow);
%                 quartDnPRE(iRow) = MedsPRE(iRow) - prctile(ratesPRE, 25);



    currentMaxRate = max(max(ratesPRE), max(ratesDBS));
    if currentMaxRate > maxRate
        maxRate = currentMaxRate;

    end

end



%% Now plot out the medians and quartiles in lower axes
        
% re-order the table based on trialOrderStr for the upper axes:
dataTable = sortrows(dataTable, 'dbsFrequency');


for iFreq = 1:numFreqs
    % first check to make sure the recording exists for each frequency
    isFrequency = (freqs(iFreq) == dataTable.dbsFrequency);
    
    if any(isFrequency)
        matfn = dataTable.ratebinFn{isFrequency};
        matpn = dataTable.ratebinPn{isFrequency};
        counts = load([matpn, '\', matfn]);
        
        ratesPRE = counts.bins.pre / binSeconds;
        ratesDBS = counts.bins.dbs / binSeconds;


        % Update descriptive stats for rate data
           MedsDBS(iFreq) = prctile(ratesDBS, 50);
        quartUpDBS(iFreq) = prctile(ratesDBS, 75) - MedsDBS(iFreq);
        quartDnDBS(iFreq) = MedsDBS(iFreq) - prctile(ratesDBS, 25);

           MedsPRE(iFreq) = prctile(ratesPRE, 50);
        quartUpPRE(iFreq) = prctile(ratesPRE, 75) - MedsPRE(iFreq);
        quartDnPRE(iFreq) = MedsPRE(iFreq) - prctile(ratesPRE, 25);
        
    end
    
end

        
        


% for iRow = 1:size(dataTable,1)
%     % load the individual file with rate bins
%     matfn = dataTable.ratebinFn{iRow};
%     matpn = dataTable.ratebinPn{iRow};
%     counts = load([matpn, '\', matfn]);
% 
%     ratesPRE = counts.bins.pre / binSeconds;
%     ratesDBS = counts.bins.dbs / binSeconds;
% 
%     
%     % Update descriptive stats for rate data
%        MedsDBS(iRow) = prctile(ratesDBS, 50);
%     quartUpDBS(iRow) = prctile(ratesDBS, 75) - MedsDBS(iRow);
%     quartDnDBS(iRow) = MedsDBS(iRow) - prctile(ratesDBS, 25);
% 
%        MedsPRE(iRow) = prctile(ratesPRE, 50);
%     quartUpPRE(iRow) = prctile(ratesPRE, 75) - MedsPRE(iRow);
%     quartDnPRE(iRow) = MedsPRE(iRow) - prctile(ratesPRE, 25);
% 
%     
%     currentMaxRate = max(max(ratesPRE), max(ratesDBS));
%     if currentMaxRate > maxRate
%         maxRate = currentMaxRate;
% 
%     end
% 
% end  

% Plot the medians and quartiles in the bottom axes
errorbar(Ax(13), freqs, MedsPRE, quartDnPRE, quartUpPRE, '-bo');
hold on;
errorbar(Ax(13), freqs, MedsDBS, quartDnDBS, quartUpDBS, '-ro');


% set all rate displays across all axes to the same scale
for iAx = 1:(length(Ax) - 1)
    Ax(iAx).YLim = [0, maxRate];

end


% Set the titles, labels etc for all axes
ylabel(Ax(1), 'PRE-dbs');
ylabel(Ax(7), 'DBS-on');

xlabel(Ax(end), 'DBS frequency (Hz)');
ylabel(Ax(end), 'Neural spk rate (Hz)')
title(Ax(end), [UnitID, ': Median and quartiles'])
Ax(end).XLim = [0, 140];



end