function deltaH_barPlot(ResultsTable)
% create two boxplot figures with boxplots for each DBS condition, one
% figure for each sweep-type

% Find out which rows are Contact-sweep and Freq-sweep
sweepType = ResultsTable.sweepType(:);
isCswp = strcmp(sweepType, 'Contact');
isFswp = strcmp(sweepType, 'Frequency');

% Choose the subset of each
R_Cswp = ResultsTable(isCswp,:);
R_Fswp = ResultsTable(isFswp,:);


%% Contact Sweep

% boxplot the deltaH's data according to DBS-combo groupings: Contact-sweep
Contact = R_Cswp.dbsElectrode(:);
[Clabels, iC] = sort(Contact);
Hdelta = R_Cswp.H_delta(:);
HdeltaMatching = Hdelta(iC);


% find means and Standard Erorrs for each contact
uniqueLabels = unique(Clabels);
nLabs = numel(uniqueLabels);
for iLab = 1:nLabs
    isLab = strcmp(uniqueLabels(iLab), Clabels);
    data = HdeltaMatching(isLab);
    Hd_avC(iLab) = mean(data); 
    Hd_semC(iLab) = std(data)/sqrt(numel(data));
    
end

% plot the figure
figure;
bar(1:nLabs, Hd_avC); hold on
errorbar(1:nLabs, Hd_avC, Hd_semC, '.');
title('DBS-induced isi-Entropy change in All trials, Contact-Sweep');
ylabel('delta-H (bits/spike)');
xlabel('Electrode Contact (all at 130Hz)');
set(gca, 'XTickLabel', uniqueLabels);


%% Frequency Sweep

% boxplot the deltaH's data according to DBS-combo groupings: Freq-sweep
Freq = cellstr(num2str(R_Fswp.dbsFrequency(:)));
[Flabels, iF] = sort(Freq);
Hdelta = R_Fswp.H_delta(:);
HdeltaMatching = Hdelta(iF);


% find means and Standard Erorrs for each contact
uniqueLabels = unique(Flabels);
nLabs = numel(uniqueLabels);
for iLab = 1:nLabs
    isLab = strcmp(uniqueLabels(iLab), Flabels);
    data = HdeltaMatching(isLab);
    Hd_avF(iLab) = mean(data); 
    Hd_semF(iLab) = std(data)/sqrt(numel(data));
    
end

% plot the figure
figure;
bar(1:nLabs, Hd_avF); hold on
errorbar(1:nLabs, Hd_avF, Hd_semF, '.');
title('DBS-induced isi-Entropy change in All trials, Frequency-Sweep');
ylabel('delta-H (bits/spike)');
xlabel('DBS-frequency (all at C0)');
set(gca, 'XTickLabel', uniqueLabels);


% figure;
% boxplot(HdeltaMatching, Flabels);
% 
% hLine = refline(0,0); 
% hLine.LineStyle = '--';
% hLine.Color = [0,0,0];


end % END function