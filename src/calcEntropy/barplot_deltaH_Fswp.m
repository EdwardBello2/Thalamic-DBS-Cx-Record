function barplot_deltaH_Fswp(ResultsTable)
% create two boxplot figures with boxplots for each DBS condition, one
% figure for each sweep-type

% Find out which rows are Contact-sweep and Freq-sweep
sweepType = ResultsTable.sweepType(:);
isCswp = strcmp(sweepType, 'Contact');
isFswp = strcmp(sweepType, 'Frequency');

% Choose the subset of each
R_Cswp = ResultsTable(isCswp,:);
R_Fswp = ResultsTable(isFswp,:);



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
% figure;
bar(1:nLabs, Hd_avF, 'k'); hold on
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