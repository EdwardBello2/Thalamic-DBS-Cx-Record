function barplot_deltaH_Cswp(ResultsTable)
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
% figure;
bar(1:nLabs, Hd_avC, 'k'); hold on
errorbar(1:nLabs, Hd_avC, Hd_semC, '.');
title('DBS-induced isi-Entropy change in All trials, Contact-Sweep');
ylabel('delta-H (bits/spike)');
xlabel('Electrode Contact (all at 130Hz)');
set(gca, 'XTickLabel', uniqueLabels);





end % END function