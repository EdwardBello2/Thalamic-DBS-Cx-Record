function boxplot_deltaH_Cswp(ResultsTable)
% create two boxplot figures with boxplots for each DBS condition, one
% figure for each sweep-type

% Find out which rows are Contact-sweep and Freq-sweep
sweepType = ResultsTable.sweepType(:);
isCswp = strcmp(sweepType, 'Contact');
isFswp = strcmp(sweepType, 'Frequency');

% Choose the subset of each
R_Cswp = ResultsTable(isCswp,:);
R_Fswp = ResultsTable(isFswp,:);

% boxplot the deltaH's data according to DBS-combo groupings: Contact-sweep
Contact = R_Cswp.dbsElectrode(:);
[Clabels, iC] = sort(Contact);
Hdelta = R_Cswp.H_delta(:);
HdeltaMatching = Hdelta(iC);

% figure;
boxplot(HdeltaMatching, Clabels);
title('DBS-induced isi-Entropy change in All trials, Contact-Sweep');
ylabel('delta-H (bits/spike)');
xlabel('Electrode Contact (all at 130Hz)');
hLine = refline(0,0); 
hLine.LineStyle = '--';
hLine.Color = [0,0,0];

% 
% 
% % boxplot the deltaH's data according to DBS-combo groupings: Freq-sweep
% Freq = cellstr(num2str(R_Fswp.dbsFrequency(:)));
% [Flabels, iF] = sort(Freq);
% Hdelta = R_Fswp.H_delta(:);
% HdeltaMatching = Hdelta(iF);
% 
% figure;
% boxplot(HdeltaMatching, Flabels);
% title('DBS-induced isi-Entropy change in All trials, Frequency-Sweep');
% ylabel('delta-H (bits/spike)');
% xlabel('DBS-frequency (all at C0)');
% hLine = refline(0,0); 
% hLine.LineStyle = '--';
% hLine.Color = [0,0,0];


end % END function


