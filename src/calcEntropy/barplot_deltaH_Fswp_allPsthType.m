function barplot_deltaH_Fswp_allPsthType(ResultsTable)
% create barplot for all PSTH types in one figure

% Find out which rows are Contact-sweep and Freq-sweep
sweepType = ResultsTable.sweepType(:);
isCswp = strcmp(sweepType, 'Contact');
isFswp = strcmp(sweepType, 'Frequency');

% Choose the subset of each
% R_Fswp = ResultsTable(isCswp,:);
R_Fswp = ResultsTable(isFswp,:);


%% Contact Sweep

% Divide table into groupings based on PSTH type
psthType = {'anti', 'ortho', 'other', 'inhib'};
[isPsthType(:,1), isPsthType(:,2), isPsthType(:,3), isPsthType(:,4)] = findPsthType(R_Fswp);

[Flabels, isFlabel] = findDBSlabels_Fswp(R_Fswp);
nLabs = numel(Flabels);



%% For Antidromics:


Hd_avF = zeros(numel(Flabels), 4);
Hd_semF = zeros(numel(Flabels), 4);

for jPsthType = 1:4
    
    R_psthType = R_Fswp(isPsthType(:,jPsthType),:);

    [Hd_avF(:,jPsthType), Hd_semF(:,jPsthType)] = findMeanSem_perLabel(R_psthType, Flabels);


end




%% Plot out groups of PSTH type with error bars:


ax = axes;
b = bar(Hd_avF,'BarWidth',1);
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,[1:nLabs]);

% Naming each of the bar groups
xticklabels(ax,Flabels);

% Figure labels
ylabel('delta-H (bits/spike)');
xlabel('DBS Frequency (all at CO)');
title('DBS-induced isi-Entropy change in All trials, Frequency-Sweep');

% Creating a legend and placing it outside the bar plot
lg = legend(psthType,'AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;

% Finding the number of groups and the number of bars in each group
ngroups = size(Hd_avF, 1);
nbars = size(Hd_avF, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Hd_avF(:,i), Hd_semF(:,i), 'k', 'linestyle', 'none');
end


end % END function

function [Hd_avC, Hd_semC] = findMeanSem(HdeltaMatching, Clabels)

uniqueLabels = unique(Clabels);
nLabs = numel(uniqueLabels);


for iLab = 1:nLabs
    isLab = strcmp(uniqueLabels(iLab), Clabels);
    data = HdeltaMatching(isLab);
    Hd_avC(iLab,1) = mean(data); 
    Hd_semC(iLab,1) = std(data)/sqrt(numel(data));
    
end



end

function [Hd_avC, Hd_semC] = findMeanSem_perLabel(R, uniqueLabels)
% Output column of

% uniqueLabels = unique(Clabels);
labels = R.dbsFrequency(:);
nLabs = numel(uniqueLabels);

Hd_avC = zeros(nLabs, 1);
Hd_semC = zeros(nLabs, 1);
for iLab = 1:nLabs
    isLab = uniqueLabels(iLab) == labels;
    data = R.H_delta(isLab);
    Hd_avC(iLab,1) = mean(data); 
    Hd_semC(iLab,1) = std(data)/sqrt(numel(data));
    
end



end