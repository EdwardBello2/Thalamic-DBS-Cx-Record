function dispDeltaH_Cswp_BoxplotGroups(deltaH, grpLabel, Clabels)
% Create a boxplot of deltaH values, grouped acording to Clabels and cell
% response type


boxplot(deltaH, grpLabel, 'Orientation', 'horizontal', ...
        'factorseparator', 1, 'factorgap', 1, 'colorgroup',grpLabel(2)', ...
        'labelverbosity','majorminor');
% title('Entropy change in Modulated Contact-Sweep trials');
xlabel('delta-H (bits/spike)');
ylabel('Electrode Contact (all at 130Hz)');
line([0, 0], ylim, 'LineStyle', '--');


% Add in markers for tracking the mean deltaH values for each group-label
y = get(gca, 'YTick');
nLabs = numel(Clabels);
x = zeros(1, nLabs);
for iLab = 1:nLabs
    isLab = strcmp(Clabels(iLab), grpLabel);
    x(iLab) =  mean(deltaH(isLab));
    
end

hold on
plot(x, y, 'dk', 'LineWidth', 2)
hold off





end