function dispDeltaH_Cswp_Boxplot(deltaH, grpLabel, Clabels, varargin)
% Create two boxplots of deltaH values, grouped acording to labels
%% INPUTS

p = inputParser;

addRequired(p, 'deltaH');
addRequired(p, 'grpLabel');
addRequired(p, 'Clabels');

defaultPres = 'default';
addParameter(p, 'Pres', defaultPres, @ischar)

parse(p, deltaH, grpLabel, Clabels, varargin{:});

Pres = p.Results.Pres;



%% CODE

switch Pres
    
    case 'default'
        boxplot(deltaH, grpLabel, 'GroupOrder', Clabels, ...
                'Orientation', 'horizontal');
            
    case 'clean'
        q3 = norminv(0.75);
        q90 = norminv(0.9);
        w90 = (q90 - q3) / (2 * q3);
        boxplot(deltaH, grpLabel, 'GroupOrder', Clabels, ...
                'Orientation', 'horizontal', 'Whisker', w90, ...
                'Symbol', '');
            
    case 'cleanWoutliers'
        q3 = norminv(0.75);
        q90 = norminv(0.9);
        w90 = (q90 - q3) / (2 * q3);
        boxplot(deltaH, grpLabel, 'GroupOrder', Clabels, ...
                'Orientation', 'horizontal', 'Whisker', w90, ...
                'Symbol', 'o', 'Jitter', 0.4, 'ExtremeMode', 'clip', ...
                'DataLim', [-0.10, 0.15]);
            
    otherwise
        error([Pres, ' is not a valid Name-Value pair arg!']);

end
            
            
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
plot(x, y, '-dk', 'LineWidth', 2)
hold off





end