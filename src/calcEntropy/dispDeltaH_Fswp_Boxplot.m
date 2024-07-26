function dispDeltaH_Fswp_Boxplot(deltaH, grpLabel, Flabels, ppar, varargin)
% Create two boxplots of deltaH values, grouped acording to labels
%% INPUTS

p = inputParser;

addRequired(p, 'deltaH');
addRequired(p, 'grpLabel');
addRequired(p, 'Flabels');
addRequired(p, 'ppar');

defaultPres = 'default';
addParameter(p, 'Pres', defaultPres, @ischar)

parse(p, deltaH, grpLabel, Flabels, ppar, varargin{:});

Pres = p.Results.Pres;


% default parameters for boxplot appearance
if ~isfield(ppar, 'boxplot')
    ppar.boxplot.dataLim = [-Inf, Inf];
    
end

% ppar.boxplot.dataLim = [-Inf, Inf];
% ppar.boxplot.dataLim = [-0.10, 0.15];



%% CODE

switch Pres
    case 'default'
        boxplot(deltaH, grpLabel, 'GroupOrder', Flabels, ...
                'Orientation', 'vertical', 'Symbol', 'o');
            
    case 'clean' % Reduce whiskers to 90 and 10 percentile, and remove outliers
        q3 = norminv(0.75);
        q90 = norminv(0.9);
        w90 = (q90 - q3) / (2 * q3);
        boxplot(deltaH, grpLabel, 'GroupOrder', Flabels, ...
                'Orientation', 'vertical', 'Whisker', w90, ...
                'Symbol', '');
            
    case 'cleanWoutliers'
        q3 = norminv(0.75);
        q90 = norminv(0.9);
        w90 = (q90 - q3) / (2 * q3);
        boxplot(deltaH, grpLabel, 'GroupOrder', Flabels, ...
                'Orientation', 'vertical', 'Whisker', w90, ...
                'Symbol', 'o', 'Jitter', 0.4, 'ExtremeMode', 'clip', ...
                'DataLim', ppar.boxplot.dataLim);
            
    otherwise
        error([Pres, ' is not a valid Name-Value pair arg!']);
        
end
    

                                  
% title('Entropy change in Modulated Contact-Sweep trials');
ylabel('delta-H (bits/spike)');
xlabel('DBS frequency (all @ C0)');
line(xlim, [0, 0], 'LineStyle', '--');


% Add in markers for tracking the mean deltaH values for each group-label
x = get(gca, 'XTick');
nLabs = numel(Flabels);
y = zeros(1, nLabs);
for iLab = 1:nLabs
    isLab = strcmp(Flabels(iLab), grpLabel);
    y(iLab) =  mean(deltaH(isLab));
    
end

hold on
plot(x, y, '-dk', 'LineWidth', 2)
hold off



end