% Display a scatter of rate-change vs entropy-change 

% Author: Ed Bello
% Created: 2019/07/11

%% pipeline Parameters and script Constants

% EXAMPLE PIPELINE PARAMETERS:
% % full path on local PC where project folder is (don't include subfolders here,
% % that's tracked within the appropriate tables)
% ppar.projRootPath = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record'; % string
% 
% % full path on local PC where tables to be loaded are kept
% ppar.tablePath = 'C:\Users\Acer\Documents\GitHub\Thalamic-DBS-Cx-Record\tables'; % string
% 
% % Table selection related parameters
% ppar.subjID = 'Kramer'; % string, 'Kramer' | 'Uva'
% ppar.trialType = 'FrequencySweep'; % string, 'ContactSweep' | 'FrequencySweep'
% ppar.neuronType = 'SU'; % string, SU | MU | 
% ppar.hzThresh = 2; % Hz, numeric
% ppar.neuronMinimumTrials = 4; % integer
% 
% % NEXfile related parameters
% ppar.preDbsTime = 60; % seconds
% ppar.dbsTime = 60; % seconds
% 
% % PSTH-related parameters
% ppar.trimPSTH = true; % TF indicating whether to remove the first and last bins
% ppar.psthTimeBeg = 0; % seconds
% ppar.psthBinWidth = 0.5 / 1000; %seconds
% ppar.psthNumBins = 15;
% ppar.pAlphaPSTH = 0.05;
% 
% % log-ISI related parameters
% ppar.binsPerDecade = 15;

clear; 

% Initialize fields in ppar struct in an accompanying script:
script_pipelineParams

 % CONSTANTS
CURR_FUNC = 'script_displayRateVsEntropy_scatter'; 

COLORS = [0, 0.4470, 0.7410;       % blue 
          0.8500, 0.3250, 0.0980;  % red-orange
          0.9290, 0.6940, 0.1250;  % yellow
          0.4940, 0.1840, 0.5560;  % purple
          0.4660, 0.6740, 0.1880;  % green
          0.3010, 0.7450, 0.9330;  % cyan
          0.6350, 0.0780, 0.1840;];% red-dark
      
MARKERS = {'o';
           'd';};


%% LOAD all relevant tables and MERGE them

% Custom function for selecting and merging all tables for this analysis
Tcombo = mergeTables_Master(CURR_FUNC, ppar);



%% Specify a selection of the above joined table for final analysis

Tselect = filterTable_Master(Tcombo, CURR_FUNC, ppar);



%% Divide Rate and Entropy data according to user-defined methods and groupings

% how many groups -- note that '' is the default empty case, to be ignored
% grps = unique(Tselect.group);
% grps(strcmp(grps, '')) = [];

 


if isfield(ppar, 'groups')
    grpLabels = ppar.groups.dbsFrequency.groupLabels;
    
else
    grpLabels = unique(Tselect.group);
    
end


% Gather up data for first group first, then next, etc...
nGroups = numel(grpLabels);
for iGrp = 1:nGroups
    isCurrGrp = strcmp(Tselect.group, grpLabels{iGrp});
    TbyGroup = Tselect(isCurrGrp, :);


    % Assign Colors to group-specific data index
    grpColor = COLORS(iGrp,:);
    grpMarker = MARKERS{iGrp};
    nRows = size(TbyGroup, 1);
    dataColor = repmat(grpColor, nRows, 1);
%     dataMarker = repmat(


    % Collect Rate values for comparison
    switch ppar.rateType
        case 'meanDBS'
            dataRates = TbyGroup.meanRateDBS;
            rateTitle = 'meanDBS-Rate';
            rateAxisLabel = 'spikes/second';

        case 'meanDiff'
            dataRates = TbyGroup.meanRateDBS - TbyGroup.meanRatePRE;
            rateTitle = '\DeltameanRate';
            rateAxisLabel = '\Delta spikes/second';

        otherwise
            error('wrong input for ppar.rateType');

    end


    % Collect Entropy values for comparison
    switch ppar.entropyType
        case 'isiDbsH'
            dataEntropies = TbyGroup.H_DBSemp_Hisi;
            entropyTitle = 'logISI DBS-Entropy';
            entropyAxisLabel = 'bits/spike';
            
        case 'isiPreH'
            dataEntropies = TbyGroup.H_PREbootAv;
            entropyTitle = 'logISI baseline-Entropy';
            entropyAxisLabel = 'bits/spike';

        case 'isiDiffH'
            dataEntropies = TbyGroup.H_DBSemp_Hisi - TbyGroup.H_PREbootAv;
            entropyTitle = 'logISI \DeltaEntropy';
            entropyAxisLabel = '\Delta bits/spike';

        case 'isiPercDiffH'
            dataEntropies = (TbyGroup.H_DBSemp_Hisi - TbyGroup.H_PREbootAv) ./ ...
                        TbyGroup.H_PREbootAv;
            entropyTitle = 'logISI %\DeltaEntropy';
            entropyAxisLabel = '%\Delta ISI-Entropy';

        otherwise
            error('wrong input for ppar.entropyType');

    end

%     % Update final rate, entropy, and color labels by group
%       rates = [rates; dataRates];
%     entropies = [entropies; dataEntropies];
%       color = [color; dataColor];

    % Update final rate, entropy, and color labels by group
      ratesByGrp{iGrp} = dataRates;
    entropiesByGrp{iGrp} = dataEntropies;
      colorByGrp{iGrp} = dataColor;
      
end



%% Scatterplot by group 

% scatterplot and include linear regression line
% regression method derived from: https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html

% groupLabels = ppar.groups.dbsFrequency.groupLabels;

f1 = figure; ax1 = axes;
hold on

nGroups = numel(grpLabels);
for iGrp = 1:nGroups
scatter(ratesByGrp{iGrp}, entropiesByGrp{iGrp}, [], colorByGrp{iGrp}, MARKERS{iGrp});

end

if isfield(ppar, 'groups'), legend(grpLabels); end % don't show a legend if there are no user-defined groups

% Turn off the plots below ability to show up in the legend, using "set"
% zero-lines should be dotted
p2 = plot(ax1.XLim, [0, 0], 'k--');
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

p3 = plot([0, 0], ax1.YLim, 'k--');
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

grid on

title([ppar.subjID, ':   ', entropyTitle, '  Vs.  ', rateTitle ]);
xlabel(rateAxisLabel);
ylabel(entropyAxisLabel);



%% Linear regression fit of chosen Rate and Entropy data (all groups lumped together)

% initialize data for scatterplot as empty doubles, to be filled
  rates = double.empty();
entropies = double.empty();
  color = double.empty();
  
for iGrp = 1:nGroups
    rates = [rates; ratesByGrp{iGrp}];
    entropies = [entropies; entropiesByGrp{iGrp}];
    
end

Y = entropies;
X = rates;
X = [ones(length(X),1) X];
b = X\Y;
yCalc = X*b;
regressionSlope = b(2)

% Turn off the plots below ability to show up in the legend, using "set"
p1 = plot(rates,yCalc, 'k');
set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


%  Test if significantly correlated and for Pearson Correlation Coeff
[rho, pVal] = corr(rates, entropies)


   