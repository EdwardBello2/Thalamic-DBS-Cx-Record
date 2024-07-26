% script for generating .csv file to be read into SAS analysis checking for
% firing rate dependence on DBS frequency:
% 1) with a One-way ANOVA for ANY differences as well as 
% 2) an orthogonal polynomial contrast hypothesis test for looking at
% linear, quadratic, and cubic effects of DBS frequency on firing rate

% first run the old script for getting firing rates:
% script_displayFRvsDBSfrequency_FRdbson_bar.m

script_displayFRvsDBSfrequency_FRdbson_bar;

% next, gather up the data into a simple table, adding in the baseline
% firing rates with "0" as the Hz used (to keep the different levels of the
% factor "dbsfrequency" numeric rather than categorical)
FRvar = 'FRdbsCorr2';

Tsas = [T(:,'dbsFrequency'), T(:,FRvar)];
Tsas.Properties.VariableNames = {'dbsfrequency', 'neuData'};
baseHz = zeros(numel(baseFR), 1);
Tbase = table(baseHz, baseFR);
Tbase.Properties.VariableNames = {'dbsfrequency', 'neuData'};

Tsas = [Tsas; Tbase];

figure; scatter(Tsas.dbsfrequency, Tsas.neuData)


% Save resulting table as .csv for analysis in SAS
savepn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataAnalysis\SAS\';
writetable(Tsas, [savepn 'FiringRates.csv']);
