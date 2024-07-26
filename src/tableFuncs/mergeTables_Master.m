function MergedTable = mergeTables_Master(callingScript, ppar)
% function that, for a given calling function, I have set up to call the
% appropriate tables and merge them in the appropriate way for that calling
% function's operations. This function is a way for me to have all the
% table loads and merges in one place, so that any updates to tables and
% changes to how they're dealt with can all be updated in one convenient
% location, as opposed to me having to look thru each and every script that
% loads and merges tables as a setup stage before operations. 


% User must specify pipeParams fields for specific ways to merge tables:

% Author: Ed Bello
% Created: 2019/07/04

% script to call for merging specific tables as an intermediate step in an
% analysis pipeline


% User must specify ppar fields for specific ways to merge tables:

%% pipelineParams used for this script:

% ppar.tablePath





switch callingScript

    case {'script_CompareRatePsthIsi', ...
          'script_CompareRatePsthIsi_specificIteration', ...
          'script_displayRateVsEntropy_scatter', ...
          'script_displayPopRates_PREandDBS_barplot', ...
          'script_displayRateChangeCategory_barplotStack', ...
          'script_plot3uniqueNeuronLogISIs_freqSwp_PREandDBSon'}
      % FULL table
      
        % Join the ROOT table with RateBins analysis table
        load([ppar.tablePath, '\', 'tableRoot'], 'tableRoot');
        RateBins = loadTable_RateBins_60sec_1secBins(ppar);

        N = join(tableRoot, RateBins, 'LeftKeys', 'objectID', 'RightKeys', 'objectID');


        % Join the two Entropy tables
        Hisi = loadTable_EntropyDirectISI_Analysis(ppar);
        Hpsth = loadTable_EntropyPsthLetter_Analysis(ppar);
       
        H = join(Hisi, Hpsth, 'LeftKeys', 'objectID', 'RightKeys', 'objectID');


        % Final join
        MergedTable = join(H, N, 'LeftKeys', 'objectID', 'RightKeys', 'objectID');

           
    case {'script_displayPhaseLockEntropyChangeVsdbsCond_boxplots', ...
          'script_displayPhaseLockNeuronPerc_pieChart', ...
          'script_displayPhaseLockPercVsdbsTrialCond_barplot', ...
          'script_displayPercPop_isiAndPhaseLckCategory_pieChart'}
      % Entrop only (PSTH and logISI)
      
        load([ppar.tablePath, '\', 'tableRoot'], 'tableRoot');

        % Join the two Entropy tables
        Hisi = loadTable_EntropyDirectISI_Analysis(ppar);
        Hpsth = loadTable_EntropyPsthLetter_Analysis(ppar);
       
        H = join(Hisi, Hpsth, 'LeftKeys', 'objectID', 'RightKeys', 'objectID');
        

        % Final join
        MergedTable = join(H, tableRoot, 'LeftKeys', 'objectID', 'RightKeys', 'objectID');
    
    case {'script_displayPhaseLockEntropyChangeVsdbsCond_boxplots_v30s', ...
          'script_displayAverageNeuronEntropyByCondGrp_boxplots'}
        % Load the table
        % "EntropyDirectISI_Analysis_SU_2hzThresh_30pre_30dbs_ordH1_15binsPD_1000boots.mat"
        % instead of the old 60-second version of entropy calculation
        load([ppar.tablePath, '\', 'tableRoot'], 'tableRoot');
        
        % Join the two Entropy tables
        Hisi = loadTable_EntropyDirectISI_Analysis_v30s(ppar);
        Hpsth = loadTable_EntropyPsthLetter_Analysis_v30s(ppar);
%        
        H = join(Hisi, Hpsth, 'LeftKeys', 'objectID', 'RightKeys', 'objectID');

        % Final join
        MergedTable = join(H, tableRoot, 'LeftKeys', 'objectID', 'RightKeys', 'objectID');
    
%     case {'script_displayAverageNeuronEntropyByCondGrp_boxplots'}
%         % load root table
%         load([ppar.tablePath, '\', 'tableRoot'], 'tableRoot');
%         
%         % load unique Neuron table
%         load([ppar.tablePath, '\', 'SortedUnits'], 'SortedUnits');
% 
%         % join them
%         MergedTable = join(tableRoot, SortedUnits, 'LeftKeys', 'Unit_objectID', ...
%                            'RightKeys', 'objectID');
        
   
    otherwise
        error('bla')
        
        
end




end

%% SUB-FUNCTIONS

function RateBins = loadTable_RateBins_60sec_1secBins(ppar)
% load table and rename variables

% load RateBins_60sec_1secBins table
load([ppar.tablePath, '\', 'RateBins_60sec_1secBins'], 'RateBins');

% rename matFile and matFileFolder variables:
varNames = RateBins.Properties.VariableNames;
varNames{strcmp(varNames, 'matFile')} = 'matFile_RateBins';
varNames{strcmp(varNames, 'matFileFolder')} = 'matFileFolder_RateBins';
RateBins.Properties.VariableNames = varNames;
        
end

function Hisi = loadTable_EntropyDirectISI_Analysis(ppar)
% load table and rename variables

% load EntropyDirectISI table:
tabName = 'EntropyDirectISI_Analysis_SU_2HzThresh_60pre_60dbs_ordH1_15binsPD_10000boots';
load([ppar.tablePath, '\', tabName], 'EntropyDirectISI_analysis');
Hisi = EntropyDirectISI_analysis;

% rename matFile and matFileFolder variables:
varNames = Hisi.Properties.VariableNames;
varNames{strcmp(varNames, 'matFile')} = 'matFile_Hisi';
varNames{strcmp(varNames, 'matFileFolder')} = 'matFileFolder_Hisi';
Hisi.Properties.VariableNames = varNames;

end

function Hisi = loadTable_EntropyDirectISI_Analysis_v30s(ppar)
% load table and rename variables

% load EntropyDirectISI table:
tabName = 'EntropyDirectISI_Analysis_SU_2hzThresh_30pre_30dbs_ordH1_15binsPD_10000boots';
load([ppar.tablePath, '\', tabName], 'EntropyDirectISI_analysis');
Hisi = EntropyDirectISI_analysis;

% rename matFile and matFileFolder variables:
varNames = Hisi.Properties.VariableNames;
varNames{strcmp(varNames, 'matFile')} = 'matFile_Hisi';
varNames{strcmp(varNames, 'matFileFolder')} = 'matFileFolder_Hisi';
Hisi.Properties.VariableNames = varNames;

end

function Hpsth = loadTable_EntropyPsthLetter_Analysis(ppar)
% load table and rename variables

% load EntropyPsthLetter table:
tabName = 'EntropyPsthLetter_Analysis_SU_2HzThresh_60pre_60dbs_10000boots';
load([ppar.tablePath, '\', tabName], 'EntropyPsthLetter_analysis');
Hpsth = EntropyPsthLetter_analysis;

end

function Hpsth = loadTable_EntropyPsthLetter_Analysis_v30s(ppar)
% load table and rename variables

% load EntropyPsthLetter table:
tabName = 'EntropyPsthLetter_Analysis_SU_2hzThresh_30pre_30dbs_0p5msBins_1000boots';
load([ppar.tablePath, '\', tabName], 'EntropyPsth_analysis');
Hpsth = EntropyPsth_analysis;

% rename matFile and matFileFolder variables:
varNames = Hpsth.Properties.VariableNames;
varNames{strcmp(varNames, 'matFile')} = 'matFile_Hpsth';
varNames{strcmp(varNames, 'matFileFolder')} = 'matFileFolder_Hpsth';
Hpsth.Properties.VariableNames = varNames;

end