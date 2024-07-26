function runPipeline(pipelineID, pipeParams)
% Master pipeline entropy-point function for the project:
% Thalamic Stim Cx Rec
%
% Detailed instructions on the pipeline parameter inputs needed for each
% pipeline coming soon...
% 
% TO-DO in general:
% - Need to make sure all sweep analysis screens out things that don't
% belong in the sweep. 


switch pipelineID
    
    case 'EntropyDirectISI'
        runEntropyDirectISI(pipeParams);
        
    case 'EntropyPSTHletter'
        runEntropyPSTHletter(pipeParams);
        
    case 'EntropyPSTHletter_byPSTHclass'
        runEntropyPSTHletter_byPSTHclass(pipeParams);
        
    case 'EntropyDirectISI_byPSTHclass'
        runEntropyDirectISI_byPSTHclass(pipeParams);
        
    case 'EvokedSpikeLatency'
        runEvokedSpikeLatency(pipeParams);
        
    case 'CompareEntropy_ISIvsPSTH'
        runCompareEntropy_ISIvsPSTH(pipeParams);
        
    case 'EntropyPSTHletter_byCell'
        runEntropyPSTHletter_byCell(pipeParams);
        
    case 'EntropyDirectISI_byCell'
        runEntropyDirectISI_byCell(pipeParams);
        
    case 'EntropyDirectISI_PhLckVsPattMd_byCell'
        runEntropyDirectISI_PhLckVsPattMd_deltaH_byCell(pipeParams);
        
    case 'EntropyDirectISIord1_byCell'
        runEntropyDirectISIord1_byCell(pipeParams);
        
    case 'EntropyDirectISIord1_PhLckVsPattMd_byCell'
        runEntropyDirectISIord1_PhLckVsPattMd_deltaH_byCell(pipeParams);
        
    case 'EntropyDirectISIord1_PhLckVsNonPhLckVsAll_deltaH_byCell'
        runEntropyDirectISIord1_PhLckVsNonPhLckVsAll_deltaH_byCell(pipeParams);
    
    case 'EntropyDirectISIord1_PhLckVsPattMdVsAll_deltaH_byCellONEplot'
        runEntropyDirectISIord1_PhLckVsPattMdVsAll_deltaH_byCellONEplot(pipeParams);
        
    case 'EntropyDirectISIord1_PhLckVsPattMdVsAll_dbsH_byCell'
        runEntropyDirectISIord1_PhLckVsPattMdVsAll_dbsH_byCell(pipeParams);
        
    case 'EntropyPSTH_percentPopPhaseLocked'
        runEntropyPSTH_percentPopPhaseLocked(pipeParams);
        
    case 'EntropyDirectISIord1_All_userVarH'
        runEntropyDirectISIord1_All_userVarH(pipeParams);

        
    otherwise
        error('!! pipelineID not recognized !!')
        
end
    








end