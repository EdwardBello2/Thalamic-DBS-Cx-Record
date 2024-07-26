function createSelectionTableIfNeeded(NEX, Sort, IntTableName, pipeParams)
% Create the table of selected data for further analysis, but only if it 
% a) does not yet exist, or b) if user specified to overwrite existing

if ~exist([pipeParams.tablepn, '\', IntTableName, '.mat']) || ...
        pipeParams.overwriteTables
    
    disp(['Did not find ', IntTableName, '.mat'])
    disp('Creating...')

    
    % Create a table of selected rows
    Selected = createSelectedDataTable(NEX, Sort, pipeParams);

    % Match each row with DBS Contact-Sweep or Freq-Sweep info, and remove any
    % rows that pertain to neither
    NEX_2analyze = updateNEXprocfiles_subjID_4isiEntropyAnalysis(Selected, TrialInfo);


    % SAVE final result, as both .mat and .xlsx
    save([pipeParams.tablepn, '\', IntTableName, '.mat'], 'NEX_2analyze');

    writetable(NEX_2analyze, [pipeParams.tablepn, '\', IntTableName, '.xlsx'])



    disp('CREATED')
    
else
    disp(['FOUND ', IntTableName, '.mat'])

end % END if ~exist

end