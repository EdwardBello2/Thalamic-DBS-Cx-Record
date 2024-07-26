function createCollisionBlockTableIfNeeded(NEX_2analyze, AnalysisTableName, pipeParams);
% Create the table of Collision-block analysis data for further analysis,
% but only if it a) does not yet exist, or b) if user specified to 
% overwrite existing table.


if ~exist([pipeParams.tablepn, '\', AnalysisTableName, '.mat']) || ...
            pipeParams.overwriteTables

    disp(['Did not find ', AnalysisTableName, '.mat'])
    disp('Creating...')
    
    % Perform PSTH-Entropy estimate on each row
    [CollisionBlockResults] = analyzeCollisionBlock_batch(NEX_2analyze, pipeParams); %<---------


    % SAVE final result, as both .mat and .xlsx
    % save the table
    save([pipeParams.tablepn, '\', AnalysisTableName, '.mat'], 'CollisionBlockResults');

    

    disp('CREATED');

else
    disp(['FOUND ', AnalysisTableName, '.mat'])

end % END if ~exist




end
