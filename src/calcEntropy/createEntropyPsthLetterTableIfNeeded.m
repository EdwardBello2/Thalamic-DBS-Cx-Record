function createEntropyPsthLetterTableIfNeeded(NEX_2analyze, AnalysisTableName, pipeParams)
% Create the table of Evoked Spike Latency data for further analysis,
% but only if it a) does not yet exist, or b) if user specified to 
% overwrite existing table.


% Check if such a table already exists
if ~exist([pipeParams.tablepn, '\', AnalysisTableName, '.mat']) || ...
            pipeParams.overwriteTables

    disp(['Did not find ', AnalysisTableName, '.mat'])
    disp('Creating...')
    
    % Perform PSTH-Entropy estimate on each row
    [H_letterResults] = calcPSTHletterEntropy_batch(NEX_2analyze, pipeParams); %<---------


    % SAVE final result, as both .mat and .xlsx
    % save the table
    save([pipeParams.tablepn, '\', AnalysisTableName, '.mat'], 'H_letterResults');

%     writetable(H_letterResults, [pipeParams.tablepn, '\', AnalysisTableName, '.xlsx']);
    

    disp('CREATED');

else
    disp(['FOUND ', AnalysisTableName, '.mat'])

end % END if ~exist


end