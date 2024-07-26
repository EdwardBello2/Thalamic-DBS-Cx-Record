function name = buildDataSelectTableName(pipeParams)
% assembles user-specified parameters for the pipeline into on long string
% to be appended to things such as table names


%%

nStr = 1; % keep track of how many strings to join into one

% See what labels are user-specified
if isfield(pipeParams, 'subjID')
    label{nStr,1} = pipeParams.subjID;
    nStr = nStr+1;
    
end

if isfield(pipeParams, 'neuTypeFilter')
    label{nStr,1} = pipeParams.neuTypeFilter;
    nStr = nStr+1;
    
end

if isfield(pipeParams, 'hzThresh')
    label{nStr,1} = [num2str(pipeParams.hzThresh), 'HzThresh'];
    nStr = nStr+1;
    
end


% Put all the labels together
name = char; % start with empty string array

for iStr = 1:(nStr - 1)
    name = [name, '_', label{iStr,1}];
    
end


end