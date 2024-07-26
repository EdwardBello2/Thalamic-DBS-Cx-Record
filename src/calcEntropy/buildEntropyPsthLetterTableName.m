function name = buildEntropyPsthLetterTableName(pipeParams)
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

if isfield(pipeParams, 'predbsTime')
    label{nStr,1} = [num2str(pipeParams.predbsTime), 'pre'];
    nStr = nStr+1;
    
end

if isfield(pipeParams, 'dbsTime')
    label{nStr,1} = [num2str(pipeParams.dbsTime), 'dbs'];
    nStr = nStr+1;
    
end


if isfield(pipeParams, 'nBoot')
    label{nStr,1} = [num2str(pipeParams.nBoot), 'boots'];
    nStr = nStr+1;
    
end




% if isfield(pipeParams, 'ordH')
%     label{nStr,1} = ['ordH', num2str(pipeParams.ordH)];
%     nStr = nStr+1;
%     
% end
% 
% if isfield(pipeParams, 'binsPD')
%     label{nStr,1} = [num2str(pipeParams.binsPD), 'binsPD'];
%     nStr = nStr+1;
%     
% end


%% Put all the labels together
name = char; % start with empty string array

for iStr = 1:(nStr - 1)
    name = [name, '_', label{iStr,1}];
    
end



end % END function