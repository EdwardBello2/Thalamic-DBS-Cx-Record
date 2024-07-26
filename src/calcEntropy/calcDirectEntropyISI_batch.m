% Function for calculating Direct Entropy estimate of 
%
% INPUTS:
% N - matlab table of values to be analyzed
% pipeParams - struct with user-spec analysis parameters

% TO-DO:
%

function [Hdir_ISI_Results] = calcDirectEntropyISI_batch(N, pipeParams)


% in seconds. Time to include spikes for analysis before DBS onset
PREDBS_TIME = pipeParams.predbsTime; 
DBS_TIME = pipeParams.dbsTime;  % also the time that DBS was on for the experiment 

% PARAMETERS for Direct Entropy estimate with log binned ISI:
ORD_H = pipeParams.ordH; % how many orders of entropy to use for linear extrapolation
BINS_PER_DECADE = pipeParams.binsPD; % of log-spaced bins for ISI histogram



%%  MAIN FOR-LOOP

nNex = size(N,1);

H_PRE = zeros(nNex,1);
H_DBS = zeros(nNex,1);
H_delta = zeros(nNex,1);

disp('Calculating Direct Entropy estimates...')
tic
for iNex = 1:nNex

    % Get each NEX file from table
    nexfn = N.Filename{iNex,1};
    nexpn = N.Pathname{iNex,1};

    nexFile = readNexFile([nexpn, '\', nexfn]);
    [spkTimes, StimTs] = parseNexFile(nexFile);

    % separate the spike times into pre-DBS and DBS-on
    dbsTimes = StimTs.DBS;
    stimPeriod = median(diff(dbsTimes));


    % get pre-DBS spikes
    isPreDBS = (spkTimes < dbsTimes(1)) & ...
               (spkTimes >= (dbsTimes(1) - PREDBS_TIME));
    spksPRE = spkTimes(isPreDBS);


    % get DBS-on spikes
    isDBSon = (spkTimes >= dbsTimes(1)) & ...
              (spkTimes < (dbsTimes(end) + stimPeriod));
    spksDBS = spkTimes(isDBSon);



    %% Determine which condition has fewer nSpks

    if numel(spksPRE) < numel(spksDBS)
        nSpks = numel(spksPRE);
        hasMore = 'DBS';

    elseif numel(spksPRE) > numel(spksDBS)
        nSpks = numel(spksDBS);
        hasMore = 'PRE';

    else % if they are equel number spikes
        nSpks = numel(spksDBS);
        hasMore = 'NONE';

    end



    %% Calculate Entropies and differences
    % depending on which condition has more spks, resample and get the av Entropy

    switch hasMore

        case 'NONE'  % if both have same nSpks:
            % Calculate ORD_H orders of ISI entropy and fit line to find
            % extrapolated H for each
            H_DBS(iNex,1) = extrapEntropyDirect(spksDBS, BINS_PER_DECADE, ORD_H);
            H_PRE(iNex,1) = extrapEntropyDirect(spksPRE, BINS_PER_DECADE, ORD_H);


        case 'DBS' % if DBS condition had more spikes:
            % Resample DBS spikes many times so that number of spikes
            % matches the PRE case, then find the average. Calc PRE case normally 
            spksDBSresampled = groupSubsequent(spksDBS, nSpks)';
            
            nResamp = size(spksDBSresampled, ORD_H);
            H_Resamp = zeros(nResamp, ORD_H);
            
            for iRes = 1:nResamp
                
                H_Resamp(iRes,:) = entropyDirect_multOrder(spksDBSresampled(:,iRes), ...
                                                        BINS_PER_DECADE, ORD_H);

            end
            
            H_av = mean(H_Resamp);
            
            H_PRE(iNex,1) = extrapEntropyDirect(spksDBS, BINS_PER_DECADE, ORD_H);
            H_DBS(iNex,1) = fitLinear_HnOrd(H_av);


        case 'PRE' % if PRE condition had more spikes:
            % Resample PRE spikes many times so that number of spikes
            % matches the DBS case, then find the average. Calc DBS case normally 
            spksPREresampled = groupSubsequent(spksPRE, nSpks)';
            
            nResamp = size(spksPREresampled, ORD_H);
            H_Resamp = zeros(nResamp, ORD_H);
            
            for iRes = 1:nResamp
                
                H_Resamp(iRes,:) = entropyDirect_multOrder(spksPREresampled(:,iRes), ...
                                                        BINS_PER_DECADE, ORD_H);

            end
            
            H_av = mean(H_Resamp);
            
            H_PRE(iNex,1) = fitLinear_HnOrd(H_av);
            H_DBS(iNex,1) = extrapEntropyDirect(spksDBS, BINS_PER_DECADE, ORD_H);

            
        otherwise
            error('variable "hasMore" is not as expected')
            
    end % END switch



    H_delta(iNex,1) = H_DBS(iNex,1) - H_PRE(iNex,1);


end % END FOR-LOOP
toc
disp('Done Calculating!')

% halelujah!
load handel;
player = audioplayer(y, Fs);
play(player);



%% Create RESULTS table and SAVE

P = table(H_PRE);
D = table(H_DBS);
A = table(H_delta);

Hdir_ISI_Results = [N, P, D, A];

% % store the params that I used to generate this result in the table itself
% Hdir_ISI_Results.Properties.UserData.preDBStime = PREDBS_TIME;
% Hdir_ISI_Results.Properties.UserData.DBStime = DBS_TIME;
% Hdir_ISI_Results.Properties.UserData.totEntropyOrders = ORD_H;
% Hdir_ISI_Results.Properties.UserData.binsPerDecade = BINS_PER_DECADE;





disp('SUCCESS!')






end











