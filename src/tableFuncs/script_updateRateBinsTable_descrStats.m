% Update RateBins table by including summary stats in table as well as link
% to files themselves

% for the log-versions, spike counts in each bin are increased by 1; this
% is to help process the many bins that have 0 spikes in them, as log(0) is
% -Inf. Since taking the log is strictly to facilitate hypothesis testing
% by data transformation/normalization of the non-normally distributed rate
% data, it's ok to do this 


clear;
subjID = 'Uva';

tabPn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
tabFnRate = ['RateBins_60sec_1secBins_', subjID, '.mat'];

% load the "RateBins" table
load([tabPn, '\', tabFnRate]);

RateBins(:,4:end) = [];



% load just first file:
nRows = size(RateBins, 1);
for iRow = 1:nRows
    
    % load data for the current row 
    ratePn = RateBins.ratebinPn{iRow};
    rateFn = RateBins.ratebinFn{iRow};
    
    % load the "bins" struct
    load([ratePn, '\', rateFn]);

    logbins.dbs = log10(bins.dbs + 1);
    logbins.pre = log10(bins.pre + 1);

    nratebinsPRE(iRow,1) = numel(bins.pre);
    nratebinsDBS(iRow,1) = numel(bins.dbs);


    % calculate descriptive statistics for data
      meanRatePRE(iRow,1) = mean(bins.pre);
      stdvRatePRE(iRow,1) = std(bins.pre);
    medianRatePRE(iRow,1) = median(bins.pre);
     qrt25RatePRE(iRow,1) = prctile(bins.pre, 25);
     qrt75RatePRE(iRow,1) = prctile(bins.pre, 75);

      meanRateDBS(iRow,1) = mean(bins.dbs);
      stdvRateDBS(iRow,1) = std(bins.dbs);
    medianRateDBS(iRow,1) = median(bins.dbs);
     qrt25RateDBS(iRow,1) = prctile(bins.dbs, 25);
     qrt75RateDBS(iRow,1) = prctile(bins.dbs, 75);


      logmeanRatePRE(iRow,1) = mean(logbins.pre);
      logstdvRatePRE(iRow,1) = std(logbins.pre);
    logmedianRatePRE(iRow,1) = median(logbins.pre);
     logqrt25RatePRE(iRow,1) = prctile(logbins.pre, 25);
     logqrt75RatePRE(iRow,1) = prctile(logbins.pre, 75);

      logmeanRateDBS(iRow,1) = mean(logbins.dbs);
      logstdvRateDBS(iRow,1) = std(logbins.dbs);
    logmedianRateDBS(iRow,1) = median(logbins.dbs);
     logqrt25RateDBS(iRow,1) = prctile(logbins.dbs, 25);
     logqrt75RateDBS(iRow,1) = prctile(logbins.dbs, 75);
    
end


% incorporate all the above stats into RateBins table
RateBins = [RateBins, table(nratebinsPRE), table(nratebinsDBS), ...
                                            table(meanRatePRE), ...
                                            table(stdvRatePRE), ...
                                            table(medianRatePRE), ...
                                            table(qrt25RatePRE), ...
                                            table(qrt75RatePRE), ...
                                            table(meanRateDBS), ...
                                            table(stdvRateDBS), ...
                                            table(medianRateDBS), ...
                                            table(qrt25RateDBS), ...
                                            table(qrt75RateDBS), ...
                                            table(logmeanRatePRE), ...
                                            table(logstdvRatePRE), ...
                                            table(logmedianRatePRE), ...
                                            table(logqrt25RatePRE), ...
                                            table(logqrt75RatePRE), ...
                                            table(logmeanRateDBS), ...
                                            table(logstdvRateDBS), ...
                                            table(logmedianRateDBS), ...
                                            table(logqrt25RateDBS), ...
                                            table(logqrt75RateDBS)];
               
                                        
save([tabPn, '\', tabFnRate], 'RateBins');








