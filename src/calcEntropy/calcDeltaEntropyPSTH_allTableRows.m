function deltaH = calcDeltaEntropyPSTH_allTableRows(R)
% Generate a numerical column vector of delta-H values from the table of
% PSTH-letter Entropy values. 




nRows = size(R, 1);
H_PreAv = zeros(nRows, 1);
for iRow = 1:nRows
    H_PreAv(iRow,1) = mean(R.H_PREbootdistr{iRow,1});
    
end
deltaH = R.H_DBS(:) - H_PreAv;



end



