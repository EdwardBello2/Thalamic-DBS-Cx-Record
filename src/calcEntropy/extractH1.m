function H1 = extractH1(H_multOrd)
% Takes in cell array with each row containing a row-vector of multiorder
% Entropy estimates. Outputs a corresponding numerical column-vector of
% first-order Entropy only

nRows = size(H_multOrd, 1);
H1 = zeros(nRows, 1);
for iRow = 1:nRows
    H1(iRow,1) = H_multOrd{iRow,1}(1,1);
    
end

end
