function extrapHinfinite = fitLinear_HnOrd(nHs)
% fit a line to the n Entropies input. Output the y-intercept. nHs should
% be in ascending order, such that nHs(1) corresponds to 1st order entropy,
% nHs(n) corresponds to nth order entropy. 
%
% for reference to this technique, see Strong et al. 1998 or Dorval 2008
% Basically extrapolate from empirical estimates of Entropy orders to an
% "infinite order" entropy (ISIs) or infinitely long dataset (word-method)


nHs = fliplr(nHs);

nOrdH = 1:length(nHs);
dimrecip = 1./(fliplr(nOrdH));


% plot the preDBS Entropy estimate:
coeff = polyfit(dimrecip(end-1:end), nHs, 1);
fitLine = polyval(coeff, [0, dimrecip]);

extrapHinfinite = fitLine(1);


end