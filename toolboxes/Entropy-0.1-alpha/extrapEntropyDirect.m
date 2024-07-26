function extrapHinfinite = extrapEntropyDirect(spkTimes, binPD, ordH)
% calculate ordH orders of Entropy for spike times and find the y-intercept
% of the extrapolated line (See Strong et al 1998)





nOrdH = 1:ordH;
for iOrdH = nOrdH
    nHs(iOrdH) = entropyDirect(spkTimes, binPD, iOrdH);
end

extrapHinfinite = fitLinear_HnOrd(nHs);



end
