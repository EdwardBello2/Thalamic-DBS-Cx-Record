function dispModPercent_Csweep(Clabels, isClab, isClabMod, alpha)

Ctot = sum(isClab)';
Cmod = sum(isClabMod)';

Cpercent(:,1) = 100 * (Cmod ./ Ctot);
Cpercent(:,2) = 100 * (1 - (Cmod ./ Ctot));


% figure; 
b = barh(Cpercent, 'stacked');
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [1, 1, 1];
set(gca, 'yticklabel', Clabels)
set(gca, 'XLim', [0, 100]);



end