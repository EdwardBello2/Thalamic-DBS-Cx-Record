function dispModPercent_ISIvsPSTH_Cswp(Clabels, isClab, isClabPhLck, isClabPattM)

Ctot = sum(isClab, 1)';
Cphlck = sum(isClabPhLck, 1)';
CptMd = sum(isClabPattM, 1)';

Cpercent(:,1) = 100 * (Cphlck ./ Ctot);
Cpercent(:,2) = 100 * (1 - ((Cphlck + CptMd)  ./ Ctot));
Cpercent(:,3) = 100 * (CptMd ./ Ctot);


% figure; 
b = barh(Cpercent, 'stacked');
b(1).FaceColor = [0, 0, 0];
b(2).FaceColor = [1, 1, 1];
b(3).FaceColor = [0.5, 0.5, 0.5];
set(gca, 'yticklabel', Clabels)
set(gca, 'XLim', [0, 100]);



end