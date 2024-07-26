function dispModPercent_ISIvsPSTH_Fswp(Flabels, isFlab, isFlabPhLck, isFlabPattM)


Ftot = sum(isFlab, 1)';
Fphlck = sum(isFlabPhLck, 1)';
FptMd = sum(isFlabPattM, 1)';

Fpercent(:,1) = 100 * (Fphlck ./ Ftot);
Fpercent(:,2) = 100 * (1 - ((Fphlck + FptMd)  ./ Ftot));
Fpercent(:,3) = 100 * (FptMd ./ Ftot);



% figure; 
b = bar(Fpercent, 'stacked');
b(1).FaceColor = [0, 0, 0];
b(2).FaceColor = [1, 1, 1];
b(3).FaceColor = [0.5, 0.5, 0.5];
set(gca, 'xticklabel', Flabels)
set(gca, 'YLim', [0, 100]);







end