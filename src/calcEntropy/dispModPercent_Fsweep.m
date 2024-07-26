function dispModPercent_Fsweep(Flabels, isFlab, isFlabMod, alpha)


Ftot = sum(isFlab)';
Fmod = sum(isFlabMod)';

Fpercent(:,1) = 100 * (Fmod ./ Ftot);
Fpercent(:,2) = 100 * (1 - (Fmod ./ Ftot));


% figure; 
b = bar(Fpercent, 'stacked');
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [1, 1, 1];
set(gca, 'xticklabel', Flabels)
set(gca, 'YLim', [0, 100]);






end