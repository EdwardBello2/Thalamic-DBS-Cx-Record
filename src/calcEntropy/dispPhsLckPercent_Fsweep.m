function dispPhsLckPercent_Fsweep(Flabels, isFlab, isFlabMod, pipeParams)


Ftot = sum(isFlab)';
Fmod = sum(isFlabMod)';

Fpercent(:,1) = 100 * (Fmod ./ Ftot);
% Fpercent(:,2) = 100 * (1 - (Fmod ./ Ftot));


% figure; 
b = bar(Fpercent);
b(1).FaceColor = [0.5, 0.5, 0.5];
% b(2).FaceColor = [1, 1, 1];
set(gca, 'xticklabel', Flabels)
set(gca, 'YLim', pipeParams.PhsLckPercentYLim);






end