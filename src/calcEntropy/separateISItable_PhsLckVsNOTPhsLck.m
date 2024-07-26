function [ISI_PhsLckNeu, ISI_NOTPhsLckNeu, PSTH_PhsLckNeu] = separateISItable_PhsLckVsNOTPhsLck(ISI, PSTH, alpha)



% ISI_NOTPhsLckNeu = extractModNeuronRows(ISI, alpha);
PSTH_PhsLckNeu = extractModNeuronRows(PSTH, alpha);


% Get Rows from ISI table pertaining to Phase-locked neurons ONLY
[~, idxPhsLck] = ismember(PSTH_PhsLckNeu.objectID, ISI.objectID);
ISI_PhsLckNeu = ISI(idxPhsLck,:);

% Get Rows from ISI table NOT pertaining to Phase-locked neurons ONLY
[~, idxNOTPhsLck] = setdiff(ISI.objectID, PSTH_PhsLckNeu.objectID);
ISI_NOTPhsLckNeu = ISI(idxNOTPhsLck,:);


% % Get Final version of PatternMod-ONLY neurons with common neurons removed
% objectID_both = intersect(ISI_NOTPhsLckNeu.objectID, ...
%                           R_PSTH_PhsLckNeu_Cswp.objectID);
% [~, idxBoth] = ismember(objectID_both, ISI_NOTPhsLckNeu.objectID);
% ISI_NOTPhsLckNeu(idxBoth, :) = [];




end