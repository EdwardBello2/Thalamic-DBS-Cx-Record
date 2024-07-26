function [ISI_PhsLckNeu, ISI_PatModNeu] = separateISItable_PhsLckVsPatMod(ISI, PSTH, alpha)



ISI_PatModNeu = extractModNeuronRows(ISI, alpha);
R_PSTH_PhsLckNeu_Cswp = extractModNeuronRows(PSTH, alpha);


% Get Rows from ISI table pertaining to Phase-locked neurons ONLY
[~, idxPhsLck] = ismember(R_PSTH_PhsLckNeu_Cswp.objectID, ISI.objectID);
ISI_PhsLckNeu = ISI(idxPhsLck, :);


% Get Final version of PatternMod-ONLY neurons with common neurons removed
objectID_both = intersect(ISI_PatModNeu.objectID, ...
                          R_PSTH_PhsLckNeu_Cswp.objectID);
[~, idxBoth] = ismember(objectID_both, ISI_PatModNeu.objectID);
ISI_PatModNeu(idxBoth, :) = [];




end