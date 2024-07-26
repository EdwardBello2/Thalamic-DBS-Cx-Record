% script to inspect cell-by-cell frequency response in Entropy changes

cellIDs = unique(R_ISI_Fswp.Unit_objectID);

cellNum = 3;
nCells = numel(cellIDs);
for cellNum = 1:nCells
    idx = strcmp(cellIDs{cellNum,1}, R_ISI_Fswp.Unit_objectID);
    R = R_ISI_Fswp;

    CellTab = R(idx,:);
    deltaHC = calcDeltaEntropyISIord1_allTableRows(CellTab);
    preHC = extractH1pre(CellTab);
    Entropy_AllNeu = 100*(deltaHC./preHC);


    figure; scatter(CellTab.dbsFrequency, Entropy_AllNeu)
    title(cellIDs{cellNum,1});

end