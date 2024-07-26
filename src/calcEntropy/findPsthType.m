function [isAnti, isOrtho, isOther, isInhib] = findPsthType(R)

isAnti = strcmp('anti', R.psthType(:)) | ...
         strcmp('anti-ortho', R.psthType(:));

isOrtho = strcmp('ortho', R.psthType(:));


isInhib = strcmp('inhib', R.psthType(:));

isOther = ~(isAnti | isOrtho | isInhib);


end