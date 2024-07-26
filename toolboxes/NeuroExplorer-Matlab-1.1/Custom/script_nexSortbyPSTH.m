
addpath('C:\Users\bello043\OneDrive\Tools\NEXtools')

nexpath = 'D:\PROJECTS\Thalamic DBS Cx Record\DataProcessing\170825session01\';
session = '170825session01'
ch = '16'
block = '08'

nexFile = readNexFile([nexpath session '_SPKch' ch '_Block' block '_sorted.nex']);



edgeLeft = 1/1000; % seconds
winWidth = 5/1000; % seconds
prethreshTime = 5.3000e-04; % seconds

nexFileNEW = nexSortbyPSTH(nexFile,edgeLeft,winWidth,prethreshTime)

savepath = [nexpath session '_SPKch' ch '_Block' block '_PSTHsorted'];
writeNexFilev2(nexFileNEW, [savepath '.nex']);