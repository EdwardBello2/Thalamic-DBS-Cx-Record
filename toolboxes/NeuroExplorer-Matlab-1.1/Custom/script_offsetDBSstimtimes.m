

fullfilename = 'D:\DataProcessing\170801session01\170801session01_Block02_SPKch05_sorted.nex';

nexFile = readNexFile(fullfilename)


timestamps = nexFile.events{4,1}.timestamps;


offset = 33/24414.0625;
timestamps = timestamps+offset;

nexFile.events{4,1}.timestamps = timestamps;



writeNexFile(nexFile,fullfilename)