% script for comparing %deltaH in LFS vs. HFS 


close all; clear;

% first run the analysis script below:
script_runEntropyDirectISIord1_All_userVarH;
% outputs 2 variables Entropy_AllNeu and Entropy_grpLabel


% Choose the two DBS frequency labels from the data for comparison
% strings: 10 | 20 | 30 | 50 | 100 | 130
LFSfreq1 = '10'; 
LFSfreq2 = '20';

HFSfreq1 = '100';
HFSfreq2 = '130';


% Sift out all frequencies but the two choices
isLFS1 = strcmp(LFSfreq1, Entropy_grpLabel);
isLFS2 = strcmp(LFSfreq2, Entropy_grpLabel);

isHFS1 = strcmp(HFSfreq1, Entropy_grpLabel);
isHFS2 = strcmp(HFSfreq2, Entropy_grpLabel);

% isInclude = isLFS1 | isLFS2 | isHFS1 | isHFS2;
% 
% H_include = Entropy_AllNeu(isInclude);
% Freq_include = Entropy_grpLabel(isInclude);


% combine lows and combine highs
isLFS = isLFS1 | isLFS2;
isHFS = isHFS1 | isHFS2;


% Divide Entropy data between the two groups
Hlfs = Entropy_AllNeu(isLFS);
Hhfs = Entropy_AllNeu(isHFS);

[hTtest, pTtest, ci, statsTtest] = ttest2(Hlfs, Hhfs)

[pRanksum, hRanksum, statsRanksum] = ranksum(Hlfs, Hhfs)


