

%% Update metadata on SweepAnalysis Table


% parse objectID into date, session, and block

% read in one Primary Key
S = SweepAnalysisTrials4Paper2018;

clear YY MM DD sessionNum blockNum
nRows = size(S,1);



for iRow = 1:nRows
                primKey = S.objectID{iRow,1};
             YY{iRow,1} = primKey(1:2);
             MM{iRow,1} = primKey(3:4);
             DD{iRow,1} = primKey(5:6);
     sessionNum{iRow,1} = str2num(primKey(7:8));
    blockNum{iRow,1} = str2num(primKey(10:11));
    
end

S = [S, table(YY), table(MM), table(DD), table(sessionNum), table(blockNum)];

descr = 'Table with details about the DBS trials associated wtih a given date and session, as well as the "block" number (i.e. block 1 is the first DBS trial recorded for the current session, block 2 the second DBS trial, etc). "ContactSweep" and "FrequencySweep" refer to whether  or not (1 or 0) that current block is considered a Contact- or Frequency-Sweep for the purposes of analysis in this current project. Note that some blocks were not either.';
S.Properties.Description = descr;



S.Properties.VariableDescriptions{1,1} = 'Primary Key';
S.Properties.VariableDescriptions{1,2} = 'Electrode ID on the Numed lead, C0 is deepest';
S.Properties.VariableDescriptions{1,3} = 'Num of DBS pulses per second';
S.Properties.VariableDescriptions{1,4} = 'is Trial part of the ContactSweep analysis?';
S.Properties.VariableDescriptions{1,5} = 'is Trial part of the FrequencySweep analysis?';
S.Properties.VariableDescriptions{1,6} = '';
S.Properties.VariableDescriptions{1,7} = '';
S.Properties.VariableDescriptions{1,8} = '';
S.Properties.VariableDescriptions{1,9} = 'Session for that day';
S.Properties.VariableDescriptions{1,10} = 'Order of presentation for this Trial within the Session';

S.Properties.VariableUnits{1,1} = '';
S.Properties.VariableUnits{1,2} = '';
S.Properties.VariableUnits{1,3} = 'Hz';
S.Properties.VariableUnits{1,4} = '';
S.Properties.VariableUnits{1,5} = '';
S.Properties.VariableUnits{1,6} = '';
S.Properties.VariableUnits{1,7} = '';
S.Properties.VariableUnits{1,8} = '';
S.Properties.VariableUnits{1,9} = '';
S.Properties.VariableUnits{1,10} = '';


SweepAnalysisTrials4Paper2018 = S;

open SweepAnalysisTrials4Paper2018


%% Update table metadata for NexProcfiles

N = NEXprocfiles;

descr = 'Table containing the filenames and locations of all NEX files pertaining to each individual DBS trial, and each row pertains to the neuronal activity of at least one neuron, in response to one trial of DBS; This table functions as the highest "Parent" in the "Hierarchical Database" model that Im adopting here with matlab tables, with Foreign Keys pertaining to other "Child" tables.';

subjID = 'Kramer';





N.Properties.Description = descr;
N.Properties.VariableDescriptions{1,1} = 'Primary Key';
N.Properties.VariableDescriptions{1,2} = ['Foreign Key for the corresponding row in SortedUnits_', subjID, ' ; relationship of this table to SortedUnits is many-to-one.'];
N.Properties.VariableDescriptions{1,3} = ['Foreign Key for the corresponding row in SweepAnalysisTrials4Paper2018_', subjID, ' ; relationship of this table to SweepAnalysisTrials4Paper2018 is many-to-one.'];
N.Properties.VariableDescriptions{1,4} = 'Unique filename of the .nex file';
N.Properties.VariableDescriptions{1,5} = 'Filepath to the .nex file on your local machine; make sure this is correct.';

NEXprocfiles = N;

open NEXprocfiles



%% Update table metadata for SortedUnits

U = SortedUnits;

descr = 'Table with descriptions of the individual neurons that were recorded on TDT and then spike-sorted on Plexon Offline Sorter; each row pertains to one unique neuron, with no reference to how many DBS trials (or which ones) they were.';


U.Properties.Description = descr;
U.Properties.VariableDescriptions{1,1} = 'Primary Key';
U.Properties.VariableDescriptions{1,2} = 'Date recorded in YYMMDD';
U.Properties.VariableDescriptions{1,3} = 'Channel ID on the recording microelectrode array; if a traditional single channel microelectrode was used, then should just be 1';
U.Properties.VariableDescriptions{1,4} = 'Letter assigned to this unit within Plexon Offline Sorter';
U.Properties.VariableDescriptions{1,5} = 'SU (clear single unit) | MU/noise (spikes were noisy or unable to be reliably called SU) | Mov_spikes (spike only showed up during NHP movement episodes) | MixU (waveforms clearly belong to more than one SU, but unable to reliably cluster them apart)';
U.Properties.VariableDescriptions{1,6} = 'Comments on Spike waveform shape';

SortedUnits = U;



%% Update table metadata for RateBins table


subjID = 'Kramer';
pn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
tabFn = ['RateBins_60sec_1secBins_', subjID, '.mat'];
load([pn, '\', tabFn]);

R = RateBins;

% R(:,2:3) = [];

R.Properties.Description = 'Table with details about the DBS trial rate-binning step, including location of the data file.';

R.Properties.VariableNames{1,2} = 'ratebinFn';
R.Properties.VariableNames{1,3} = 'ratebinPn';

R.Properties.VariableDescriptions{1,1} = 'Primary Key';
R.Properties.VariableDescriptions{1,2} = 'Unique filename of the .mat file holding rate bin data for PRE and DBS epocs for current trial-block.';
R.Properties.VariableDescriptions{1,3} = 'Filepath to the .mat file on your local machine; make sure this is correct.';

RateBins = R;

save([pn, '\', tabFn], 'RateBins');



%% Add RateBins Primary Key to NEXprocfiles as a Foreign Key

clear;
subjID = 'Kramer';

tabPn = 'L:\My Drive\PROJECTS\Thalamic DBS Cx Record\DataProcessing\Tables';
tabFnRate = ['RateBins_60sec_1secBins_', subjID, '.mat'];
tabFnNEX = ['NEXprocfiles_', subjID];
load([tabPn, '\', tabFnRate]);
load([tabPn, '\', tabFnNEX]);

Rt = RateBins(:,1);
Rt.Properties.VariableNames{1,1} = 'RateBins_objectID';

NEXprocfiles = [NEXprocfiles, Rt];
NEXprocfiles.Properties.VariableDescriptions{1,6} = ['Foreign Key for the corresponding row in RateBins_XXsec_YsecBins_', subjID, ' ; relationship of this table to RateBins_XXsec_YsecBins is one-to-one.'];


save([tabPn, '\', tabFnNEX], 'NEXprocfiles');











