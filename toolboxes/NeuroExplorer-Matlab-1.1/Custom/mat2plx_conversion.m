%% Run Filippo's code (emailed on 03/26/15) on acquired data

[fn,pn]=uigetfile('mat','Pick file','C:\data\kiwi\Anod repeat');
% pn = 'D:\Data\Kramer\New folder\'; % pathname of data
% fn = '13_04_19k_001.mat'; % specific filename to convert
load([pn fn]);

% Plot Spiking
d=F.d_sub_C3;%define the time series to be converted from .mat to .plx
fig=figure;
plot([1:length(d)]/F.s_rate,F.d_sub_C3);%initial plot view; sanity check 
spike_amplitude = 150;
ylim([-spike_amplitude,spike_amplitude])

%% Convert
n=[];
dsn=d*(254*255/(2*spike_amplitude));% % rescale subtracted data to better use 16 bit encoding in binary file
n(1:2:2*length(dsn)-1)=floor(double(dsn)/255); % first 8 bit of every sample (big part of the value)
n(2:2:2*length(dsn))=rem(double(dsn),255); % second 8 bit of every sample (little endian coding)
fid=fopen([pn,fn(1:end-4),sprintf('.2plx')],'w+'); % save binary file output to plexon
fwrite(fid,n,'int8');
fclose(fid);



%% Once sorted on Plexon Offline Spike Sorter, get timestamps
addpath('I:\NRTL\tACS-SI\Joe_code_04_02_2015')
nexFilename = '15_05_29u_0020.nex';
nexFile = readNexFile([pn nexFilename]); % nex file is presumably in the same directory as the original data
recording_filename = strsplit(nexFilename, '.'); recording_filename = recording_filename{1};
N_NEURONS = length(nexFile.neurons);

% Correct times based on "Filippo Ratio" (i.e. 22.321428298950195/22 )
frequency_correction_ratio = 22.321428298950195 / 22;
Fs = 2.75e3 * frequency_correction_ratio;
CAI_001_TimeBegin_corrected = CAI_001_TimeBegin / frequency_correction_ratio;
CAI_001_TimeEnd_corrected = CAI_001_TimeEnd / frequency_correction_ratio;

neuron_timestamps = cell(N_NEURONS,1);
neuron_labels = cell(N_NEURONS,1);
neuron_label_base = [recording_filename '_unit_'];
for i=1:N_NEURONS
    neuron_timestamps{i} = nexFile.neurons{i}.timestamps / frequency_correction_ratio;
    neuron_labels{i} = strrep([neuron_label_base sprintf('%02d',i)],'_','\_');
end
