%%DBS EEG MSIT
%%Nithya Ramakrishnan 7/29/2019
%%Baylor College of Medicine 
%% LOAD
%%Setting all paths. Just change the path to point to the file that will be
%%in the same directory as this pipeline
addpath('/raid5/rcho/TOOLS/eeglab14_1_1b/')
addpath('/data/rcho/TripolarEEG/SCRIPTS/matlab_tools_bcm/')
addpath(genpath('/data/rcho/TripolarEEG/SCRIPTS/new/'))
addpath('/raid5/rcho/TOOLS/NIC/MATLAB/')
%Channels 1-10 and 23-32 is EEG 
%EOG is channels 12, 14 and 16 
%ECG is channels 14 and 16 
% Put into EEGLab style environment
load('aDBS003_MSIT_2019_06_17_14_19_06_synced_ephys_behav.mat')
for j=1:20;
ephysdata(j,:)=data.EEG{:,j};
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG.setname = 'tripolarcontinuous';
EEG.filename = '';
EEG.filepath = '';
EEG.subject = '';
EEG.group = '';
EEG.condition = '';
EEG.session = [];
EEG.comments = '';
EEG.nbchan = 20;
EEG.trials = 1;
EEG.pnts = length(ephysdata);
EEG.srate = 30000;
EEG.xmin = 0;
EEG.xmax = EEG.pnts*(1/EEG.srate);
EEG.times = linspace(0,EEG.xmax,EEG.pnts);
EEG.data = ephysdata;
EEG.icaact = [];
EEG.icawinv = [];
EEG.icaweights = [];
EEG.icachansind = [];
EEG.icasphere = [];
EEG.chanlocs = [];
EEG.urchanlocs = [];
EEG.chaninfo = [];
EEG.ref = '';
EEG.event = [];
EEG.urevent = [];
EEG.eventdescription = [];
EEG.epoch = [];
EEG.epochdescription = [];
EEG.reject = [];
EEG.stats = [];
EEG.specdata = [];
EEG.splinefile = [];
EEG.icasplinefile = [];
EEG.dipfit = [];
EEG.history = '';
EEG.saved = 'no';
EEG.etc = [];
eeglab redraw
%manually add event markers (need to add the code version of this later)
%File -> import event info -> From Matlab Array or ASCII file
%Make sure to write in input field "latency type duration" and for Number
%of file header lines =1, time unit 1E-3, align event latencies to data
%events=NaN and Auto adjust new events sampling rate= yes 


%% Resampling
EEG = pop_resample(EEG , 1000);
eeglab redraw

%% CREATE CHANLOCS

%for chanlocs, load bvef file
% REF-FCz
% 1-Fp1
% 2-AFz
% 3-Fp2
% 4-F7
% 5-F3
% 6-Fz
% 7-F4
% 8-F8
% 9-FT9
% 10-FC5
% 11-FC1
% 12-FC2
% 13-FC6
% 14-FT10
% 15-T7
% 16-C3
% 17-Cz
% 18-C4
% 19-T8
% 20-CPz
load('/data/rcho/TripolarEEG/DBSstudychanlocs.mat') % path to channel locations
chanlocsT = chanlocs; % electrodes used pt1
chanorder = [71,11,72,44,40,82,41,45,76,52,48,49,53,77,119,14,37,15,120,36]; % cut list to actual electrodes used
chanlocsT = chanlocsT(chanorder); % sort and cut down the chanlocs
EEG.chanlocs = chanlocsT; % add to eeglab
clear chanlocs chanlocsT chanorder
eeglab redraw

%% Filtering

% apply resampling and then filtering and notch filtering
EEG.data=fqfilter(EEG.data,[3,80],1000,'pass',[],5);
EEG.data=fqfilter(EEG.data,[58,62],1000,'N',[],5);
eeglab redraw

% %% RE-REFERENCE
% EEG = pop_reref(EEG, []); % common average reference
% eeglab redraw

%Plot with cleaned data 
for i=1:20
    [pxx1(i,:),f] = pwelch(EEG.data(i,:),5096,[],[4:55],1000);
end
chanlocs=EEG.chanlocs
figure; scmatrix_tripolar(chanlocs, pxx1, [],[],10*log10(f))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove bad channels 

%%%%%%%%%%%%%%%%%%%%%%%%
%manual epoching (also code later)
%Tools -> Extract epochs -> select locking event (12) 
% start and end in seconds -1 2
%baseline correct to -200 ms before stimulus onset 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reject bad epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ARTIFACT REMOVAL via ICA
EEGold = pop_reref(EEG,[]);
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 0);
pop_topoplot(EEG,0, [1:20] ,'tripolarcontinuous resampled',[4 5] ,0,'electrodes','on');
pop_eegplot( EEG, 0, 1, 1);
pop_prop( EEG, 0, [1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20], NaN, {'freqrange' [2 50] });
pop_prop( EEG, 0, [1   14   16], NaN, {'freqrange' [2 50] });

EEG = pop_subcomp( EEG, [1  14  16 ], 0); % add bad components 
eeglab redraw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ready for plotting ERP 