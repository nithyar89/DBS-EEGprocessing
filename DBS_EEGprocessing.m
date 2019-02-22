%%DBS EEG initial analysis 
%%Nithya Ramakrishnan 1/29/2018
%%Baylor College of Medicine 
%% LOAD
%%Setting all paths. Just change the path to point to the file that will be
%%in the same directory as this pipeline
%cd ..BROWN_share/
addpath('/raid5/rcho/TripolarEEG/BROWN_share/eeglab14_1_1b/')
addpath(genpath('/raid5/rcho/TripolarEEG/BROWN_share/analysis-tools-master'))

%%%add file pattern for the file that you want load 
[contents] = get_directory_contents('/data/rcho/TripolarEEG/DBS_study/DBS002/2018-10-17_13-57-27', '100', 'CH');
contents = natsortfiles(contents); % sort to natural order
pathway= '/data/rcho/TripolarEEG/DBS_study/DBS002/2018-10-17_13-57-27/';

%%%Loads the data
for iter = 1:32
    a = [pathway,contents{iter}];
    if iter == 1
        [data, timestamps, ~] = load_open_ephys_data(a);
        ephysdata = zeros(32,length(data));
        ephysdata(iter,:) = data';
        %         clear data
    else
        [data, ~, ~] = load_open_ephys_data(a);
        ephysdata(iter,:) = data';
    end
end
%Channels 1-10 and 23-32 is EEG 
%EOG is channels 12, 14 and 16 
%ECG is channels 14 and 16 
% Put into EEGLab style environment
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG.setname = 'tripolarcontinuous';
EEG.filename = '';
EEG.filepath = '';
EEG.subject = '';
EEG.group = '';
EEG.condition = '';
EEG.session = [];
EEG.comments = '';
EEG.nbchan = 32;
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
%% Resampling to 1000Hz
EEG = pop_resample(EEG , 1000);
eeglab redraw
%% DETERMINE USELESS CHANNELS
%defing channels not recording brain activity
uc = 11:22;
EEG.data = EEG.data([1:uc(1)-1,uc(end)+1:end],:);
EEG.nbchan = size(EEG.data,1);
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
% Band Pass Filter - 0.1-100Hz and notch filtering to remove line noise
EEG.data=fqfilter(EEG.data,[0.1,100],1000,'pass',[],5);
EEG.data=fqfilter(EEG.data,[58,62],1000,'N',[],5);
eeglab redraw

%% RE-REFERENCE
%average reference
EEG = pop_reref(EEG, []); % common average reference
eeglab redraw

%% SPECTRA
%plotting the psd 
for i=1:20
    [pxx(i,:),f] = pwelch(EEG.data(i,:),5096,[],[1:55],1000);
end
for i=1:20
    if i==1
        figure;plot(f,10.*log10(pxx(i,:)));hold on;
    else
        plot(f,10.*log10(pxx));
    end
end
%%For scalp plot
chanlocs=EEG.chanlocs;
figure; scmatrix_tripolar(chanlocs, pxx,[],[],10*log10(f))



