
set(0,'DefaultFigureWindowStyle','docked'); % fix matlab's figure positioning bug

%% or load from MAT file
%{
%load('C:\Users\Jakob Voigts\Downloads\example_ephys\ancalogon_2019-11-13T07_51_08_TT13_100sec.mat');
load('C:\Users\Jakob Voigts\Downloads\example_ephys\ancalogon_2019-11-13T07_51_08_TT13_600sec.mat');
fs = 30000; % 30KHz
v_temp=v_temp';
%}


%% load raw data

% raw data available on
% usb sticks
%
addpath('C:\Users\Jakob Voigts\Documents\GitHub\analysis-tools\');
%datapath='C:\Users\Jakob Voigts\Desktop\pondababa_2021-09-18_09-33-02\Record Node 111\'
datapath='C:\Users\Jakob Voigts\Downloads\example_ephys\2022-06-12_16-24-58\Record Node 117\experiment1\recording1\structure.oebin'
datapath='C:\Users\Jakob Voigts\Desktop\ephys data spike gals\mouse_task\ephys\2022-06-15_15-32-19\Record Node 104\experiment2\recording1\structure.oebin'

D=load_open_ephys_binary(datapath,'continuous',1,'mmap');



tt=3;
tt=1;

chunksize=50; % seconds
treshold=40;

fs = D.Header.sample_rate;
chunksize_samples=D.Header.sample_rate*chunksize;
chunk_boundaries=unique([[0:chunksize_samples:numel(D.Timestamps)],numel(D.Timestamps)]);
pad = 1*D.Header.sample_rate; % padding each chunk in sec

% output variables for accumulation
spike_times=[];
spike_waveforms=[];
spike_amps=[];

for chunk=1:numel(chunk_boundaries)-1
    fprintf('processing  TT %d, chunk %d/%d, ',tt,chunk,numel(chunk_boundaries)-1);
    chunk_padded=min(max(1,[chunk_boundaries(chunk)-pad , chunk_boundaries(chunk+1)+pad]),numel(D.Timestamps));

    v_t=D.Data.Data.mapped(:,chunk_padded(1):chunk_padded(2))'; % grab time slice
    v_t=v_t-repmat(median(v_t')',1,size(v_t,2)); % median remove

    v_trode=v_t(:,[1:4]+((tt-1)*4))';

    v_trode=v_trode.*0.195; % go from raw to uV


    %% filter

    clf; hold on;
    [b,a] = butter(3, [300 3000]/(fs/2)); % choose filter (normalize bp freq. to nyquist freq.)

    data_bp=filtfilt(b,a,double(v_trode')); %apply filter in one direction

    plot(data_bp); drawnow;

    %% find treshold crossings

    crossed= min(data_bp,[],2)<-treshold; % trigger if _any_ channel crosses in neg. direction

    spike_onsets=find(diff(crossed)==1);
    spike_onsets(spike_onsets<pad)=[]; spike_onsets(spike_onsets>pad+chunksize_samples)=[]; % remove doubled ones

    length_sec=size(v_trode,1)/fs;
    fprintf('got %d candidate events in %dmin of data, ~%.2f Hz\n',numel(spike_onsets),round(chunksize/60),numel(spike_onsets)/chunksize);


    %% extract spike waveforms and make some features

    spike_window=[1:32]-5; % grab some pre-treshold crossign samples

    spike_onsets(spike_onsets<10)=[];


    spike_times=[spike_times;spike_onsets+chunk_boundaries(chunk)];

    tmp_amps=zeros(numel(spike_onsets),4);


    tmp_waveforms=zeros(numel(spike_onsets),4*numel(spike_window)); % pre-allocate memory
    for i=1:numel(spike_onsets)
        this_spike=(data_bp(spike_onsets(i)+spike_window,:));

        tmp_waveforms(i,:)= this_spike(:);% grab entire waveform
        tmp_amps(i,:)=min(this_spike);
    end;
    spike_waveforms=[spike_waveforms;tmp_waveforms];
    spike_amps=[spike_amps,tmp_amps'];

end
spike_amps=spike_amps';

%% plot peak to peak amplitudes
clf; hold on;
subs=10;
plot(spike_amps(1:subs:end,3),spike_amps(1:subs:end,4),'.');
%daspect([1 1 1]);

%% save for simpleclust
disp('saving')
mua=[];
mua.ts_spike= ([1:(32*4)]-5)./D.Header.sample_rate;
mua.ncontacts=4;
mua.waveforms = spike_waveforms;
mua.ts = double(spike_times./D.Header.sample_rate);
mua.Nspikes = size(spike_waveforms,1);
mua.sourcechannel = tt;
save(sprintf('extracted_spikes_TT%d.mat',tt), 'mua');
disp('saved')