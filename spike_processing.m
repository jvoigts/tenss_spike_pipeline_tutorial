

%% load raw data

% raw data available on
% https://drive.google.com/drive/folders/1CwFcErgp3F3D6I2TB_hTtW1JAQB21TAC?usp=sharing
%
datapath='/home/jvoigts/Dropbox (MIT)/tenss/tenss_2017_lectures_data/Ringo/2017-06-14_12-56-12/'


data_raw=[];
for ch=[1:4]+12 % grab 4 channels of raw data from one tetrode
    fname=sprintf('100_CH%d.continuous',ch)
    [data, timestamps, info]=load_open_ephys_data_faster(fullfile(datapath,fname));
    data_raw(:,end+1) = data;
end;

data_raw=data_raw.*info.header.bitVolts;
fs = info.header.sampleRate;

%data_raw=data_raw(1:30000,:); % cut away some data for faster testing

%% plot

plotlim=100000;
figure(1);
clf;
hold on;
plot(data_raw(1:plotlim,:));


%% filter

clf; hold on;
[b,a] = butter(3, [300 3000]/(fs/2)); % choose filter (normalize bp freq. to nyquist freq.)

data_bp=filter(b,a,data_raw); %apply filter in one direction

plot(data_bp(1:plotlim,:));

% find treshold crossings
treshold=20;
crossed= min(data_bp,[],2)<-treshold; % trigger if _any_ channel crosses in neg. direction

spike_onsets=find(diff(crossed)==1);

length_sec=size(data,1)/fs;
fprintf('got %d candidate events in %dmin of data, ~%.2f Hz\n',numel(spike_onsets),round(length_sec/60),numel(spike_onsets)/length_sec)

for i=1:numel(spike_onsets)
    if(spike_onsets(i)<plotlim)
        plot([1 1].*spike_onsets(i),[-1 1].*treshold*2,'k--')
    end;
end;


%% extract spike waveforms and make some features

spike_window=[1:32]-5; % grab some pre-treshold crossign samples

spikes=[];
spikes.waveforms=zeros(numel(spike_onsets),4*numel(spike_window)); % pre-allocate memory
spikes.peakamps=zeros(numel(spike_onsets),4);


for i=1:numel(spike_onsets)
    this_spike=(data_bp(spike_onsets(i)+spike_window,:));
    
    spikes.waveforms(i,:)= this_spike(:);% grab entire waveform
    spikes.peakamp(i,:)=min(this_spike); % grab 4 peak amplitudes
end;


%% plot peak to peak amplitudes
clf; hold on;
plot(spikes.peakamp(:,2),spikes.peakamp(:,4),'.');
daspect([1 1 1]);