

%% load raw data, this time load all channels

set(0,'DefaultFigureWindowStyle','docked'); % fix matlab's figure positioning bug

% raw data available on
% usb sticks
%
datapath='C:\Users\Jakob Voigts\Downloads\example_ephys\2022-06-12_16-24-58\Record Node 117\experiment2\recording1\continuous\Rhythm_FPGA-111.0\continuous.dat'
D=load_open_ephys_binary(datapath,'continuous',1,'mmap');
%%

data_raw=[];
for ch=[1:4*6] % grab a few channels of raw data from one tetrode
    fname=sprintf('100_%d.continuous',ch)
    [data, timestamps, info]=load_open_ephys_data_faster(fullfile(datapath,fname));
    data_raw(:,end+1) = data;
end;
disp('done loading');
%
data_raw=data_raw.*info.header.bitVolts;
fs = info.header.sampleRate;

%%data_raw=data_raw(1:30000,:); % cut away some data for faster testing

%% plot raw data
clf; hold on;
% data_raw(50000+[0:2],:)=-80; % add fake spike?

plotlim=200000;
spacing=200;
figure(1);
clf;
hold on;

median_data_window = median(data_raw([1:plotlim],:)')';

for i=1:size(data_raw,2)
    plot(data_raw([1:plotlim],i)+(i-1)*spacing,'r');
    plot(data_raw([1:plotlim],i)-median_data_window+(i-1)*spacing,'k');
end;

%% make median for entire recording

Tmax=size(data_raw,1);
median_all_data=zeros(1,Tmax);
stepsize=200000;
for t=1:stepsize:Tmax
    t/Tmax
    w=[1:stepsize]-1+t;
    median_all_data(w) = median(data_raw(w,:)');
end;


