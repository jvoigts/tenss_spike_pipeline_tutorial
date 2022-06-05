



fname='C:\Users\voigtsj\Downloads\intan_2019-11-13T07_51_08.raw'
fid_in = fopen(fname);

data_size=30000*60*30;
offs_samples=30000*200;
fseek(fid_in, 2*nchan*offs_samples, 'bof');

v_temp = double(fread(fid_in, [nchan,data_size], 'uint16'))-2^16/2; %(fid, [],  'uint16', 64); % read all chanels per chunk
ref=(repmat( median(v_temp(:,:)) ,[nchan 1]));
v_temp=v_temp-ref;
v_temp=v_temp(use_chans,:);
t=whos('v_temp');
t.bytes./(1024^2)
%plot(v_temp')

save('ancalogon_2019-11-13T07_51_08_TT13_30min.mat','v_temp');