clear all; close all; clc;
fold1 = '\\ad.monash.edu\home\User009\escu0001\Documents\SUMMER RESEARCH SCHOLARSHIP 2020_2021\neuropixels\neuropixels data\2020-11-13_19-10-06\Record Node 106\';
fold2 = 'experiment2\recording1\'; %exp1, 2 or 3, always recording1
analysis_tools_folder = '\\ad.monash.edu\home\User009\escu0001\Documents\GitHub\analysis-tools';
addpath(analysis_tools_folder); just run this in command window and
savepath upon first run.

jsonFile = [fold1 fold2 'structure.oebin'];
sampling_rate = 30000;

if sampling_rate ==30000
   sr = 1;
   min_peak_distance = 24000;
elseif sampling_rate == 2500
   sr = 2;
   min_peak_distance = 24000/12;

end

D = load_open_ephys_binary(jsonFile, 'continuous', sr); % 1 - sampled at 30 kHz; 2 - sampled at 2.5 kHz

nCh = size(D.Data,1);

% plot signal in frequency domain 
        data = D.Data;
        fs = sampling_rate;
        ts = 1/fs;
        t = 0:ts:1-ts;

        nfft = length(data);
        nfft2 = 2.^nextpow2(nfft);

        fy = fft(data, nfft2);
        fy = fy(1:nfft2/2);

        xfft = fs.*(0:nfft2/2-1)/nfft2;
        subplot(2,1,1)
        plot(t,data(1:30000));
        subplot(2,1,2)
        plot(xfft, abs(fy/max(fy)));
% FIRST LOW-PASS ANTI ALIASING  
for c = 1:10
filt1(c,:) = lowpass(data(c, :), 10000, sampling_rate);
end
% check frequency response
fy2 = fft(filt1, nfft2);
fy2 = fy2(1:nfft2/2);

xfft2 = fs.*(0:nfft2/2-1)/nfft2;
figure
subplot(2,1,1)
plot(t,data(1:30000));
subplot(2,1,2)
plot(xfft2, abs(fy2/max(fy2)));

% down sample after anti - aliasing
for c = 1:10
downsample1(c,:) = downsample(filt1(c,:), 30);
end
t_down = downsample(t, 30);
% frequency analysis again
fs_down = 1e3;
nfft = length(downsample1);
nfft2 = 2.^nextpow2(nfft);

for c = 1:10
fy3(c,:) = fft(downsample1(c,:), nfft2);
fy3(c,:) = fy3(c,1:nfft2/2);
end
xfft2 = fs_down.*(0:nfft2/2-1)/nfft2;
figure
subplot(2,1,1)
plot(t_down,downsample1(1:1000));
subplot(2,1,2)
plot(xfft2, abs(fy3/max(fy3)));

% notch-filter the 50Hz noise
nyquist = fs_down/2;
Wn = [40/nyquist, 60/nyquist];
[b,a] = butter(3, Wn, 'stop');
notched = filter(b,a,downsample1);

check frequency again
fy4 = fft(notched, nfft2);
fy4 = fy4(1:nfft2/2);

xfft2 = fs_down.*(0:nfft2/2-1)/nfft2;
figure
subplot(2,1,1)
plot(t_down,notched(1:1000));
subplot(2,1,2)
plot(xfft2, abs(fy4/max(fy4)));

% find peaks in 10 channels
pre = 0.1;
post = 0.5;
for c = 1:10
    [peaks, indices] = findpeaks(-D.Data(c,:), 'MinPeakDistance', min_peak_distance);
        for n = 2:length(peaks)-1
            data_peaks(c,n,:) = D.Data(c,indices(n)-pre*sampling_rate:indices(n)+post*sampling_rate);  
        end
end


