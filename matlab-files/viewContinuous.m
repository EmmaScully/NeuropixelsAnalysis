% Rat47 - 13 Nov. 

fold1 = '\\ad.monash.edu\home\User009\escu0001\Documents\SUMMER RESEARCH SCHOLARSHIP 2020_2021\neuropixels\neuropixels data\2020-11-13_19-10-06\Record Node 106\';
fold2 = 'experiment2\recording1\'; %exp1, 2 or 3, always recording1
analysis_tools_folder = '\\ad.monash.edu\home\User009\escu0001\Documents\GitHub\analysis-tools';
%addpath(analysis_tools_folder); just run this in command window and
%savepath upon first run.

jsonFile = [fold1 fold2 'structure.oebin'];
sampling_rate = 30000;

if sampling_rate ==30000
   sr = 1;
elseif sampling_rate == 2500
   sr = 2;
end

D = load_open_ephys_binary(jsonFile, 'continuous', sr); % 1 - sampled at 30 kHz; 2 - sampled at 2.5 kHz

nCh = size(D.Data,1);

%% low-pass filter, downsample to 1 kHz median shift
% As R=30, decimate suggests doing this twice (e.g. with R=5 and R=6)
% Using default filter, cutoff is 320 Hz (0.8*0.8*(30000/2)/30)
nCh = size(D.Data,1);
lfp = nan(nCh,ceil(size(D.Data,2)/30));

for a = 1:nCh
    tmp = decimate(D.Data(a,:),6);
     lfp(a,:) = decimate(tmp,5);
end

%% Additional filter to <100 Hz
cutoff = 100;
Wn = cutoff/(1000/2); % Normalized cutoff frequency
[BB,AA] = butter(8,Wn,'low');
for a = 1:nCh
    tmp = filter(BB,AA,lfp(a,:));
    lfp(a,:) = tmp - median(tmp); % median shift
end
    
    
lfp_med = median(lfp);

% seems like every 826 samples there is a peak.
magicPer = 826; % magic number of samples per period (totally guessed)
tt = 1:magicPer;
ind = 1:magicPer:length(lfp_med);
ind2 = ind(1:end-1)' + tt;

for a = 1:nCh
    tmp = lfp(a,:);
    lfp_align(a,:) = mean(tmp(ind2));
end
lfp_med_align = mean(lfp_med(ind2));


%% Alternative approach with raw data
raw = zeros(size(D.Data));
for a = 1:nCh
    raw(a,:) = D.Data(a,:)-median(D.Data(a,:));
end
raw_med = median(raw);
magicPer2 = 2064; %24773/12;
tt2 = 1:magicPer2;
ind = 1:magicPer2:length(raw);
ind2 = ind(1:end-1)' + tt2;
raw_align = zeros(nCh,magicPer2);
for a = 1:nCh
    tmp = raw(a,:);
    raw_align(a,:) = sgolayfilt(mean(tmp(ind2)),3,25); %was 301
end
% for filter - 100 Hz = 10 ms period = 300 samples
% raw_align_med = mean(raw_med(ind2));
% sg = sgolayfilt(raw_align_med,5,111);

figure
plot(tt2/2500, mean(raw_align))
xlabel('s')

figure
nSkip = 15;
plot(tt2/2500, raw_align(1:nSkip:nCh,tt2)+(1:nSkip:nCh)'*5)

%%
xlabel('s')

%% PLOTTING
figure
plot(lfp_med)
diff(find(lfp_med>200))

figure
plot(lfp_med_align)
title('AllCh Median LFP aligned on a guess')
ylabel('uV')
xlabel('Time (ms)')

figure
plot(tt,lfp_align(1:50:nCh,tt))