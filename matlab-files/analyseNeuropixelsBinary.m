% Rat47 - 13 Nov. 

fold1 = '\\ad.monash.edu\home\User009\escu0001\Documents\SUMMER RESEARCH SCHOLARSHIP 2020_2021\neuropixels\neuropixels data\2020-11-13_19-10-06\Record Node 106\';
fold2 = 'experiment2\recording1\'; %exp1, 2 or 3, always recording1
analysis_tools_folder = '\\ad.monash.edu\home\User009\escu0001\Documents\GitHub\analysis-tools';
%addpath(analysis_tools_folder); just run this in command window and
%savepath upon first run.

jsonFile = [fold1 fold2 'structure.oebin'];
sampling_rate = 2500;

if sampling_rate == 30000
   sr = 1;
elseif sampling_rate == 2500
   sr = 2;
end

D = load_open_ephys_binary(jsonFile, 'continuous', sr); % 1 - sampled at 30 kHz; 2 - sampled at 2.5 kHz

nCh = size(D.Data,1);
%% finding artefacts 
integers = zeros(length(D.Data),1);
for a=1:length(D.Data)
if D.Data(1,a) < -4000
    integers(a) = a;
end
end
diff_integers = diff(find(integers));
diff_integers = diff_integers(diff_integers>1);
rate_of_flash = floor(mean(diff_integers));
maxflash = max(diff_integers);
%% average across timestamps
indx = 1:rate_of_flash:length(D.Data);

test = average_flashes(3, maxflash, D.Data, -4000);

%% aligning on first flash
processed = zeros(size(D.Data));
for c = 1:nCh
    for a = 1:length(D.Data)
        if D.Data(c,a) < -4000
         processed(c,1:length(D.Data)-a+1) = D.Data(c,a:end); 
         break
        end
    end
end

%% plot pre processed 
plot(processed(1,:));
%% testing function
number_of_flashes = ceil(length(processed)/maxflash);
blank_matrix = zeros(number_of_flashes,maxflash);
average_matrix = zeros(4, maxflash);
for c = 1:1
    e = 1;
        for b = 1:number_of_flashes
            for d = 1:maxflash %from 1 - 24776
                blank_matrix(b,d) = processed(c, e);
                e = e+1
            end
        end 
    
end
    average_matrix(c,:) = mean(blank_matrix);



%% Spectral analysis
x = D.Data;
n = length(x);
fs = 30000;
dt = 1/fs;
t = (0:n-1)/fs;
y = fft(x);
f = (0:n-1)*(fs/n);
power = abs(y).^2/n;


plot(f, power);

%% median filtering to remove artefacts 
x_median = medfilt1(x(1,:),10);
figure;
plot(t, x_median);

%% Band-pass between 300 - 3000 Hz
bp = zeros(size(D.Data));
for a = 1:nCh
bp(a,:) = bandpass(D.Data(a,:), [300 3000], 30000);
end


%% Down-sampling 
%Dec = zeros(size(D.Data));
for a = 1:nCh
temp(a,:) = decimate(bp(a,:), 6);
bp_decimate(a,:) = decimate(temp(a,:), 2); %the signal is now at 2.5kHz
end

%% plot the band passed signal and original
%figure
%plot(D.Data);
%hold on
figure
plot(bp_decimate(1,:));
%% bandpass again
bp2 = zeros(size(bp_decimate));
for a = 1:nCh
bp2(a,:) = bandpass(bp_decimate(a,:), [300, 3000], 2500);
end

%% plot the band passed signal and original
%figure
%plot(D.Data);
%hold on
hold on
plot(bp2(1,:));
%% 
magicPer2 = 24773; %24773/12;
tt2 = 1:magicPer2;
ind = 1:magicPer2:length(D.Data);
ind2 = ind(1:end-1)' + tt2;


for a = 1:nCh
    tmp = D.Data(a,:);
    Dec_align(a,:) = sgolayfilt(tmp,3,25); %was 301, removes noise 
end

%% 
figure
plot(1:length(Dec_align), mean(Dec_align))
xlabel('s')
hold on
plot(tt2, D.Data(1,:));
%% Plotting original vs down sampled
figure
t = 1:length(Dec);
plot(t, Dec(1,:));

figure
t = 1:length(D.Data);
plot(t, D.Data(1,:));


%% Alternative approach with raw data
raw = zeros(size(D.Data));
for a = 1:nCh
    raw(a,:) = D.Data(a,:) -median(D.Data(a,:));
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


%% Plots
figure
plot(tt2/30000, mean(raw_align))
xlabel('s')

figure
nSkip = 15;
plot(tt2/30000, raw_align(1:nSkip:nCh,tt2)+(1:nSkip:nCh)'*5)
xlabel('s')


% %% PLOTTING
% figure
% plot(lfp_med)
% diff(find(lfp_med>200))
% 
% figure
% plot(lfp_med_align)
% title('AllCh Median LFP aligned on a guess')
% ylabel('uV')
% xlabel('Time (ms)')
% 
% figure
% plot(tt,lfp_align(1:50:nCh,tt))