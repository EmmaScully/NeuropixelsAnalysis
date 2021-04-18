%% try to extract spikes: one channel version
clear all; close all; clc;

fold1 = 'D:\UNI\SUMMER RESEARCH SCHOLARSHIP\NEUROPIXELS DATA\2020-11-13_19-10-06\Record Node 106\';
fold2 = 'experiment2\recording1\'; %exp1, 2 or 3, always recording1
analysis_tools_folder = 'D:\UNI\SUMMER RESEARCH SCHOLARSHIP\analysis-tools';
%% UNCOMMENT THIS SECTION ON FIRST RUN
% addpath(analysis_tools_folder);
% savepath;
%% 
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

%% variables
time_period = 2; % seconds
samples = sampling_rate*time_period;
nCh = size(D.Data,1);
t = linspace(0,time_period,samples);
plot(t,D.Data(1,1:samples))
%% Step 1: Find standard deviation of the signal
SD = std(D.Data(1,:));
threshold = -4*SD;
hold on
plot([1,time_period], [threshold, threshold]);

%% Step 2: Find parts of signal that cross this threshold
% find aeras below threshold
indices = D.Data(1, :)<threshold;
locater = indices.*D.Data(1, :);
%indices = find(indices);
%locater = D.Data(1, indices);

%extract peaks from these areas
[a, peak_locals] = findpeaks(-locater, 'MinPeakDistance', 2);
a = -a;
peaks_in_original = D.Data(1,indices(peak_locals));

%% convert back to time
indices_in_original = indices(peak_locals);
convert_to_seconds = indices_in_original/sampling_rate;
hold on
plot(convert_to_seconds, peaks_in_original, 'o');
% for n = 1:length(peak_locals)
%     spikes(n, :) = (
% end