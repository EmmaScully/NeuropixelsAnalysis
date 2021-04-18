clear all; close all; clc;

f1 = 30e3*2*pi;
f2 = 20e3*2*pi;
f3 = 10e3*2*pi;
f4 = 1e3*2*pi;
f5 = 100*2*pi;

t = 0:0.00001:0.01;
x1 = sin(f1*t);
x2 = sin(f2*t);
x3 = sin(f3*t);
x4 = sin(f4*t);
x5 = sin(f5*t);

plot(t, x1);
hold on
plot(t, x2);
plot(t, x3);
plot(t, x4);
plot(t, x5);

combined = x1+x2+x3+x4+x5;
figure
plot(t, combined);
%% ffts 


y_combined = fft(combined);
N = length(combined);
fs = 30e3;
fnyquist = fs/2;
bin_vals = [0:N-1];
fax_hz = bin_vals*fs/N;
N_2 = ceil(N/2);
plot(fax_hz(1:N_2), abs(y_combined(1:N_2)));

%% 
X_mags = abs(fftshift(fft(combined)));
bin_vals = 0 : (N-1);
N_2 = ceil(N/2);
fax_Hz = (bin_vals-N_2)*fs/N;
plot(fax_Hz, X_mags)
xlabel('Frequency (Hz)')
ylabel('Magnitude');
title('Double-sided Magnitude spectrum (Hertz)');
axis tight