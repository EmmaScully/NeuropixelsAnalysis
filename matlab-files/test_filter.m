clear all;
close all;
clc;

f1 = 10e3*2*pi;
f5 = 100*2*pi;

fs = 30000;
ts = 1/fs;
t = 0:ts:0.1-ts;

x1 = sin(f1*t);
x5 = sin(f5*t);


x = x1+x5;
plot(t, x);

nfft = length(x);
nfft2 = 2.^nextpow2(nfft);

fy = fft(x, nfft2);
fy = fy(1:nfft2/2);

xfft = fs.*(0:nfft2/2-1)/nfft2;
subplot(2,1,1)
plot(t,x);
subplot(2,1,2)
plot(xfft, abs(fy/max(fy)));

%% low-pass filtering

low = lowpass(x, 1000, fs);
figure
subplot(2,1,1);
plot(t, low);

fy_low = fft(low, nfft2);
fy_low = fy_low(1:nfft2/2);

xfft_low = fs.*(0:nfft2/2-1)/nfft2;
subplot(2,1,2)
plot(xfft_low, abs(fy_low/max(fy_low)));



%% high-pass filtering

high = highpass(x, 1000, fs);
figure
subplot(2,1,1);
plot(t, high);

fy_high = fft(high, nfft2);
fy_high = fy_high(1:nfft2/2);

xfft_high = fs.*(0:nfft2/2-1)/nfft2;
subplot(2,1,2)
plot(xfft_high, abs(fy_high/max(fy_high)));