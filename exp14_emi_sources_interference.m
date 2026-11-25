%% Exp 14 : Simulation of EMI Sources and Interference Analysis using MATLAB
% Needs: MATLAB + Signal Processing Toolbox (butter, filtfilt)
% Aim  : model 1 kHz signal + 3 EMI sources, FFT analysis, band-pass mitigation, SNR improvement
clc; clear; close all;
rng(1);                                   % repeatable noise (remove for new noise each run)

fs = 10000; t = 0:1/fs:0.1-1/fs;          % sampling frequency, 1000 samples
signal    = sin(2*pi*1000*t);             % desired 1 kHz signal
emi_sine  = 0.5*sin(2*pi*2500*t);         % sinusoidal EMI (switching supply)
emi_noise = 0.3*randn(size(t));           % broadband EMI (sigma = 0.3)
emi_imp   = zeros(size(t)); emi_imp(1:20:end) = 1;   % impulsive EMI (every 20 samples)
emi_total = emi_sine + emi_noise + emi_imp;
received  = signal + emi_total;

%% Input SNR : simulated and theoretical
SNR_in  = 10*log10(mean(signal.^2)/mean(emi_total.^2));
Ps = 0.5; Pe = 0.5^2/2 + 0.3^2 + 1/20;    % signal power, EMI power (tone + noise + impulses)
SNR_th  = 10*log10(Ps/Pe);
fprintf('Input SNR  simulated = %.2f dB, theoretical = %.2f dB, difference = %.2f dB (limit 1.5 dB)\n', ...
        SNR_in, SNR_th, abs(SNR_in-SNR_th));

%% Time-domain plots
figure('Name','Time domain');
subplot(3,1,1); plot(t,signal);    xlim([0 0.01]); grid on; title('Desired Signal (1 kHz)');
subplot(3,1,2); plot(t,emi_total); xlim([0 0.01]); grid on; title('Combined EMI');
subplot(3,1,3); plot(t,received);  xlim([0 0.01]); grid on; title('Received Signal with EMI'); xlabel('Time (s)');

%% Frequency domain (FFT)
N = length(received); f = (-N/2:N/2-1)*(fs/N);
R = abs(fftshift(fft(received)))/N;
figure('Name','Spectrum'); plot(f,R); xlim([0 5000]); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude'); title('Frequency Spectrum of Received Signal');

pos = f >= 0; fp = f(pos); Rp = R(pos);
[~,i1] = min(abs(fp-1000));  win = fp>2300 & fp<2700;
[pk,ik] = max(Rp(win)); fwin = fp(win); f_emi = fwin(ik);
fprintf('1 kHz signal peak   : %.3f at %.0f Hz\n', Rp(i1), fp(i1));
fprintf('EMI tone peak       : %.3f at %.1f Hz (frequency error %.2f %%, limit 1 %%)\n', pk, f_emi, 100*abs(f_emi-2500)/2500);
noise_band = fp>3300 & fp<4900 & mod(round(fp),500)~=0;      % avoid impulse-train harmonics (multiples of fs/20 = 500 Hz)
fprintf('Noise floor         : %.1f dB (re 1)\n', 20*log10(mean(Rp(noise_band))));
fprintf('Impulse harmonics   : at multiples of fs/20 = %d Hz (2.5 kHz tone coincides with the 5th)\n', fs/20);

%% EMI mitigation: 4th-order Butterworth band-pass (700-1300 Hz), zero phase
[b,a] = butter(4,[700 1300]/(fs/2),'bandpass');
cleaned = filtfilt(b,a,received);
SNR_out = 10*log10(mean(signal.^2)/mean((cleaned-signal).^2));
Improvement = SNR_out - SNR_in;
amp_out = max(abs(cleaned(100:end-100)));                     % ignore filter edge effects
fprintf('Output SNR = %.2f dB, improvement = %.2f dB (target >= 10 dB)\n', SNR_out, Improvement);
fprintf('1 kHz amplitude after filter = %.3f (must be 0.9-1.1)\n', amp_out);

figure('Name','Filtering');
subplot(2,1,1); plot(t,received,t,cleaned,'LineWidth',1.2); xlim([0 0.01]); grid on;
legend('Received','Filtered'); title('Before and after EMI mitigation');
subplot(2,1,2); [H,w] = freqz(b,a,2048,fs); plot(w,20*log10(abs(H))); grid on; ylim([-80 5]);
xlabel('Frequency (Hz)'); ylabel('|H| (dB)'); title('Band-pass filter response');
