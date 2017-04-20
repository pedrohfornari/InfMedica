%% Informatica medica trabalho experimental 1
 % Pedro Henrique kappler Fornari
 % 13104320
 % Manipulacoes e filtragens em sinal de ECG e EEG
clear all;
close all;

%% ECG details
fs = 128;
gain = 200;

%% Select database to open
d = what('sinais_ecg');
[s,v] = listdlg('ListString','Select a file:',...
                'SelectionMode','single',...
                'ListString',d.mat);
% Import data
ecg = load(d.mat{s});

%% Select ecg and scale it in mv
ecg = ecg.val(1,:)/gain;
 
% create time vector
tecg = (0:1/fs:(size(ecg, 2)-1)/fs);           
 
% create random variable to plot 5 seconds
x = round(1000000*rand());
 %plot
figure('Name', 'ECG raw signal');
subplot(2, 1, 1);
plot(tecg, ecg);
title('raw complete ECG');
subplot(2, 1, 2);
hold on
plot(tecg(x:(x+(5*fs))), ecg(x:(x+(5*fs))));
axis([x/fs ((x/fs) + 5) min(ecg) max(ecg)]);
grid on
legend('ECG Signal');
xlabel('Seconds');
ylabel('Voltage(mV)');
title('raw ECG 5 seconds');
hold off
 %% Nomalize signal
%maxecg = max(ecg);
%minecg = min(ecg);
%ecg_norm = (ecg-minecg)/(maxecg-minecg);
ecg_norm = mat2gray(ecg);
%  
figure('Name', 'Normalized ECG')
plot(tecg(x:(x+(5*fs))), ecg(x:(x+(5*fs))));
hold on
plot(tecg(x:(x+(5*fs))), ecg_norm(x:(x+(5*fs))));
axis([x/fs ((x/fs) + 5) -1.25 3.5]);
grid on
legend('ECG Signal','ECG normalized Signal')
xlabel('Seconds')
ylabel('Voltage(mV)')
title('raw ECG normalized 5 seconds')
hold off
 
%% Resample signal
%ups = 1024;
%dns = fs;
%ecgr = resample(ecg, ups, dns);
newFs = 1000;

% create resampled time vector
rtecg = 0:(1/newFs):(size(ecg,2)-1)/fs;
%resample by spline
recg = spline(tecg, ecg, rtecg);
%  
%  % plot resampled signal
figure('Name', 'ECG after and before interpolation');
subplot(2, 1, 1)
plot(tecg(x:(x+(5*fs))), ecg(x:(x+(5*fs))));
axis([x/fs (x/fs + 5) min(ecg) max(ecg)]);
grid on
xlabel('Seconds');
ylabel('Voltage(mV)');
title('ECG to 128Hz');
subplot(2, 1, 2)
hold on
plot(rtecg(x:(x+(5*newFs))), recg(x:(x+(5*newFs))));
axis([x/newFs (x/newFs + 5) min(recg) max(recg)]);
grid on
legend('ECG resampled signal')
xlabel('Seconds');
ylabel('Voltage(mV)');
title('Ressampled ECG to 1kHz');
hold off

%% Design Notch filter for 60 Hz
wo = 60/(newFs/2);
bw = wo/35;
[b60Hz, a60Hz]  = iirnotch(wo, bw);
fvtool(b60Hz, a60Hz,'Color','white');
 
 %% Design High pass filter to exclude base line
Fstop = 1e-4;
Fpass = 1e-2;
Astop = 100;
Apass = 0.15;

[b0Hz, a0Hz] = butter(5, 1/newFs, 'high'); 
fvtool(b0Hz, a0Hz,'Color','white');

%% Add some frequencies to the signal so we can test the filters
recgTest = recg + sin(2*pi*60*rtecg) + 5;
 
% plot spectrums
Nfft = 1e5;
freq = (newFs/2*linspace(0,1,Nfft/2+1));

figure ('Name', 'Dirty ECG spectrum')% Dirty signal
ecgfft = (1/length(recgTest))*fft(recgTest, Nfft);
plot(freq, 2*abs(ecgfft(1:Nfft/2+1)));
 
recgTest = filtfilt(b60Hz, a60Hz, recgTest);
recgTest = filtfilt(b0Hz, a0Hz, recgTest);

figure ('Name', 'ECG spectrum without 60Hz and base line')% Clean up 60Hz interference and remove base line
ecgfft = (1/length(recgTest))*fft(recgTest, Nfft);
plot(freq, 2*abs(ecgfft(1:Nfft/2+1))); 

 %% Filter the hole ressampled ecg
 figure('Name', 'ECG filtered') %real signal filtered
 recg = filtfilt(b60Hz, a60Hz, recg);
 recg = filtfilt(b0Hz, a0Hz, recg);
 plot(rtecg(x:(x+(5*newFs))), recg(x:(x+(5*newFs))));
 axis([x/newFs (x/newFs + 5) min(recg) max(recg)]);
 %% Find and mark R and S peaks in the hole ECG
 % Get R peaks
 Rtreshold = max(recg)*2/3;
 [~,locs_Rwave] = findpeaks(recg,'MinPeakHeight',Rtreshold,'MinPeakDistance',250);
 
 % Get Q peaks
 inverted_recg = -recg;
 Qtreshold = max(inverted_recg)*1/3;
 %plot(rtecg(x:(x+(5*newFs))), inverted_recg(x:(x+(5*newFs))));
 [~,locs_Qwave] = findpeaks(inverted_recg,'MinPeakHeight',Qtreshold,'MinPeakDistance',250);
 
%  RR_period = 0;
%  cnt = 0;
%  for i = 1:size(locs_Rwave, 2)-2
%      instant_period(i) = locs_Rwave(i+1)-locs_Rwave(i);
%     % if instant_period(i) < 1100
%       %   if instant_period(i) > 450
%             RR_period = RR_period + instant_period(i);
%             RR_vec(i) = RR_period;
%             cnt = cnt + 1;
%      %    end
%     % end
%  end
%  RR_period_mean = RR_period/cnt/newFs/60;
%  RR_freq = 1/RR_period_mean;
 RR_freq = size(locs_Rwave, 2)/size(recg,2)*60*1000;
% Get RR period
figure
hold on 
plot(rtecg, recg);
plot(locs_Rwave/newFs,recg(locs_Rwave),'rv','MarkerFaceColor','r');
plot(locs_Qwave/newFs,recg(locs_Qwave),'rs','MarkerFaceColor','b');
axis([x/newFs ((x/newFs) + 10) min(recg) max(recg)]);
legend('ECG Signal','R-waves','Q-waves')
xlabel('Seconds')
ylabel('Voltage(mV)')
str = ['RR frequency = ' num2str(RR_freq) ' b/min'];
dim = [.2 .5 .3 .4];
annotation('textbox',dim,'String',str, 'FitBoxToText', 'on');
title('R-wave and Q-wave 10 seconds')
hold off

figure
hold on 
plot(rtecg, recg);
plot(locs_Rwave/newFs,recg(locs_Rwave),'rv','MarkerFaceColor','r');
plot(locs_Qwave/newFs,recg(locs_Qwave),'rs','MarkerFaceColor','b');
grid on
legend('ECG Signal','R-waves','Q-waves')
xlabel('Seconds')
ylabel('Voltage(mV)')
title('R-wave and Q-wave Complete')
hold off

 %% Design passband filter to get QRS complex

[bband, aband] = butter(5, [10 25]*2/newFs, 'bandpass'); 
qrsComplex = filtfilt(bband, aband, recg);
figure('Name', 'QRS complex filter');
plot(rtecg, qrsComplex);
axis([x/newFs ((x/newFs) + 10) min(qrsComplex) max(qrsComplex)]);

%% Apply wavelet transform to get the QRS complex
[C, L] = wavedec(recg, 7, 'db5'); 
figure('Name', 'Wavelet transform to get QRS complex filter');
alpha = C;
betha = C;
delta = C;
teta = C;
C(1:L(3)) = 0;
C(L(5):L(9)) = 0;

WqrsComplex = waverec(C, L, 'db5');
plot(rtecg, WqrsComplex);
axis([x/newFs ((x/newFs) + 10) min(WqrsComplex) max(WqrsComplex)]);
%% Get delta waves

delta(L(2):L(9)) = 0;
%plot (delta);
deltaWaves = waverec(delta, L, 'db5');

%% Get alpha waves

alpha(1:L(2)) = 0;
alpha(L(4):L(9)) = 0;
%plot(alpha);
alphaWaves = waverec(alpha, L, 'db5');

%% Get betha waves

betha(1:L(3)) = 0;
betha(L(6):L(9)) = 0;
%plot (betha);
bethaWaves = waverec(betha, L, 'db5');

%% Get teta waves

teta(1:L(1)) = 0;
teta(L(3):L(9)) = 0;
%plot(teta);
tetaWaves = waverec(teta, L, 'db5');
%% Plot Waves
figure('Name', 'frequency reconstructed waves');
subplot(4, 1, 1)
plot(rtecg, deltaWaves);
axis([x/newFs ((x/newFs) + 10) min(deltaWaves) max(deltaWaves)]);
title('Delta Waves');

subplot(4, 1, 2)
plot(rtecg, alphaWaves);
axis([x/newFs ((x/newFs) + 10) min(alphaWaves) max(alphaWaves)]);
title('Alpha Waves');

subplot(4, 1, 3)
plot(rtecg, bethaWaves);
axis([x/newFs ((x/newFs) + 10) min(bethaWaves) max(bethaWaves)]);
title('Betha Waves');

subplot(4, 1, 4)
plot(rtecg, tetaWaves);
axis([x/newFs ((x/newFs) + 10) min(tetaWaves) max(tetaWaves)]);
title('Teta Waves');


%% Segmentation of ECG signal into complete periods using the peaks as reference

ecgseg = cell(size(locs_Rwave, 2), 1);
seg_offset = round(size(recg,2)/(2*size(locs_Rwave, 2)));

for i = 1:size(locs_Rwave, 2)
    ecgseg{i} = recg((locs_Rwave(i)-seg_offset):(locs_Rwave(i)+seg_offset));
end

figure
hold on;
for i = 1:size(ecgseg)
    plot(ecgseg{i});
end
hold off;