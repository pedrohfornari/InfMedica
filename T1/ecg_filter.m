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
 x = round(100000*rand());
 %plot
 plot(tecg(x:(x+(5*fs))), ecg(x:(x+(5*fs))));
 
 %% Nomalize signal
 %maxecg = max(ecg);
 %minecg = min(ecg);
 %ecg_norm = (ecg-minecg)/(maxecg-minecg);
 ecg_norm = mat2gray(ecg);
 
 figure
 plot(tecg(x:(x+(5*fs))), ecg(x:(x+(5*fs))));
 hold on
 plot(tecg(x:(x+(5*fs))), ecg_norm(x:(x+(5*fs))));
 hold off
 
 %% Resample signal
 %ups = 1024;
 %dns = fs;
 %ecgr = resample(ecg, ups, dns);
 newFs = 1000;
 
 % create resampled time vector
 rtecg = 0:(1/newFs):(size(ecg,2)-1)/fs;
 %resample by spline
 recg = spline(tecg, ecg_norm, rtecg);
 
 % plot resampled signal
 figure
 plot(rtecg(x:(x+(5*newFs))), recg(x:(x+(5*newFs))));
 
 %% Design Notch filter for 0 and 60 Hz
 notchspec_0hz = fdesign.notch('N,F0,Q',2,0.0000001,10,newFs);
 notchfilt = design(notchspec_0hz,'SystemObject',true);
fvt= fvtool(notchfilt,'Color','white');
 
