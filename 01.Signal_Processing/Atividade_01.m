% Universidade Federal de Santa Catarina
% Centro Tecnológico
% Departamento de Engenharia Elétrica e Eletrônica
% Disciplina: Introdução à Informática Médica (EEL7307)
% Professora: Christine Fredel Boos
% Alunos: Lucas Pereira Luiz - 13101258
%         Pedro Henrique Kappler Fornari - 13104320
% Semestre letivo: 2017/1
%
% Atividade_01 : Processamento de Sinais Digitais
%

clc;
clear;
close all;


%% ----- 1 -----

% Gets the directory and the filelist under the directory:
local_path = 'F:\Inform_Med\Atividade_01\';
filename = what([local_path,'\Base de dados - ECG e EEG\sinais_ecg']);

[file_id,~] = listdlg('PromptString','Select a file:','SelectionMode',...
                'single','ListString',filename.mat);

% Chooses between the 1st or 2nd line of the ECG matrix:
ecg_num = questdlg('The ECG files have two ECGs each. Choose one of them:',...
                     'ECG Selection','1','2','1');

% Imports the ECG data:
ECG = importdata([filename.path,'\',filename.mat{file_id}]);
ECG = ECG(str2num(ecg_num),:);

% Sets the time scale (based on the 128Hz sampling freq.:
time = (0:1/128:(size(ECG,2)-1)/128);


% Divides by the 200 gain to scale to the mV unit:
ECG = ECG/200;

% Plots the whole ECG raw data:
figure ('name','Whole ECG data','NumberTitle','off');
plot(time,ECG);
ylabel('Amplitude (mV)');
xlabel('Time (s)');
title('Whole ECG raw data');


%% ----- 2 -----

% Gets a pseudo-random starting point for the 5s plot interval:
x = round(rand*(1000000-1280));

% Plots a 5s interval of the ECG raw data:
figure('Name','5s ECG raw data','NumberTitle','off');
plot(time(x:(x+640)),ECG(x:(x+640)));
title('ECG raw data (divided by the gain) in a 5s interval');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
xlim([time(x) time(x+640)]);


%% ----- 3 -----

% Normalizes the ECG data between 0 and 1 values:
ECG_norm = (ECG -min(ECG))/max(ECG -min(ECG));

% Plots a 5s interval of the ECG normalized data:
figure('Name','5s ECG raw and normalized data','NumberTitle','off');
plot(time(x:(x+640)),ECG(x:(x+640)),'k');
hold on;
plot(time(x:(x+640)),ECG_norm(x:(x+640)),'r--');
legend('Raw','Normalized','Location','northeast');
title('ECG raw and normalized data in a 5s interval');
xlim([time(x) time(x+640)]);
ylabel('Amplitude (mV)');
xlabel('Time (s)');
hold off;


%% ----- 4 -----

% Resamples the ECG signal to 1kHz:
%ECG_1k = resample(ECG,1000,128);
time_1k = (0:1/1000:(size(ECG,2)-1)/128);
ECG_1k = spline(time,ECG,time_1k);


%% ----- 5 -----

% ...60Hz notch ...

% Notch filter design
[notch60_b, notch60_a] = iirnotch(60/(1000/2), 1/(1000/2), 1);

% Signal fitering
notched_60 = filtfilt(notch60_b,notch60_a,ECG_1k);

% Spectrum (128Hz)
Fs1 = 128;                       % Sampling frequency
T1 = 1/Fs1;                      % Sample time
L1 = size(ECG,2);                % Length of signal
t1 = (0:L1-1)*T1;                % Time vector

NFFT1 = 2^nextpow2(L1);          % Next power of 2 from length of L2
spec_ecg = fft(ECG,NFFT1)/L1;
f1 = Fs1/2*linspace(0,1,NFFT1/2+1);

% Spectrum (1kHz)
Fs2 = 1000;                     % Sampling frequency
T2 = 1/Fs2;                     % Sample time
L2 = size(ECG_1k,2);            % Length of signal
t2 = (0:L2-1)*T2;               % Time vector

NFFT2 = 2^nextpow2(L2);         % Next power of 2 from length of L2
spec_ecg_1k = fft(ECG_1k,NFFT2)/L2;
f2 = Fs2/2*linspace(0,1,NFFT2/2+1);

% 60Hz filter response
[notch60_h,notch60_w] = freqz(notch60_b,notch60_a,'whole',size(spec_ecg_1k,2));

% Spectrum of the filtered signal @60Hz
spec_notched_60 = fft(notched_60,NFFT2)/L2;

% Plots
figure('Name','ECG @ 128Hz and @ 1kHz','NumberTitle','off');
subplot(2,1,1);
plot(time,ECG,'k');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
hold on;
plot(time_1k,ECG_1k,'r-.');
title('ECG @ 128Hz and @ 1kHz');
plot(time_1k,notched_60,'Color',[0 0.5 0]','LineStyle','--');
legend('@ 128Hz','@ 1kHz','@ 1kHz filtered in 60Hz','Location','northeast');
xlim([time(x) time(x+640)]);
hold off;

subplot(2,1,2);

plot(f1,(2*abs(spec_ecg(1:NFFT1/2+1))),'k--');

title('Spectrum of the ECG sampled @ 128Hz and @ 1kHz');
xlabel('Frequency (Hz)');
ylabel('|Magnitude|') % left y-axis

hold on;

plot(f2,(2*abs(spec_ecg_1k(1:NFFT2/2+1))),'r-.');

plot(f2,(2*abs(spec_notched_60(1:NFFT2/2+1))),'Color',[0 0.5 0],'LineStyle',':');

xlim([0 500]);
plot(notch60_w*1000/(2*pi),abs(notch60_h)/10,'Color',[0.7 0 0.7]);

legend('@ 128Hz','@ 1kHz','@ 1kHz filtered in 60Hz','60Hz notch filter (\div10)','Location','northeast');
hold off;


% 60Hz Notch frequency response
figure('Name','60Hz Notch Filter Frequency Response','NumberTitle','off');
ax = plotyy(notch60_w*1000/(2*pi),db(abs(notch60_h)),notch60_w*1000/(2*pi),angle(notch60_h));
xlabel('Frequency (Hz)');
ylabel(ax(1), 'Magnitude (dB)');
ylabel(ax(2), 'Phase (rad)');
hold on;
title('60Hz Notch Filter Frequency Response');
xlim([0 500]);
legend('Magnitude (dB)','Phase (rad)');

% ...HighPass...

% HighPass filter design

[highpass_a,highpass_b]  = butter(4, 0.6/1000, 'high');

% Signal fitering
hpFiltered = filtfilt(highpass_a,highpass_b,ECG_1k);

% HP filter response
[hpFilt_h,hpFilt_w] = freqz(highpass_a,highpass_b,'whole',size(spec_ecg_1k,2));

% Spectrum of the HP filtered signal
spec_hpFiltered = fft(hpFiltered,NFFT2)/L2;

% Plots
figure('Name','ECG @ 128Hz and @ 1kHz','NumberTitle','off');

% Time subplot
subplot(2,1,1);

plot(time,ECG,'k');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
hold on;
title('ECG @ 128Hz and @ 1kHz');
plot(time_1k,ECG_1k,'r-.');

plot(time_1k,hpFiltered,'Color',[0 0.5 0],'LineStyle','--');

legend('@ 128Hz','@ 1kHz','@ 1kHz HP filtered','Location','northeast');
xlim([time(x) time(x+640)]);
hold off;

% Spectrum subplot
subplot(2,1,2);

plot(f1,(2*abs(spec_ecg(1:NFFT1/2+1))),'k--');
xlabel('Frequency (Hz)');
ylabel('|Magnitude|');
title('Spectrum of the ECG sampled @ 128Hz and @ 1kHz');
hold on;

plot(f2,(2*abs(spec_ecg_1k(1:NFFT2/2+1))),'r-.');

plot(f2,(2*abs(spec_hpFiltered(1:NFFT2/2+1))),'Color',[0 0.5 0],'LineStyle',':');

plot(hpFilt_w*1000/(2*pi),abs(hpFilt_h)/10,'Color',[0.7 0 0.7]);

xlim([0 500]);
legend('@ 128Hz','@ 1kHz','@ 1kHz HP filtered','HP filter (\div10)','Location','northeast');
hold off;

%HighPass frequency response
figure('Name','HighPass Filter Frequency Response','NumberTitle','off');
ax = plotyy(hpFilt_w*1000/(2*pi),db(abs(hpFilt_h)),hpFilt_w*1000/(2*pi),angle(hpFilt_h));
xlabel('Frequency (Hz)');
ylabel(ax(1), 'Magnitude (dB)');
ylabel(ax(2), 'Phase (rad)');
hold on;
title('HighPass Filter Frequency Response');
xlim([0 500]);
legend('Magnitude (dB)','Phase (rad)');


% ...Using both filters...

ECG_step5 = filtfilt(highpass_a,highpass_b,notched_60);

spec_fullyFilt= fft(ECG_step5,NFFT2)/L2;
 
% Plots
figure('Name','ECG noisy and filtered @ 1kHz','NumberTitle','off');
plot(time_1k,ECG_1k,'k');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
hold on;
title('ECG noisy and filtered @ 1kHz');
plot(time_1k,hpFiltered,'r--');
legend('Noisy Signal','Filtered Signal','Location','northeast');
xlim([time(x) time(x+640)]);
hold off;


%% ----- 6 -----

% Finding the peaks
[peak_value,peak_time] = findpeaks(ECG_step5,'MINPEAKHEIGHT',1,'MINPEAKDISTANCE',250);
ECG_bpm = size(peak_time,2)/size(time_1k,2)*60000;

% Plot
figure('Name','ECG R-peaks identification','NumberTitle','off');
subplot(2,1,1);
plot(time_1k,ECG_step5);
ylabel('Amplitude (mV)');
xlabel('Time (s)');
hold on;
plot(time_1k(peak_time),peak_value,'v','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off;
title('Full ECG R-peaks identificaition');
legend('Full ECG','R-peaks','Location','southeast','Orientation','horizontal');

subplot(2,1,2);
plot(time_1k,ECG_step5);
hold on;
ylabel('Amplitude (mV)');
xlabel('Time (s)');
plot(time_1k(peak_time),peak_value,'v','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','r');
hold off;
xlim([time(x) time(x+1280)]);
title('5s period ECG R-peaks identificaition');
legend('Full ECG','R-peaks','Location','southeast','Orientation','horizontal');
annotation('textbox', 'Position',[0.14,0.13,0.125,0.03],...
           'String',['Frequency = ',num2str(ECG_bpm),' bpm'],...
           'BackgroundColor','white');
%% ----- 7 -----

[highpass2_b,highpass2_a]  = butter(5, 10*2/1000, 'high');
[lowpass2_b,lowpass2_a]  = butter(8, 25*2/1000, 'low');

% HP2 & LP2 filtering == BP filtering
ECG_step7 = filtfilt(lowpass2_b,lowpass2_a,filtfilt(highpass2_b,highpass2_a,ECG_step5));

% HP2 & LP2 filter response
[hpFilt2_h,hpFilt2_w] = freqz(highpass2_b,highpass2_a,size(spec_ecg_1k,2));
[lpFilt2_h,lpFilt2_w] = freqz(lowpass2_b,lowpass2_a,size(spec_ecg_1k,2));

% Spectrum of the "BP" filtered signal
spec_ECG_step7 = fft(ECG_step7,NFFT2)/L2;


% Plots
figure('Name','Bandpass Filtering (QRS)','NumberTitle','off');

% Time subplot
subplot(2,1,1);
plot(time_1k,ECG_1k,'k');
title('ECG @ 1k');
hold on;

plot(time_1k,ECG_step7,'r');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
legend('@ 128Hz','@ 1kHz BP filtered','Location','northeast');
xlim([time(x) time(x+640)]);
hold off;

% Spectrum subplot
subplot(2,1,2);
plot(f2,(2*abs(spec_ecg_1k(1:NFFT2/2+1))),'k');

title('Spectrum of the ECG @ 1kHz and BP filtered @ 1kHz');
xlabel('Frequency (Hz)');
ylabel('|Magnitude|')
hold on;

plot(f2,(2*abs(spec_ECG_step7(1:NFFT2/2+1))),'r--');

plot(hpFilt2_w*1000/(2*pi),abs(hpFilt2_h)/10,'Color',[0.7 0.0 0.7]);
plot(lpFilt2_w*1000/(2*pi),abs(lpFilt2_h)/10,'Color',[0.0 0.7 0.7]);

xlim([0 500]);
legend('@ 128Hz','@ 1kHz BP filtered','HP filter (\div10)','LP filter (\div10)','Location','northeast');
hold off;

% HighPass & LowPass frequency response
figure('Name','HighPass Filter Frequency Response','NumberTitle','off');
ax = plotyy(hpFilt2_w*1000/(2*pi),db(abs(hpFilt2_h)),hpFilt2_w*1000/(2*pi),angle(hpFilt2_h));
axis([0 500 -200 1]);
title('HighPass Filter Frequency Response');
legend('Magnitude (dB)','Phase (rad)');
xlabel('Frequency (Hz)');
ylabel(ax(1), 'Magnitude (dB)');
ylabel(ax(2), 'Phase (rad)');

figure('Name','LowPass Filter Frequency Response','NumberTitle','off');
ax = plotyy(lpFilt2_w*1000/(2*pi),db(abs(lpFilt2_h)),lpFilt2_w*1000/(2*pi),angle(lpFilt2_h));
title('LowPass Filter Frequency Response');
axis([0 500 -200 1]);
legend('Magnitude (dB)','Phase (rad)');
xlabel('Frequency (Hz)');
ylabel(ax(1), 'Magnitude (dB)');
ylabel(ax(2), 'Phase (rad)');



%% ----- 8 -----

[C_ecg,L_ecg]=wavedec(ECG_step5,6,'sym6');

C_ecg_new = [C_ecg(1:L_ecg(1))*0 C_ecg((L_ecg(1)+1):(L_ecg(1)+L_ecg(2)+...
    L_ecg(3))) 0*C_ecg((L_ecg(1)+L_ecg(2)+L_ecg(3)+1):end)]; %7.125Hz to 31.25 Hz

ECG_step8 = waverec(C_ecg_new,L_ecg,'sym6');

spec_ECG_step8 = fft(ECG_step8,NFFT2)/L2;



% Plots
figure('Name','Bandpass Filtering','NumberTitle','off');

% Time subplot
subplot(2,1,1);
plot(time_1k,ECG_1k,'k');

title('ECG @ 1k');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
hold on;

plot(time_1k,ECG_step8,'b');

plot(time_1k,ECG_step7,'r-.');

legend('Interpolated','Wavelet','BP filtered','Location','northeast');
xlim([time(x) time(x+640)]);
hold off;



% Spectrum subplot
subplot(2,1,2);
plot(f2,(2*abs(spec_ecg_1k(1:NFFT2/2+1))),'k');

title('Spectrum of the ECG @ 1kHz');
xlabel('Frequency (Hz)');
ylabel('|Magnitude|');
hold on;

plot(f2,(2*abs(spec_ECG_step8(1:NFFT2/2+1))),'b');

plot(f2,(2*abs(spec_ECG_step7(1:NFFT2/2+1))),'r-.');

xlim([0 500]);
legend('Interpolated','Wavelet','BP filtered','Location','northeast');
hold off;

%% ----- 9 -----
% EEG wavelets

% Gets the directory and the filelist under the directory:
%local_path = 'F:\Inform_Med\Atividade_01\';
filename = what([local_path,'\Base de dados - ECG e EEG\sinais_eeg']);

[file_id,~] = listdlg('PromptString','Select a file:','SelectionMode',...
                'single','ListString',filename.mat);


legend_tick = (1:64)*250;
legend_str = {'Fc5.','Fc3.','Fc1.','Fcz.','Fc2.','Fc4.','Fc6.','C5..',...
              'C3..','C1..','Cz..','C2..','C4..','C6..','Cp5.','Cp3.',...
              'Cp1.','Cpz.','Cp2.','Cp4.','Cp6.','Fp1.','Fpz.','Fp2.',...
              'Af7.','Af3.','Afz.','Af4.','Af8.','F7..','F5..','F3..',...
              'F1..','Fz..','F2..','F4..','F6..','F8..','Ft7.','Ft8.',...
              'T7..','T8..','T9..','T10.','Tp7.','Tp8.','P7..','P5..',...
              'P3..','P1..','Pz..','P2..','P4..','P6..','P8..','Po7.',...
              'Po3.','Poz.','Po4.','Po8.','O1..','Oz..','O2..','Iz..'};
          
% Imports the EEG data:
EEG = importdata([filename.path,'\',filename.mat{file_id}]);
EEG = EEG((1:64),:);

% Sets the time scale (based on the 160Hz sampling freq.):
time_eeg = (0:1/160:(size(EEG,2)-1)/160);

Fs_eeg = 160;                     % Sampling frequency
T_eeg = 1/Fs_eeg;                     % Sample time
L2_eeg = size(EEG,2);            % Length of signal
t_eeg = (0:L2_eeg-1)*T_eeg;               % Time vector

NFFT_eeg = 2^nextpow2(L2_eeg);         % Next power of 2 from length of L2_eeg
f_eeg = Fs_eeg/2*linspace(0,1,NFFT_eeg/2+1);


for i=1:size(EEG,1)
    [C_eeg(i,:),L_eeg(i,:)]=wavedec(EEG(i,:),4,'sym6');

    C_eeg_d(i,:) = [C_eeg(i,1:L_eeg(i,1)) ...
        0*C_eeg(i,(L_eeg(i,1)+1):end)]; %(0Hz to 5Hz)
    
    C_eeg_t(i,:) = [C_eeg(i,1:L_eeg(i,1))*0 ...
        C_eeg(i,(L_eeg(i,1)+1):(L_eeg(i,1)+L_eeg(i,2))) ...
        0*C_eeg(i,(L_eeg(i,1)+L_eeg(i,2)+1):end)]; %(5Hz to 10Hz)
    
    C_eeg_a(i,:) = [C_eeg(i,1:L_eeg(i,1))*0 ...
        C_eeg(i,(L_eeg(i,1)+1):(L_eeg(i,1)+L_eeg(i,2)+L_eeg(i,3))) ...
        0*C_eeg(i,(L_eeg(i,1)+L_eeg(i,2)+L_eeg(i,3)+1):end)]; %(5Hz to 20Hz)
    
    C_eeg_b(i,:) = [C_eeg(i,1:(L_eeg(i,1)+L_eeg(i,2)))*0 ...
        C_eeg(i,(L_eeg(i,1)+L_eeg(i,2)+1):(L_eeg(i,1)+L_eeg(i,2)+L_eeg(i,3)+L_eeg(i,4))) ...
        0*C_eeg(i,(L_eeg(i,1)+L_eeg(i,2)+L_eeg(i,3)+L_eeg(i,4)+1):end)]; %(10Hz to 40Hz)

    EEG_d(i,:) = waverec(C_eeg_d(i,:),L_eeg(i,:),'sym6');
    spec_eeg_d(i,:) = fft(EEG_d(i,:),NFFT_eeg)/L2_eeg;
    
    EEG_t(i,:) = waverec(C_eeg_t(i,:),L_eeg(i,:),'sym6');
    spec_eeg_t(i,:) = fft(EEG_t(i,:),NFFT_eeg)/L2_eeg;
    
    EEG_a(i,:) = waverec(C_eeg_a(i,:),L_eeg(i,:),'sym6');
    spec_eeg_a(i,:) = fft(EEG_a(i,:),NFFT_eeg)/L2_eeg;
    
    EEG_b(i,:) = waverec(C_eeg_b(i,:),L_eeg(i,:),'sym6');
    spec_eeg_b(i,:) = fft(EEG_b(i,:),NFFT_eeg)/L2_eeg;
end

spec_eeg_d =  fft(sum(EEG_d),NFFT_eeg)/L2_eeg;
spec_eeg_t =  fft(sum(EEG_t),NFFT_eeg)/L2_eeg;
spec_eeg_a =  fft(sum(EEG_a),NFFT_eeg)/L2_eeg;
spec_eeg_b =  fft(sum(EEG_b),NFFT_eeg)/L2_eeg;
spec_eeg = fft(sum(EEG),NFFT_eeg)/L2_eeg;


figure('Name','Raw EEG data','NumberTitle','off');
hold on;
for i=1:64
    plot(time_eeg(4*160:14*160),250*(65-i)+EEG(i,4*160:14*160));
    ylim([0 65*250]);
    xlim([4 14]);
end
set(gca, 'YTick', legend_tick, 'YTickLabel', fliplr(legend_str),'FontSize',8.0,'FontName','Courier New');
title('EEG raw waves'' signals');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
hold off;


figure('Name','EEG signals spectrum','NumberTitle','off');
subplot(2,2,1);
hold on;
plot(f_eeg,(2*abs(spec_eeg(1:NFFT_eeg/2+1))),'k');
title('EEG \delta waves'' spectrum');
plot(f_eeg,(2*abs(spec_eeg_d(1:NFFT_eeg/2+1))),'r');
legend('Raw','\delta Waves','Location','northeast');
xlabel('Frequency (Hz)');
ylabel('|Magnitude|');
hold off;

subplot(2,2,2);
hold on;
plot(f_eeg,(2*abs(spec_eeg(1:NFFT_eeg/2+1))),'k');
title('EEG \theta waves'' spectrum');
plot(f_eeg,(2*abs(spec_eeg_t(1:NFFT_eeg/2+1))),'r');
legend('Raw','\theta Waves','Location','northeast');
xlabel('Frequency (Hz)');
ylabel('|Magnitude|');
hold off;

subplot(2,2,3);
hold on;
plot(f_eeg,(2*abs(spec_eeg(1:NFFT_eeg/2+1))),'k');
title('EEG \alpha waves'' spectrum');
plot(f_eeg,(2*abs(spec_eeg_a(1:NFFT_eeg/2+1))),'r');
legend('Raw','\alpha Waves','Location','northeast');
xlabel('Frequency (Hz)');
ylabel('|Magnitude|');
hold off;

subplot(2,2,4);
hold on;
plot(f_eeg,(2*abs(spec_eeg(1:NFFT_eeg/2+1))),'k');
title('EEG \beta waves'' spectrum');
plot(f_eeg,(2*abs(spec_eeg_b(1:NFFT_eeg/2+1))),'r');
legend('Raw','\beta Waves','Location','northeast');
xlabel('Frequency (Hz)');
ylabel('|Magnitude|');
hold off;


figure('Name','EEG delta signals','NumberTitle','off');
hold on;
for i=1:64
    plot(time_eeg(4*160:14*160),250*(65-i)+EEG_d(i,4*160:14*160));
    ylim([0 65*250]);
    xlim([4 14]);
end
set(gca, 'YTick', legend_tick, 'YTickLabel', fliplr(legend_str),'FontSize',8.0,'FontName','Courier New');
title('EEG \delta waves'' signals');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
hold off;


figure('Name','EEG theta signals','NumberTitle','off');
hold on;
for i=1:64
    plot(time_eeg(4*160:14*160),250*(65-i)+EEG_t(i,4*160:14*160));
    ylim([0 65*250]);
    xlim([4 14]);
end
set(gca, 'YTick', legend_tick, 'YTickLabel', fliplr(legend_str),'FontSize',8.0,'FontName','Courier New');
title('EEG \theta waves'' signals');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
hold off;

figure('Name','EEG alpha signals','NumberTitle','off');
hold on;
for i=1:64
    plot(time_eeg(4*160:14*160),250*(65-i)+EEG_a(i,4*160:14*160));
    ylim([0 65*250]);
    xlim([4 14]);
end
set(gca, 'YTick', legend_tick, 'YTickLabel', fliplr(legend_str),'FontSize',8.0,'FontName','Courier New');
title('EEG \alpha waves'' signals');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
hold off;


figure('Name','EEG beta signals','NumberTitle','off');
hold on;
for i=1:64
    plot(time_eeg(4*160:14*160),250*(65-i)+EEG_b(i,4*160:14*160));
    ylim([0 65*250]);
    xlim([4 14]);
end
set(gca, 'YTick', legend_tick, 'YTickLabel', fliplr(legend_str),'FontSize',8.0,'FontName','Courier New');
title('EEG \beta waves'' signals');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
hold off;



%% ----- 10 -----

seg_size = round(size(ECG_step5,2)/(2*size(peak_value,2)))-1;

ECG_seg = zeros((size(peak_time,2)-2),seg_size*2+1);
for i=1:(size(peak_time,2)-2)
    ECG_seg(i,:) = ECG_step5((peak_time(i+1)-seg_size):(peak_time(i+1)+seg_size))';
end

figure('Name','ECG signal segmentation','NumberTitle','off');
hold on;
for j=1:100
    plot(ECG_seg(j,:));
end
title('100 first segments of the ECG @ 1k');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
hold off;

sinais_ECG = ECG_seg(1:100,:);

save([local_path,'\sinais_ECG.mat'],'sinais_ECG');