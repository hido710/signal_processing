%% Initialize
clc;
clear;
close all;

%% time delay 구하기

%%% Setting ---------------------------------------------------------------
[x_wav,Fs] = audioread('./x_2x8.wav');
[nsample, nch] = size(x_wav);

SS = 343.3;
rad = pi/180;
load('location_sensor');
load('location_source'); % Ideal Source 위치

mic_anchor = [1.5 2 1];    % mic array 중심
L = 19;                    % -90~90 10도 간격
d_azimuth = 10;            % 10도 간격

%%% virtual source location ----------------------------------------------
s_tmp = zeros(L, 3);

for dl = 1 : L             % [-90 -80 ... 80 90]
    s_tmp(dl, :) = [sin((dl-1)*d_azimuth*rad) -cos((dl-1)*d_azimuth*rad) 0] + mic_anchor;
end

figure(); hold on;
for ch = 1 : nch
    plot3(locationSensor{ch}(1), locationSensor{ch}(2), locationSensor{ch}(3), 'bo')
end
for i = 1 : length(locationSource)
    plot3(locationSource{i}(1), locationSource{i}(2), locationSource{i}(3), 'r*')
end
for dl = 1 : L
    plot3(s_tmp(dl,1), s_tmp(dl,2), s_tmp(dl,3), 'go')
end
xlabel('x axis (m)'); ylabel('y axis (m)'); zlabel('z axis (m)');
axis([1.5 3.5 1 3]);
title('Location of Sensor and Virtual Sources');
grid; hold off;

%%% time delay -----------------------------------------------------------
MicDist = zeros(nch, L);  % mic와 가상 source 간의 distance
MicDelay = zeros(nch, L); % mic와 가상 source 간의 time delay

for dl = 1 : L
    for ch = 1 : nch
        MicDist(ch, dl) = norm(locationSensor{ch} - s_tmp(dl, :));
        MicDelay(ch, dl) = MicDist(ch,dl)/SS;
    end
end

%% E(l) 계산하기

%%% STFT Setting ---------------------------------------------------------
winL = 512;  % number of sample to windowing
nfft = winL;
nshift = 128;
shiftdiv = winL/nshift;

if shiftdiv == 2
    win = sin(pi*([0:1:winL-1]'+0.5)/winL); %1/2 shift
end
if shiftdiv == 4
    win = sqrt(2/3)*hanning(winL,'periodic');%1/4shift
end

nover = winL - nshift; % number of sample to window overlap

nhfft = nfft/2+1; % half + 1 of fftsize
nframe = floor((nsample - nover)/nshift); %number of frame on STFT matrix

x_frame = zeros(winL,nch);
X = zeros(nch,nhfft,nframe); % frame 단위 fft 저장
X_ideal = zeros(nch,nhfft,nframe); % frame 단위 fft 저장

Power = zeros(L, nframe);
azimuth = zeros(1, nframe);

%%% Ideal source ----------------------------------------------------------
[x_ideal,Fs] = audioread('./s_2x8_source_1.wav'); % white 없는 ideal 신호

% STFT
for dfr = 1 : nframe   % frame by frame
    x_frame = x_ideal((dfr-1)*nshift+1 : (dfr-1)*nshift+winL, :);
    for ch = 1 : nch
        Xtmp = fft([x_frame(:,ch).* win; zeros(nfft-winL,1)]);  % 양방향 fft
        X_ideal(ch,:,dfr) = Xtmp(1:nhfft);                      % 단방향 fft
    end
end

figure(2); sgtitle('Steered Response Power');

subplot(1, 2, 2);
ch_spec = 4; % ch4에 대해서만 살펴봄
imagesc(10*log10(abs(squeeze(X_ideal(ch_spec,:,:))).^2));
colormap jet; colorbar;
title('STFT of 4th mic');
xlabel('Frame Index'); ylabel('Freq bin');
axis xy; axis([1 nframe 1 nhfft]); caxis([-80 0]);

%%% Source SRP -----------------------------------------------------------
for dfr = 1 : nframe   % frame by frame
    % STFT
    x_frame = x_wav((dfr-1)*nshift+1 : (dfr-1)*nshift+winL, :);
    for ch = 1 : nch
        Xtmp = fft([x_frame(:,ch).* win; zeros(nfft-winL,1)]); % 양방향 fft
        X(ch,:,dfr) = Xtmp(1:nhfft);                           % 단방향 fft
    end
    
    %%% E(l) -------------------------------------------------------------
    for dl = 1 : L
        power_tmp2 = 0;
        for df = 1 : nhfft
            power_tmp1 = 0;
            for ch = 1 : nch
                exp_tmp = exp(1j*(2*pi*(df-1)/nfft*Fs*MicDelay(ch, dl)));
                power_tmp1 = power_tmp1 + X(ch, df, dfr).*exp_tmp;
            end
            power_tmp2 = power_tmp2 + power_tmp1.^2;
        end
        Power(dl, dfr) = power_tmp2;
    end
    
    %%% argmax{E(l)} -----------------------------------------------------
    [Max, dl_argmax] = max(Power(:,dfr), [], 'linear');
    azimuth(dfr) = -90 + (dl_argmax-1)*d_azimuth;
    fprintf('[frame %5d] azimuth : %d\n', dfr, azimuth(dfr));
        
    %%% Power map and Spectrogram --------------------------------------------
    % Power    
    subplot(1, 2, 1);
    plot([-90:10:90], (abs(Power(:,dfr)/Max)).');
    hold on; grid on;
    plot(azimuth(dfr), 1, '*r');
    xlabel('azimuth (degree)'); ylabel('Steered Response Power');
    axis([-90 90 0 1]);
    title('Power');    
    hold off;
    
    % Spectrogram
    figure(2);
    subplot(1, 2, 2);
    xlim([dfr dfr+200]); % 출력할 단위 spectrogram
    title(['azimuth : ' num2str(azimuth(dfr))]);
    
end



