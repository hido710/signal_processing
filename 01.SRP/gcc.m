%% Clear all
clc;
clear;
close all;

%% Plot location sources and sensors
load('location_sensor');
load('location_source');

len_sen = length(locationSensor);
len_sour = length(locationSource);

figure(); hold on;
for i = 1:len_sen
    plot3(locationSensor{i}(1), locationSensor{i}(2), locationSensor{i}(3), 'bo')
end
for j = 1:len_sour
    plot3(locationSource{j}(1), locationSource{j}(2), locationSource{j}(3), 'r*')
end
xlabel('x axis (m)'); ylabel('y axis (m)'); zlabel('z axis (m)');
title('Location of Sensors and Sources');
grid; hold off;

%% Calculate Ideal Sample delay

SS = 343.3;
Fs = 16000; % [x_wav,Fs] = audioread('--.wav');에서 얻어냄, s1, s2의 Fs 같음

MicDist1=zeros(len_sen,3);
MicDist2=zeros(len_sen,3);
dist1=zeros(len_sen,1);
dist2=zeros(len_sen,1);

for i = 1 : len_sen
    MicDist1(i,:) = locationSensor{i} - locationSource{1};
    dist1(i) = sqrt(sum(MicDist1(i,:).*MicDist1(i,:)));
    MicDist2(i,:) = locationSensor{i} - locationSource{2};
    dist2(i) = sqrt(sum(MicDist2(i,:).*MicDist2(i,:)));
end

Ideal1 = zeros(len_sen,1);
Ideal2 = zeros(len_sen,1);
for i = 1 : len_sen
    Ideal1(i) = dist1(i)/SS *Fs;
    Ideal2(i) = dist2(i)/SS *Fs;
end

fprintf('Ideal Delay _s1 = %f\n', Ideal1(1)- Ideal1(8));
fprintf('Ideal Delay _s2 = %f\n', Ideal2(1)- Ideal2(8));

%% GCC를 이용하여 sample delay 추정하기
for opt = 1:2
    if opt == 1
        [x_wav,Fs] = audioread('./s_2x8_source_1.wav');
        [nsample, nch] = size(x_wav);
    else
        [x_wav,Fs] = audioread('./s_2x8_source_2.wav');
        [nsample, nch] = size(x_wav);
    end
    
    %%% STFT Setting ---------------------------------------------------
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
    X = zeros(nch,nhfft);
    
    Dtmp = zeros(1, nframe);
    R = zeros(winL, 1, nframe);
    
    X_buffer1 = zeros(nch,nhfft,nframe); % frame 단위 fft 저장
    X_buffer2 = zeros(nch,nhfft,nframe);
    
    if opt == 1
        G1 = zeros(1, nhfft, nframe);
    else
        G2 = zeros(1, nhfft, nframe);
    end
    G = zeros(1, nhfft);
    
    %%% STFT ------------------------------------------------------------
    for dfr = 1 : nframe
        x_frame = x_wav((dfr-1)*nshift+1 : (dfr-1)*nshift+winL, :);
        for ch = 1 : nch
            Xtmp = fft(x_frame(:,ch).* win);  % 양방향 fft
            X(ch,:) = Xtmp(1:nhfft).';        % 단방향 fft
        end
        
        %%% GCC
        if opt == 1
            X_buffer1(:,:,dfr) = X;
            G1(:,:, dfr) = X(1,:) .* conj(X(8,:));
            G = G1(:,:, dfr);
        else
            X_buffer2(:,:,dfr) = X;
            G2(:,:, dfr) = X(1,:) .* conj(X(8,:));
            G = G2(:,:, dfr);
        end
        
        %%% ISTFT
        R(:,:,dfr) = ifft([(G.'); (G(1, end-1:-1:2)')], 'symmetric');
        
        [Max, locs] = max(R(:,:,dfr), [], 'linear');
        if Max ~= 0
            % 실제 delay는 0부터 시작하는데 locs는 1부터 시작하므로 1 빼줌
            Dtmp(dfr) = locs-1;
        end
    end
    
    if opt == 1
        % frame 길이에 따라서 다를 수 있음. non-speech 구간때문에 약간의 오차 존재
        Delay1 = mode(Dtmp); % 최빈값
    else
        % -delay_max ~ delay_max까지만 유효함 -> ring buffer로 돌기 때문에 음수 값 얻으려면 -winL 해야함
        Delay2 = mode(Dtmp) - winL; % 최빈값
    end
end

fprintf("Sample Delay_s1 = %f\n", Delay1);
fprintf("Sample Delay_s2 = %f\n", Delay2);

%% time domain Interpolation를 GCC에 적용하여 sample delay 추정하기
R1_inp = zeros(16*nfft, 1, nframe);
R2_inp = zeros(16*nfft, 1, nframe);
Dtmp1 = zeros(1,nframe);
Dtmp2 = zeros(1,nframe);

for dfr = 1:nframe
    if mod(nfft/2, 2) == 0 % N = even
        R1_inp(:,:,dfr) = ifft([(G1(1, 1:end-1, dfr).'); zeros(nfft/2*15, 1); G1(1, end, dfr); zeros(nfft/2*15, 1); (G1(1, end-1:-1:2, dfr)')], 'symmetric');
        R2_inp(:,:,dfr) = ifft([(G2(1, 1:end-1, dfr).'); zeros(nfft/2*15, 1); G2(1, end, dfr); zeros(nfft/2*15, 1); (G2(1, end-1:-1:2, dfr)')], 'symmetric');        
    else                   % N = odd
        R1_inp(:,:,dfr) = ifft([(G1(1, : , dfr).'); zeros(nfft*15, 1); (G1(1, end-1:-1:2, dfr)')], 'symmetric');
        R2_inp(:,:,dfr) = ifft([(G2(1, : , dfr).'); zeros(nfft*15, 1); (G2(1, end-1:-1:2, dfr)')], 'symmetric');
    end
    
    [Max1, locs1] = max(R1_inp(:,:,dfr), [], 'linear');
    if Max1 ~= 0
        Dtmp1(dfr) = (locs1-1)/16; % interpolation 때문에 16배 길어진 길이 조정
    end
    [Max2, locs2] = max(R2_inp(:,:,dfr), [], 'linear');
    if Max2 ~= 0
        Dtmp2(dfr) = (locs2-1)/16 - winL;
    end
end

Delay1_inp = mode(Dtmp1); % 최빈값
Delay2_inp = mode(Dtmp2); % 최빈값

fprintf("Interpolation Delay_s1 = %f\n", Delay1_inp);
fprintf("Interpolation Delay_s2 = %f\n", Delay2_inp);

%% SRP를 이용하여 각도 추정하기

% mic array 중심
mic_anchor = [1.5 2 1];
source_arc = zeros(181, 3);
rad = pi/180;
L = 181; % Source를 L points로 dicretize

for dl = 1 : L
    source_arc(dl, :) = [sin(dl*rad) -cos(dl*rad) 0] + mic_anchor;
end

% tau(m,l) = mth mic와 lth 가상 source간의 거리, sample delay
ArcDist = zeros(nch, L);
ArcDelay = zeros(nch, L);

for dl = 1 : L
    for ch = 1 : nch
        ArcDist(ch, dl) = norm(locationSensor(ch) - source_arce(dl, :));
        ArcDelay(ch, dl) = ArcDist(ch,dl)/SS * Fs;
    end
end

% E(l)
Power = zeros(L, nframe);
X_buffer = zeros(nch,nhfft,nframe);

% 수식 (6.21) discrete ver.
for opt = 1 : 2
    if opt == 1
        X_buffer = X_buffer1;
    else
        clear X_buffer;
        X_buffer = X_buffer2;
    end
    
    for dfr = 1 : nframe
        for dl = 1 : L
            power_tmp2 = 0;
            for df = 1 : nhfft
                power_tmp1 = 0;
                for ch = 1 : nch
                    power_tmp1 = power_tmp1 + X_buffer(ch, df, dfr).*exp(j*(2*pi*df/nfft*ArcDelay(ch, dl)));
                end
                power_tmp2 = power_tmp2 + power_tmp1.^2;
            end
            Power(dl, dfr) = power_tmp2;
        end
    end
    
    if opt == 1
        Power1 = Power;
    else
        Power2 = Power;
    end
end








