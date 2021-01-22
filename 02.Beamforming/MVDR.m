%% Initialize
clear; close all;  clc;

%% Data Setting
[x_wav,Fs] = audioread('./x_2x8.wav');
[nsample, nch] = size(x_wav);

[s_wav,~] = audioread('./s_2x8_source_1.wav');
[u_wav,~] = audioread('./s_2x8_source_2.wav');

%%% STFT Setting ----------------------------------------------------------
winL = 512;  % number of sample to windowing
nfft = winL;
nshift = 128;
shiftdiv = winL/nshift;

nover = winL - nshift; % number of sample to window overlap
nhfft = nfft/2+1; % half + 1 of fftsize
nframe = floor((nsample - nover)/nshift); %number of frame on STFT matrix

cmin = -150; cmax = -40;

if shiftdiv == 2
    win = sin(pi*([0:1:winL-1]'+0.5)/winL); %1/2 shift
end
if shiftdiv == 4
    win = sqrt(2/3)*hanning(winL,'periodic');%1/4shift
end

x_frame = zeros(winL,nch);
Xtmp = zeros(nch, nhfft);
X = zeros(nch, nhfft, nframe);      % X(t,f)

%%% STFT of x(t) ----------------------------------------------------------
for dfr = 1 : nframe
    x_frame = x_wav((dfr-1)*nshift+1 : (dfr-1)*nshift+winL, :);
    for ch = 1 : nch
        Xtmp = fft([x_frame(:,ch).* win; zeros(nfft-winL,1)]);
        X(ch,:,dfr) = Xtmp(1:nhfft);
    end
end

%%% Steering Vector -------------------------------------------------------
load('location_sensor.mat');
load('location_source.mat');

SS = 343.3;              % Speed of Sound

MicDist = zeros(1, nch); % target source인 1st source와 mic array의 거리
for ch = 1 : nch
   MicDist(ch) = norm(locationSensor{ch} - locationSource{1}, 2);    
end

H = zeros(nch, nhfft);   % steering vector
for df = 1 : nhfft
    for ch = 1 : nch        
        H(ch, df) = exp(-1i*2*pi*(df-1)/nfft*Fs*MicDist(ch)/SS);
    end
    H(:,df) = H(:,df)./ (MicDist).';
end

%%% Spectrogram of Mic signal and Target source ---------------------------
figure(1); sgtitle('Target');

subplot(1, 2, 1);
spectrogram(x_wav(:,1), winL, nover, winL, Fs, 'yaxis');
colormap jet; colorbar;
caxis([cmin cmax]); % 얘를 맞춰줘야 확인하기 쉽다~
title('x(t)');

subplot(1, 2, 2);
spectrogram(s_wav(:,8), winL, nover, winL, Fs, 'yaxis');
colormap jet; colorbar;
caxis([cmin cmax]);
title('s(t)');

%% Delay and Sum Beamformer (DS)
Y_ds = zeros(1, nhfft, nframe);     % Y(t,f) for DS beamformer
W_ds = zeros(nch, nhfft);           % W(f)   for DS beamformer

%%% Filter Matrix ---------------------------------------------------------
W_ds = 1/nch*H;

%%% Beamformer Output -----------------------------------------------------
y_ds = zeros(nsample, 1);
y_tmp = zeros(winL, 1);

for dfr = 1 : nframe
    for df = 1 : nhfft
        Y_ds(1, df, dfr) = (W_ds(:, df))' * X(:, df, dfr);
    end
    % ISTFT
    y_tmp = ifft([(Y_ds(1, :, dfr)).'; (Y_ds(1, end-1:-1:2, dfr))'], 'symmetric');
    y_ds((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) = y_ds((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) + y_tmp(1:winL).*win;
end

%%% Spectrogram -----------------------------------------------------------
figure(2); sgtitle('DS Beamformer');

subplot(1, 1, 1);
spectrogram(y_ds, winL, nover, winL, Fs, 'yaxis');
colormap jet; colorbar;
caxis([cmin cmax]);
title('y_{ds}(t)');

%% Minimum Variance Distortionless Response Beamformer (MVDR)
Y_mvdr_u = zeros(nch, nhfft, nframe);   % Y(t,f) for MVDR beamformer with.R_u
Y_mvdr_x = zeros(nch, nhfft, nframe);   % Y(t,f) for MVDR beamformer with.R_x
W_mvdr_u = zeros(nch, nhfft);         % W(f) for MVDR beamformer with.R_u
W_mvdr_x = zeros(nch, nhfft);         % W(f) for MVDR beamformer with.R_x
delta = 1e-3;

%%% Autocorrelation -------------------------------------------------------
% STFT of noise u(t)
u_frame = zeros(winL,nch);
Utmp = zeros(nch, nhfft);
U = zeros(nch, nhfft, nframe);      % U(t,f)
for dfr = 1 : nframe
    u_frame = u_wav((dfr-1)*nshift+1 : (dfr-1)*nshift+winL, :);
    for ch = 1 : nch
        Utmp = fft([u_frame(:,ch).* win; zeros(nfft-winL,1)]);
        U(ch,:,dfr) = Utmp(1:nhfft);
    end
end

% Autocorrelation of noise u(t)
R_tmp = zeros(nch, nch);
R_u = zeros(nch, nch, nhfft);

for df = 1 : nhfft
    for dfr = 1 : nframe
        R_tmp = U(:, df, dfr) * (U(:, df, dfr))';
        R_u(:, :, df) = R_u(:, :, df) + R_tmp;
    end
    R_u(:, :, df) = 1/nframe*(R_u(:, :, df)) + delta*eye(nch);
end

% Autocorrelation of mic signal x(t)
R_tmp = 0 * R_tmp;
R_x = zeros(nch, nch, nhfft);

for df = 1 : nhfft
    for dfr = 1 : nframe
        R_tmp = X(:, df, dfr) * (X(:, df, dfr))';
        R_x(:, :, df) = R_x(:, :, df) + R_tmp;
    end
    R_x(:, :, df) = 1/nframe*(R_x(:, :, df)) + delta*eye(nch);
end

%%% Filter Maxtrix --------------------------------------------------------
for df = 1 : nhfft
   W_mvdr_u(:, df) = (R_u(:,:,df)\H(:,df)) ./ (H(:,df)'*(R_u(:,:,df)\H(:,df)));
   W_mvdr_x(:, df) = (R_x(:,:,df)\H(:,df)) ./ (H(:,df)'*(R_x(:,:,df)\H(:,df)));
end

%%% Beamformer Output -----------------------------------------------------
% with R_u
y_tmp2 = zeros(winL, 1);
y_mvdr_u = zeros(nsample, 1);

for dfr = 1 : nframe
    for df = 1 : nhfft
        Y_mvdr_u(1, df, dfr) = (W_mvdr_u(:, df))' * X(:, df, dfr);
    end
    % ISTFT
    y_tmp2 = ifft([(Y_mvdr_u(1, :, dfr)).'; (Y_mvdr_u(1, end-1:-1:2, dfr))'], 'symmetric');
    y_mvdr_u((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) = y_mvdr_u((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) + y_tmp2(1:winL).*win;
end

% with R_x
y_tmp2 = 0 * y_tmp2;
y_mvdr_x = zeros(nsample, 1);

for dfr = 1 : nframe
    for df = 1 : nhfft
        Y_mvdr_x(1, df, dfr) = (W_mvdr_x(:, df))' * X(:, df, dfr);
    end
    % ISTFT
    y_tmp2 = ifft([(Y_mvdr_x(1, :, dfr)).'; (Y_mvdr_x(1, end-1:-1:2, dfr))'], 'symmetric');
    y_mvdr_x((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) = y_mvdr_x((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) + y_tmp2(1:winL).*win;
end

%%% Spectrogram -----------------------------------------------------------
figure(3); sgtitle('MVDR Beamformer');
subplot(1, 2, 1);
spectrogram(y_mvdr_u, winL, nover, winL, Fs, 'yaxis');
colormap jet; colorbar;
caxis([cmin cmax]);
title('y_{mvdr}(t) with R_{u}');

subplot(1, 2, 2);
spectrogram(y_mvdr_x, winL, nover, winL, Fs, 'yaxis');
colormap jet; colorbar;
caxis([cmin cmax]);
title('y_{mvdr}(t) with R_{x}');

% % Audiowrite
% audiowrite('y_ds.wav', y_ds, Fs);
% audiowrite('y_mvdr_u.wav', y_mvdr_u, Fs);
% audiowrite('y_mvdr_x.wav', y_mvdr_x, Fs);
