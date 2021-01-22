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

%%% STFT of x(t) ----------------------------------------------------------
x_frame = zeros(winL,nch);
Xtmp = zeros(nch, nhfft);
X = zeros(nch, nhfft, nframe);      % X(t,f)
for dfr = 1 : nframe
    x_frame = x_wav((dfr-1)*nshift+1 : (dfr-1)*nshift+winL, :);
    for ch = 1 : nch
        Xtmp = fft([x_frame(:,ch).* win; zeros(nfft-winL,1)]);
        X(ch,:,dfr) = Xtmp(1:nhfft);
    end
end

%%% STFT of u(t) ----------------------------------------------------------
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
caxis([cmin cmax]);
title('x(t)');

subplot(1, 2, 2);
spectrogram(s_wav(:,8), winL, nover, winL, Fs, 'yaxis');
colormap jet; colorbar;
caxis([cmin cmax]);
title('s(t)');

%% Minimum Variance Distortionless Response Beamformer (MVDR)
%%% Recursive Least-Squares (RLS) 미사용 ----------------------------------
Y_mvdr = zeros(nch, nhfft, nframe);      % Y(t,f) for MVDR beamformer with.R_u
W_mvdr = zeros(nch, nhfft);              % W(f) for MVDR beamformer with.R_u
delta = 1e-3;

% Autocorrelation of noise u(t)
R_tmp = zeros(nch, nch);
R_u = zeros(nch, nch, nhfft);

tic
for df = 1 : nhfft
    for dfr = 1 : nframe
        R_tmp = U(:, df, dfr) * (U(:, df, dfr))';
        R_u(:, :, df) = R_u(:, :, df) + R_tmp;
    end
    R_u(:, :, df) = 1/nframe*(R_u(:, :, df)) + delta*eye(nch);
end

% Filter Matrix
for df = 1 : nhfft
   W_mvdr(:, df) = (R_u(:,:,df)\H(:,df)) ./ (H(:,df)'*(R_u(:,:,df)\H(:,df)));
end

% Beamformer Output
y_tmp1 = zeros(winL, 1);
y_mvdr = zeros(nsample, 1);
for dfr = 1 : nframe
    for df = 1 : nhfft
        Y_mvdr(1, df, dfr) = (W_mvdr(:, df))' * X(:, df, dfr);
    end
    % ISTFT
    y_tmp1 = ifft([(Y_mvdr(1, :, dfr)).'; (Y_mvdr(1, end-1:-1:2, dfr))'], 'symmetric');
    y_mvdr((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) = y_mvdr((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) + y_tmp1(1:winL).*win;
end
fprintf('RLS 미사용 시, ');
toc

%%% Recursive Least-Squares (RLS) 사용 ------------------------------------
Y_mvdr_rls = zeros(nch, nhfft, nframe);   % Y(t,f) for MVDR beamformer with.R_u
W_mvdr_rls = zeros(nch, nhfft);           % W(f) for MVDR beamformer with.R_u

% Autocorrelation of noise u(t) _ Inversion
R_u_inv = zeros(nch, nch, nhfft); % nframe 행 만들어서 따로 저장하지 않고 update만 함
for df = 1 : nhfft
    R_u_inv(:,:,df) = eye(nch);   % 초기 값
end
tmp1 = zeros(nch, nch);
gamma = 0.985; % forgetting factor

% Update R_u, W and Y
tic
for dfr = 1 : nframe
    for df = 1 : nhfft
        tmp1 = gamma^(-2) * R_u_inv(:,:,df) * U(:,df,dfr) * (U(:,df,dfr))' * R_u_inv(:,:,df);
        tmp2 = 1 + gamma^(-1) * (U(:,df,dfr))' * R_u_inv(:,:,df) * U(:,df,dfr);
        
        R_u_inv(:,:,df) = gamma^(-1) * R_u_inv(:,:,df) - tmp1/tmp2;
        W_mvdr_rls(:, df) = (R_u_inv(:,:,df) * H(:,df)) ./ (H(:,df)'*R_u_inv(:,:,df)*H(:,df));
        Y_mvdr_rls(1, df, dfr) = (W_mvdr_rls(:, df))' * X(:, df, dfr);
    end
end

% Beamformer Output
y_tmp2 = zeros(winL, 1);
y_mvdr_rls = zeros(nsample, 1);

for dfr = 1 : nframe
    for df = 1 : nhfft
        Y_mvdr_rls(1, df, dfr) = (W_mvdr_rls(:, df))' * X(:, df, dfr);
    end
    % ISTFT
    y_tmp2 = ifft([(Y_mvdr_rls(1, :, dfr)).'; (Y_mvdr_rls(1, end-1:-1:2, dfr))'], 'symmetric');
    y_mvdr_rls((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) = y_mvdr_rls((dfr-1)*nshift+1 : (dfr-1)*nshift+winL) + y_tmp2(1:winL).*win;
end
fprintf('\nRLS 사용 시, ');
toc

%%% Spectrogram -----------------------------------------------------------
figure(2); sgtitle('MVDR Beamformer');

subplot(1, 2, 1);
spectrogram(y_mvdr, winL, nover, winL, Fs, 'yaxis');
colormap jet; colorbar;
caxis([cmin cmax]);
title('RLS X');

subplot(1, 2, 2);
spectrogram(y_mvdr_rls, winL, nover, winL, Fs, 'yaxis');
colormap jet; colorbar;
caxis([cmin cmax]);
title('RLS O');
