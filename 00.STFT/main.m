%% Readme
% main.m version 1.0.0
% STFT-ISTFT batch
%% Clear all

% clc;
clear all;
% close all;

%% Load Wav file
[x,Fs] = audioread('s_3x8_source_1.wav');
%x = x(1:Fs,1);
[nsample, nch] = size(x);

%% STFT batch
winL = 512
nfft = winL
nshift = 128
X = STFT_batch(x,winL,nfft,nshift);
nFrame = size(X,3);

%% Show Spectrogram
for ch = 1:nch
    figure;
    %imagesc(10*log10(abs(squeeze(X(ch,:,:))).^2));
    %xlabel('Frame Index')
    %ylabel('Freq bin')
    
    imagesc(linspace(0,nshift*nFrame/Fs,nFrame),linspace(0,Fs/2,nfft/2+1),10*log10(abs(squeeze(X(ch,:,:))).^2));
    xlabel('sec')    
    ylabel('hz')
    
    
    colormap jet;
    caxis([-110 25])
    title(['STFT of x ' num2str(ch)  'th sensor (nFrame ' num2str(winL) ', nshift ' num2str(nshift) ', nfft ' num2str(nfft) ')'])
    axis xy
    
end

%% Process
Y = X;

%% ISTFT batch
y = ISTFT_batch(Y,winL,nshift,length(x));

%%
for ch = 1:nch
    figure;
    plot(x(:,ch)-y(:,ch))
end


