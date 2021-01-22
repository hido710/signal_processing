%% Readme
% STFT_batch version 1.0.0
% pair with ISTFT_batch version 1.0.0
%% Input
% x: [samples x channels]
% - time domain amplitude

% winL: number of sample to windowing
% nfft: number of fft point
% nshift: number of sample to window shift
%% Output
% X: [channels x nhfft x frames]
% - STFT result
%%
function X = STFT_batch(x,winL,nfft,nshift)

shiftdiv = winL/nshift;

if shiftdiv == 2
    win = sin(pi*([0:1:winL-1]'+0.5)/winL); %1/2 shift
end
if shiftdiv == 4
    win = sqrt(2/3)*hanning(winL,'periodic');%1/4shift
end

nover = winL - nshift; % number of sample to window overlap
nhfft = nfft/2+1; % half +1 of fftsize

nch = size(x,2);
nframe = floor((length(x(:,1)) - nover)/nshift); %number of frame on STFT matrix
X = zeros(nch,nhfft,nframe); %STFT memory allocate

for ch = 1 : nch
    for t_idx = 1 : nframe
        Xtmp = fft([x((t_idx-1)*nshift + 1 : (t_idx-1)*nshift + winL, ch).*win ; zeros(nfft - winL,1)]);
        X(ch,:,t_idx) = Xtmp(1:nhfft);
    end
end

end