%% Readme
% ISTFT_batch version 1.0.0
% pair with STFT_batch version 1.0.0
%% Input
% X: [channels x nhfft x frames]
% - STFT result
% winL: number of sample to windowing
% nshift: number of sample to window shift
% nsamples: mic signal's time domain length
%% Output
% x: [nsamples x nch]
% ISTFT result

%%
function x = ISTFT_batch(X,winL,nshift,nsamples)

[nch, nhfft, nframe] = size(X);
x = zeros(nsamples,nch);
shiftdiv = winL/nshift;
nfft = (nhfft-1)*2;
if shiftdiv == 2
    win = sin(pi*([0:1:winL-1]'+0.5)/winL); %1/2 shift
end
if shiftdiv == 4
    win = sqrt(2/3)*hanning(nfft,'periodic');%1/4shift
end


for ch = 1:nch
    for t_idx = 1 : nframe
        %xTmp = real(ifft([transpose(X(ch,:,t_idx)) ; (X(ch,end-1:-1:2,t_idx)')]));
        xTmp = ifft([(X(ch,:,t_idx).'); (X(ch,end-1:-1:2,t_idx)')],'symmetric');
        x((t_idx-1)*nshift + 1 : (t_idx-1)*nshift + winL, ch) = x((t_idx-1)*nshift + 1 : (t_idx-1)*nshift + winL, ch) + xTmp(1:winL).*win;
        
    end
end
end