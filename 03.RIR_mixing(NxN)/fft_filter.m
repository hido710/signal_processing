function [y]=fft_filter(b,a,X)

% [ 주파수 도메인에서 수행하는 Filtering 함수]
% ===================================================================================
% [입력 인수] 
%
% 1. b  : numerator coefficient vector
% 2. a  : denominator coefficient, 여기선 그냥 1로 고(추후 보완 예정)
% 3. X  : input data in vector


% [출력 인수]
% 1. y  : Output data in vector

if size(b,1) < size(b,2)    
    b=b';
end

Size_X =size(X);
if size(X,1) < size(X,2)
    X=X';
end

length_X = length(X);  % Input data 길이
length_b = length(b);  % b의 길이

% FFT point 결정
for i=1:10000,
	if length_b <= 2^(i+1) & length_b > 2^i
		nfft = 2^(i+2);                     % Tap수보다 2배 많게 결정
		break;
	end
end

% b를 Time -> Freq domain로 변환
freq_b = fft(b(:,1),nfft);

% Error check
if i==10000,
	error('The length of the filter is too long!!!');
end

%%% 변수 및 벡터 설정

pad_length = nfft-length_b;    % padding size
pad = zeros(pad_length,1);     % FFT 변환 시 zero padding(overlap save)
x_block  = zeros(nfft,1);      % 입력 data X의 STFT 사용될 벡터 설정

%%% Mixing 관련
% STFT시 짜투리 부분까지 처리하고자 입력신호 X의 길이를 확장(zero padding)
length_tmp_mix = length_X+pad_length;
tmp_X = zeros(length_tmp_mix,1);        % 길이가 확장된 입력 신호
tmp_X(1:length_X,1) = X(1:length_X,1);  % 원래의 입력 X부분은 그대로 입력, 나머지는 zero padding
tmp_mix = zeros(length_tmp_mix,1);      % 길이가 확장된 혼합신호

%Mixing with room impulse response
for n=1:floor(length_X/length_b)               % overlap save 
    x_block = (fft([tmp_X(length_b*(n-1)+1:length_b*n); pad]));
    tmp_mix(length_b*(n-1)+1:length_b*n+pad_length) = tmp_mix(length_b*(n-1)+1:length_b*n+pad_length) + real(ifft(x_block.*freq_b));
end

% 원래의 길이에 맞게 조정
y = tmp_mix(1:length_X);
if Size_X(1) < Size_X(2)    %원래의 사이즈대로 맞춰서 return
    y=y';
end
















       
            

         
