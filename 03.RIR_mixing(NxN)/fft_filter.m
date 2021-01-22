function [y]=fft_filter(b,a,X)

% [ ���ļ� �����ο��� �����ϴ� Filtering �Լ�]
% ===================================================================================
% [�Է� �μ�] 
%
% 1. b  : numerator coefficient vector
% 2. a  : denominator coefficient, ���⼱ �׳� 1�� ��(���� ���� ����)
% 3. X  : input data in vector


% [��� �μ�]
% 1. y  : Output data in vector

if size(b,1) < size(b,2)    
    b=b';
end

Size_X =size(X);
if size(X,1) < size(X,2)
    X=X';
end

length_X = length(X);  % Input data ����
length_b = length(b);  % b�� ����

% FFT point ����
for i=1:10000,
	if length_b <= 2^(i+1) & length_b > 2^i
		nfft = 2^(i+2);                     % Tap������ 2�� ���� ����
		break;
	end
end

% b�� Time -> Freq domain�� ��ȯ
freq_b = fft(b(:,1),nfft);

% Error check
if i==10000,
	error('The length of the filter is too long!!!');
end

%%% ���� �� ���� ����

pad_length = nfft-length_b;    % padding size
pad = zeros(pad_length,1);     % FFT ��ȯ �� zero padding(overlap save)
x_block  = zeros(nfft,1);      % �Է� data X�� STFT ���� ���� ����

%%% Mixing ����
% STFT�� ¥���� �κб��� ó���ϰ��� �Է½�ȣ X�� ���̸� Ȯ��(zero padding)
length_tmp_mix = length_X+pad_length;
tmp_X = zeros(length_tmp_mix,1);        % ���̰� Ȯ��� �Է� ��ȣ
tmp_X(1:length_X,1) = X(1:length_X,1);  % ������ �Է� X�κ��� �״�� �Է�, �������� zero padding
tmp_mix = zeros(length_tmp_mix,1);      % ���̰� Ȯ��� ȥ�ս�ȣ

%Mixing with room impulse response
for n=1:floor(length_X/length_b)               % overlap save 
    x_block = (fft([tmp_X(length_b*(n-1)+1:length_b*n); pad]));
    tmp_mix(length_b*(n-1)+1:length_b*n+pad_length) = tmp_mix(length_b*(n-1)+1:length_b*n+pad_length) + real(ifft(x_block.*freq_b));
end

% ������ ���̿� �°� ����
y = tmp_mix(1:length_X);
if Size_X(1) < Size_X(2)    %������ �������� ���缭 return
    y=y';
end
















       
            

         
