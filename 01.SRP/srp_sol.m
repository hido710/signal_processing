%% Readme
% main.m version 2.0.0

%% Clear all
clc;
clear;
close all;

%% Load Wav file
[x_wav,Fs] = audioread('./s_2x8_source_1.wav'); % [-180 120 0 -150 90 30 -120 -90 60 150 -60 -30]
[nsample, M] = size(x_wav);

%% STFT Setting
winL = 512  % number of sample to windowing
nfft = winL
nshift = 128

%%% Window Setting for STFT, stft_inverse.pdf 참고하여 작성-------------%%%
if winL/nshift == 2
    win = sin(pi*([0:1:winL-1]'+0.5)/winL); %1/2 shift
end
if winL/nshift == 4
    win = sqrt(2/3)*hanning(winL,'periodic');%1/4shift
end
%%% -------------------------------------------------------------------%%%
nover = winL - nshift; % number of sample to window overlap
nframe = floor((nsample - nover)/nshift); %number of frame on STFT matrix
nhfft = nfft/2+1;

x_frame = zeros(winL,M);
X = zeros(M,nhfft);
X_buffer = zeros(M,nhfft,nframe);


%% Particle Delay
SS = 343.3   % velocity of sound
anchor_ch = 1  % 기준 sensor의 channel index 

MicDist=zeros(M,3);
load('location_sensor');
MicDist(1,:) = locationSensor{1}'-[1.5 2 1].';
MicDist(2,:) = locationSensor{2}'-[1.5 2 1].';
MicDist(3,:) = locationSensor{3}'-[1.5 2 1].';
MicDist(4,:) = locationSensor{4}'-[1.5 2 1].';

theta = [-180 -150 -120 -90 -60 -30 0 30 60 90 120 150]; % azimuth candidates
phi = 90; % elevation candidate

Lt = length(theta); % number of azimuth candidates 
Lp = length(phi); % number of elevation candidates

%%%------------Plot location sensors and sources 작성-------------------%%%
figure;
hold on
grid minor
for ch = 1:M
    plot3(MicDist(ch,1),MicDist(ch,2),MicDist(ch,3),'bo-')
end

useDist = 1
if useDist == 1
    distsource = 1
end

st_L = zeros(M,nhfft,Lp,Lt);
if useDist == 0
    for it = 1:Lt
        for ip = 1:Lp
            for ch = 1:M
                GapDist = ( MicDist(anchor_ch,:) - MicDist(ch,:) ) * [cos(theta(it)/180*pi)*sin(phi(ip)/180*pi), sin(theta(it)/180*pi)*sin(phi(ip)/180*pi), cos(phi(ip)/180*pi)]';
                for f_idx = 1:nhfft
                    st_L(ch,f_idx,ip,it) = exp(1i*2*pi*(f_idx-1)/nfft*GapDist*Fs/SS); 
                end
            end
        end        
    end
elseif useDist == 1
    for it = 1:Lt
        for ip = 1:Lp
            locationSource = distsource.*[sin(phi(ip)/180*pi)*cos(theta(it)/180*pi) sin(phi(ip)/180*pi)*sin(theta(it)/180*pi) cos(phi(ip)/180*pi)];
            %if (mod(it,10) == 0) && (mod(ip,10)==0)
                plot3(locationSource(1),locationSource(2),locationSource(3),'ro')
            %end
            
            for ch = 1:M                    
                    GapDist = norm(locationSource-MicDist(ch,:), 2); 
                    for f_idx = 1:nhfft
                        st_L(ch,f_idx,ip,it) = exp(1i*2*pi*(f_idx-1)/nfft*GapDist*Fs/SS); 
                    end
            end
        end
       
    end
end


%% Show Spectrogram
%%%--------- STFT구현하기, STFT Guide Task3 참고하여 작성----------------%%%
for t_idx = 1:nframe
    
    %% STFT
    x_shift = x_wav((t_idx-1)*nshift + 1:t_idx*nshift,:);
    for t_idx2 = nshift+1:winL
        x_frame(t_idx2-nshift,:) = x_frame(t_idx2,:); 
    end
    x_frame(nover+1:winL,:) = x_shift; 
    for ch = 1 : M
        Xtmp = fft([x_frame(:,ch).*win ; zeros(nfft - winL,1)]); %frequency resolution만큼 채워주기 위해서 window sample data수에서 추가적인 zero padding작업. 어차피 여기서는 nfft=winL이기 때문에 노상관
        X(ch,:) = Xtmp(1:nhfft).'; 
    end
    X_buffer(:,:,t_idx) = X;
end

figure(4);
hold on
for ch = 1:1 %M
    %subplot(M,1,ch)    
    imagesc(10*log10(abs(squeeze(X_buffer(ch,:,:))).^2));
    colormap jet;
    caxis([-80 0])
    xlabel('Frame Index')
    ylabel('Freq bin')
    title(['STFT of x ' num2str(ch)  'th sensor (nFrame ' num2str(winL) ', nshift ' num2str(nshift) ', nfft ' num2str(nfft) ')'])
    axis xy
    colorbar
    axis([1 nframe 1 nhfft])
end
hold off

%% Frame-by-frame processing
Power_map = zeros(Lp,Lt);
smooth_factor = 0.0
for t_idx = 1:nframe
    
    %%%--------- STFT구현하기, STFT Guide Task3 참고하여 작성------------%%%
   %% STFT
    x_shift = x_wav((t_idx-1)*nshift + 1:t_idx*nshift,:);
    for t_idx2 = nshift+1:winL
        x_frame(t_idx2-nshift,:) = x_frame(t_idx2,:);
    end
    x_frame(nover+1:winL,:) = x_shift;
    for ch = 1 : M
        Xtmp = fft([x_frame(:,ch).*win ; zeros(nfft - winL,1)]);
        X(ch,:) = Xtmp(1:nhfft).';
    end
    X_buffer(:,:,t_idx) = X;
    %%%----------------------------------------------------------------%%%
    
     %% Processing
    inst_Power_map = zeros(Lp,Lt);
    for it = 1:Lt
        for ip = 1:Lp
            wy_all = (st_L(:,:,ip,it)).*(X./(abs(X)+eps));
            wy_all2 = sum(wy_all,1).';            
            inst_Power_map(ip,it) = real(wy_all2'*wy_all2);
        end
    end
    
    Power_map = smooth_factor.*Power_map + (1-smooth_factor).*inst_Power_map;

    %%
    Power_map_norm = Power_map./max(Power_map);
    mean_val = sum(Power_map_norm)/Lt;
    var_val = sum((Power_map_norm-mean_val).^2)/Lt;

    if var_val > 0.01 %max var = 0.25
        sil = 0;
    else
        sil = 1;
    end
    
    if Lp > 1
        [max_az, elevation_ind] = max(Power_map);
        [~ , azimuth] = max(max_az);
        elevation = elevation_ind(azimuth);
        
        azimuthList = theta(azimuth);
        elevationList = phi(elevation);
        
        figure(2);
        imagesc(theta, phi, Power_map);
        xlabel('azimuth(degree)');
        ylabel('elevation(degree)');
        title('Steered Response Power Map')
        colormap jet
        colorbar
        
     else
        [~, az_ind] = max(Power_map);        
        
        azimuthList = theta(az_ind);
        elevationList = phi*ones(length(azimuthList),1);
        
        figure(2);
        %hold on
        plot(theta, Power_map);
        hold on
        plot(azimuthList, Power_map(az_ind), '*');
        hold off
        grid on
        grid minor
        xlabel('azimuth(degree)');
        ylabel('Steered Response Power');
        title('Steered Response Power Map')
        %hold off
        
        figure(3);
        plot(theta, Power_map./max(Power_map),'bo');
        hold on
        plot(azimuthList, Power_map(az_ind)./max(Power_map), '*');
        hold off
        grid on
        grid minor
        xlabel('azimuth(degree)');
        ylabel('Steered Response Power');
        title('Steered Response Power Map')
        ylim([0.5 1])
        xlim([min(theta) max(theta)])
        
    end
    
    figure(4)
    xlim([t_idx t_idx+300])
    title(['sil: ' num2str(sil) ' azimuth' num2str(azimuthList)])
    
    disp(['frame: ' num2str(t_idx)])
    for i = 1:length(azimuthList)
        disp(['sil:' num2str(sil) ' azimuth: ' num2str(azimuthList(i)) ' elevation: ' num2str(elevationList(i))])
    end
    disp([num2str(mean_val) ' ' num2str(var_val)])
    %pause(0.0000001);
end
