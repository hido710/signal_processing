function H = mixingNxN(mixoutRoot, signalFileName, SNR, RT60, distSensors, distSources, azimuthAngle, elevationAngle, targetFs, upFs, roomDim, centerSensor)

% NxN mixing
disp(['SNR = ',num2str(SNR),' (dB),  RT60 = ',num2str(RT60),' (sec)']);

numSensor = size(distSensors,1);
numSource = size(distSources,2);

disp(['Number of Source = ' num2str(numSource), ', Number of Sensor = ' num2str(numSensor), ...
    ', SNR = ',num2str(SNR),' (dB),  RT60 = ',num2str(RT60),' (sec)']);

% GENERATING ROOM IMPULSE RESPONSE
k = 0.161;
[H, locationSensor, locationSource] = genNxN(upFs, roomDim, RT60, k, centerSensor, distSensors, distSources, azimuthAngle, elevationAngle, 12);

mkdir(mixoutRoot);
save([mixoutRoot, '/location_sensor.mat'],  'locationSensor');
save([mixoutRoot, '/location_source.mat'],  'locationSource');
save([mixoutRoot, '/angle_source.mat'], 'azimuthAngle', 'elevationAngle');

Signal = cell(numSource,1);
sUp = cell(numSensor, numSource);
s = cell(numSensor, numSource);

% PROCESSING
for i = 1 : numSource
    [origSrc, fsOrig] = audioread(signalFileName{i});        
    Signal{i} = resample(origSrc, upFs, fsOrig);
    if i == 1
        minLength = length(Signal{i});
        origScale = max(abs(origSrc));
    else
        minLength = min(minLength, length(Signal{i}));
        origScale = max(origScale, max(abs(origSrc)));
    end
end  


for i = 1 : numSource
    Signal{i} = Signal{i}(1:minLength);
end   

% FILTERING EACH CHANNELS
if RT60==0    
    for i = 1 : numSensor
        for j = 1 : numSource
            [filtCoeff, lagIndex] = max(H(:,i,j));
            sUp{i,j} = [zeros(lagIndex-1,1);Signal{j}(1:length(Signal{j})-lagIndex+1).*filtCoeff];
        end
    end 
else
    for i = 1 : numSensor
        for j = 1 : numSource
            sUp{i,j} = fft_filter(H(:,i,j), 1, Signal{j});
        end
    end
end

% NORMALIZING SOUND LEVEL
for i = 1 : numSensor
    for j = 1 : numSource        
        s{i,j} = resample(sUp{i,j},targetFs,upFs);

	    if i == 1 && j==1
            maxScale = max(abs(s{i, j}));
        else
            maxScale = max(maxScale, max(abs(s{i, j})));
        end
    end
end

for i = 1 : numSensor        
    for j = 1 : numSource
	    s{i, j} = s{i, j} .*origScale ./ maxScale;        
    end   
end

% SCALING NOISE FOR DESIRED SNR % MIXING CHANNELS
interference = sum(cell2mat(s(1,2:end)),2);
energyTarget = sqrt(sum(s{1,1}.^2));
energyInterference = sqrt(sum(interference.^2));
energyNormal = energyTarget/energyInterference;
x = zeros(length(s{1,1}), numSensor);

for i = 1: numSensor
    target = s{i,1};
    interference = sum(cell2mat(s(i,2:end)),2);    
    interference = interference*energyNormal/sqrt(10^(SNR/10));
    x(:, i) = target+interference;
    
    if i == 1
        maxScale = max(abs(x(:, i)));
    else
        maxScale = max(maxScale, max(abs(x(:, i))));
    end
end

% NORMALIZING SOUND LEVEL
for i = 1 : numSensor
    x(:,i) = x(:,i)./maxScale;
    s{i,1} = s{i,1}./maxScale;
    for j = 2 : numSource
        s{i,j} = s{i,j}./maxScale * energyNormal/sqrt(10^(SNR/10));
    end
end

% WRITING OUTPUT SIGNAL  
audiowrite([mixoutRoot, '/x_', num2str(numSource), 'x', num2str(numSensor), '.wav'], x, targetFs);
sMat = s(:,1).';
audiowrite([mixoutRoot, '/s_', num2str(numSource), 'x', num2str(numSensor) , '_source_', num2str(1),'.wav'], cell2mat(sMat), targetFs);
for j = 2 : numSource
    sMat = s(:,j).';
%     audiowrite([mixoutRoot, '/s_', num2str(numSource), 'x', num2str(numSensor) , '_source_', num2str(j),'.wav'], cell2mat(sMat)*energyNormal/sqrt(10^(SNR/10)), targetFs);    
    audiowrite([mixoutRoot, '/s_', num2str(numSource), 'x', num2str(numSensor) , '_source_', num2str(j),'.wav'], cell2mat(sMat), targetFs);    
end

