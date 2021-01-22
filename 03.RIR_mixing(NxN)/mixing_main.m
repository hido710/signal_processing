clear;
close all;
clc;

%% preparation

upFs                = 1024000;
distSourceRange     = [1 1];
azimuthAngleRange   = -90:10:90;
elevationAngleRange = [90];


roomDim       = [5 4 3];  % 단위 : (m)
centerSensors = [1.5 2 1];

RT60 = 0.0;


distSensors(1,:) = [0 0.14 0];
distSensors(2,:) = [0 0.10 0];
distSensors(3,:) = [0 0.06 0];
distSensors(4,:) = [0 0.02 0];
distSensors(5,:) = [0 -0.02 0];
distSensors(6,:) = [0 -0.06 0];
distSensors(7,:) = [0 -0.10 0];
distSensors(8,:) = [0 -0.14 0];


generateRIR(upFs, roomDim, RT60, centerSensors, distSensors, distSourceRange, azimuthAngleRange, elevationAngleRange);


%% mixing

numSensor = 8;
numSource = 2;
SNR = 0;
targetFs = 16000;

signalFileName{1} = './data/clean/speech_female1.wav';
signalFileName{2} = './data/clean/white_16k_30sec.wav';

distSources = [1 1];

azimuthAngle  = [30 -80];
elevationAngle  = [90 90];

if (numSource == size(distSources,2) && numSource == size(azimuthAngle,2) && numSource == size(elevationAngle,2) && numSource == length(signalFileName))
    disp('numSource correct')
else
    disp('numSource not correct. Check numSource,distSources,azimuthAngle,elevationAngle')
    pause
end

if(numSensor == size(distSensors,1))
    disp('numSensor correct')
else
    disp('numSensor not correct. Check distSensors, distSensors')
    pause
end



mixoutRoot = [ 'data/mixture_', num2str(numSource), 'x', num2str(numSensor) , '_RT60_' num2str(RT60), '_SNR_', num2str(SNR)];
disp('execute mixingNxN');
mixingNxN(mixoutRoot, signalFileName, SNR, RT60, distSensors, distSources, azimuthAngle, elevationAngle, targetFs, upFs, roomDim, centerSensors);

disp('done')

