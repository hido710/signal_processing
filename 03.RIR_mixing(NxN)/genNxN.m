function [H, locationSensor, locationSource]=genNxN(fs, roomDim, RT60, k, centerSensor, distSensors, distSources, azimuthAngle, elevationAngle, n)

% [[[ Example ]]]
%
% H=gen_2x2(44100,[5 4 3],0.5,[2 1.5 1.1],0.08,1,1,0,30,12)
%
% fs=44100;
% room_dim = [5 4 3];
% RT60 = 0.5;
% k = 0.161;
% center_sensors = [2 1.5 1.1];
% dist_sensors = .08;
% dist_source1 = 1;
% dist_source2 = 1;
% angle_source1 = 0;
% angle_source2 = 30;
% n = 12;

c = 343;
v = roomDim(1) * roomDim(2) * roomDim(3);
s = ( roomDim(1)*roomDim(2) + roomDim(2)*roomDim(3) + roomDim(3)*roomDim(1) ) * 2;
if RT60==0
    R = 0;
else
    R = 1-(k*v/s/RT60);
end
%=======================
%R = .66; % RT60 = .3 sec
%R = .79; % RT60 = .5 sec
%R = .87; % RT60 = .8 sec

numSensor = size(distSensors,1);
numSource = size(distSources,2);
hTotal = cell(numSensor, numSource);
attenu = zeros(numSensor, numSource);
locationSource = cell(numSource,1);
locationSensor = cell(numSensor,1);

roomDimStr = [num2str(roomDim(1)), 'x', num2str(roomDim(2)), 'x', num2str(roomDim(3))];
rirFolderName = ['RIR_dim_', roomDimStr, '_center_', num2str(centerSensor(1)), '_', num2str(centerSensor(2)), '_', num2str(centerSensor(3)), '_fs_', num2str(fs/1000), 'k/'];


for sensor = 1 : numSensor
    for source = 1 : numSource
        locationSource{source} = centerSensor + distSources(source)*[sin(elevationAngle(source)/180*pi)*cos(azimuthAngle(source)/180*pi) sin(elevationAngle(source)/180*pi)*sin(azimuthAngle(source)/180*pi) cos(elevationAngle(source)/180*pi)];
        locationSensor{sensor} = centerSensor + distSensors(sensor,:);
        
        rirFileName = ['h_RT_', num2str(RT60), '_sensor_', num2str(distSensors(sensor,1)), '_', num2str(distSensors(sensor,2)), '_', num2str(distSensors(sensor,3)), '_dist_' num2str(distSources(source)), '_angle_', num2str(azimuthAngle(source)),'_', num2str(elevationAngle(source)), '.mat'];
        
        if RT60 ~= 0        
            if exist([rirFolderName, rirFileName],'file')
                hTmp = load([rirFolderName, rirFileName]);
                hTmp = hTmp.h;
            else
                hTmp = rir_generator(c, fs, locationSensor{sensor}, locationSource{source}, roomDim, RT60);            
                h = hTmp;
                if ~exist(rirFolderName, 'dir')
                    mkdir(rirFolderName);
                end
                save([rirFolderName, rirFileName],'h');
            end
        else
            hTmp = rir(fs,locationSensor{sensor},n,0,roomDim,locationSource{source});            
        end
        hTotal{sensor,source} = hTmp;        
        
        if sensor == 1 && source == 1
            minLength = length(hTotal{sensor,source});
        else
            minLength = min(minLength, length(hTotal{sensor,source}));
        end
    end
end

H = zeros(minLength,numSensor,numSource);

for i = 1 : numSensor
    for j = 1 : numSource
        H(:,i,j) = hTotal{i,j}(1:minLength);
    end
end

[maxValue, ~] = max(H);
H = H/maxValue(1,1,1);

if RT60 == 0.1
    Htmp = H;
    clear H;
    H = zeros(minLength, numSensor, numSource);
    for i = 1: numSensor
        for j = 1 : numSource
            [~, pos] = max(Htmp(:, i, j));
            H(pos, i, j) = 1;
        end
    end
end


for i = 1 : numSensor
    for j = 1 : numSource
        attenu(i, j) = max(H(:, i, j));
    end
end