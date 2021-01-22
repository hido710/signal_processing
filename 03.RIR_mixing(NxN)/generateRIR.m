function generateRIR(fs, roomDim, RT60, centerSensor, distSensors, distSources, azimuthAngle, elevationAngle)

c = 343;
numDist = length(distSources);
numAzimuth = length(azimuthAngle);
numElevation = length(elevationAngle);
numSensor = size(distSensors,1);
numSource = numDist * numAzimuth * numElevation;

roomDimStr = [num2str(roomDim(1)), 'x', num2str(roomDim(2)), 'x', num2str(roomDim(3))];

rirFolderName = ['RIR_dim_', roomDimStr, '_center_', num2str(centerSensor(1)), '_', num2str(centerSensor(2)), '_', num2str(centerSensor(3)), '_fs_', num2str(fs/1000), 'k/'];
mkdir(rirFolderName);

locationSensor = zeros(numSensor,3);

for sensor = 1 : numSensor
    locationSensor(sensor,:) = centerSensor + distSensors(sensor,:);
end

index = 1;
for dist = 1 : numDist
    for azimuth = 1 : numAzimuth
        for elevation = 1 : numElevation
            
            disp(['[',num2str(index),'/',num2str(numSource),']']);
            
            for sensor = 1 : numSensor
                
                disp(['angle=[', num2str(azimuthAngle(azimuth)), ',', num2str(elevationAngle(elevation)), '] / distance=', num2str(distSources(dist)), ' / sensor:', num2str(sensor)]);
                locationSource = centerSensor + distSources(dist)*[sin(elevationAngle(elevation)/180*pi)*cos(azimuthAngle(azimuth)/180*pi) sin(elevationAngle(elevation)/180*pi)*sin(azimuthAngle(azimuth)/180*pi) cos(elevationAngle(elevation)/180*pi)];
                rirFileName = ['h_RT_', num2str(RT60), '_sensor_', num2str(distSensors(sensor,1)), '_', num2str(distSensors(sensor,2)), '_', num2str(distSensors(sensor,3)), '_dist_' num2str(distSources(dist)), '_angle_', num2str(azimuthAngle(azimuth)),'_', num2str(elevationAngle(elevation)), '.mat'];
                
                if ~exist([rirFolderName, rirFileName],'file')
                    h = rir_generator(c, fs, locationSensor(sensor,:), locationSource, roomDim, RT60);                
                    save([rirFolderName, rirFileName],'h');
                else
                    continue;
                end                
                
            end
            
            index = index + 1;
            
        end
    end
end