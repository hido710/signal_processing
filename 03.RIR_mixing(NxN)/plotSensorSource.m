clc;clear all;close all;

load('location_sensor.mat')
load('location_source.mat')

nch = size(locationSensor,1);
nS = size(locationSource,1);

% centerSensor = zeros(1,3);
% for ch = 1:nch
%     centerSensor = centerSensor + locationSensor{ch};
% end
% centerSensor = centerSensor./nch

figure;
hold on
for ch = 1:nch
    plot3(locationSensor{ch}(1),locationSensor{ch}(2),locationSensor{ch}(3),'bo')
end

for sidx = 1:nS
    plot3(locationSource{sidx}(1),locationSource{sidx}(2),locationSource{sidx}(3),'r*')
end


grid on
grid minor
xlabel('x axis (m)')
ylabel('y axis (m)')
zlabel('z axis (m)')
title('Location of Sensors and Sources')
% xlim([-1 1])
% ylim([-1 1])
% zlim([-1 1])
% zlim([0 1])
hold off