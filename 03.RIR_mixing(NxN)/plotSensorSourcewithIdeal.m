%% Initialize
clc;
clear;
close all;

%% time delay 구하기

%%% Setting ---------------------------------------------------------------
SS = 343.3;
rad = pi/180;
load('location_sensor');
load('location_source'); % Ideal Source 위치

mic_anchor = [1.5 2 1];    % mic array 중심
L = 19;                    % -90~90 10도 간격
d_azimuth = 10;            % 10도 간격

%%% virtual source location ----------------------------------------------
s_tmp = zeros(L, 3);

for dl = 1 : L             % [-90 -80 ... 80 90]
    s_tmp(dl, :) = [sin((dl-1)*d_azimuth*rad) -cos((dl-1)*d_azimuth*rad) 0] + mic_anchor;
end

figure(); hold on;
for ch = 1 : length(locationSensor)
    plot3(locationSensor{ch}(1), locationSensor{ch}(2), locationSensor{ch}(3), 'bo')
end
for i = 1 : length(locationSource)
    plot3(locationSource{i}(1), locationSource{i}(2), locationSource{i}(3), 'r*')
end
for dl = 1 : L
    plot3(s_tmp(dl,1), s_tmp(dl,2), s_tmp(dl,3), 'go')
end
xlabel('x axis (m)'); ylabel('y axis (m)'); zlabel('z axis (m)');
axis([1.5 3.5 1 3]);
title('Location of Sensor and Virtual Sources');
grid; hold off;