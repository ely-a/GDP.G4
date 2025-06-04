clc 
clear
close all

% Inputs
v = 17.1462;             % Velocity magnitude (e.g., km/s)
RA_deg = 27.2887;        % Right Ascension in degrees
Dec_deg = -1.5998;       % Declination in degrees

% Convert degrees to radians
RA = deg2rad(RA_deg);
Dec = deg2rad(Dec_deg);

% Compute velocity vector components
vx = v * cos(Dec) * cos(RA);
vy = v * cos(Dec) * sin(RA);
vz = v * sin(Dec);

% Combine into a vector
velocity_vector = [vx vy vz];

% Display result
disp('Velocity vector (vx, vy, vz):');
disp(velocity_vector);
