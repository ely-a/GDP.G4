% Main script for Neptune mission trajectory analysis
% Clear workspace and close all figures
clear all;
close all;
clc;

% Constants
c = 299792.458;       % Speed of light (km/s)
f0 = 8.4e9;           % X-band frequency (Hz)
crit_angle = 3;       % Solar conjunction critical angle (degrees)
min_elevation = 10;   % Minimum elevation for DSN visibility (degrees)

% Load spacecraft trajectory data
sc_data = load('sc_data.mat');
sc_pos = sc_data.r_out; % Nx3 array in ICRF [km]
planet_data = load('planet_data.mat');
earth_vel = planet_data.earth_vel;
earth_pos = planet_data.earth_pos;
jupiter_pos = planet_data.jupiter_pos;
neptune_pos = planet_data.neptune_pos;
sun_pos = planet_data.sun_pos;
N = size(sc_pos, 1);

% Time setup (departure: 23-Mar-2032, arrival: 01-Jun-2040)
t_start = datetime('2032-03-23');
t_end = datetime('2040-06-01');
t_flyby = datetime('2033-07-01');
time_vec = t_start + days(linspace(0, days(t_end - t_start), N))';
time_days = days(time_vec - t_start); % Time since departure [days]
time_flyby = days(t_flyby - t_start);
total_seconds = seconds(t_end - t_start);
dt_sec = total_seconds / (N-1); % Time step [s]

% Compute spacecraft velocity using finite differences
vel_sc = zeros(N, 3);
for i = 1:N
    if i == 1
        vel_sc(i,:) = (sc_pos(i+1,:) - sc_pos(i,:)) / dt_sec;
    elseif i == N
        vel_sc(i,:) = (sc_pos(i,:) - sc_pos(i-1,:)) / dt_sec;
    else
        vel_sc(i,:) = (sc_pos(i+1,:) - sc_pos(i-1,:)) / (2*dt_sec);
    end
end

% Vector from Earth to SC
r_earth_sc = sc_pos - earth_pos; % Nx3
dist_earth_sc = vecnorm(r_earth_sc, 2, 2); % Nx1

% Task 1: Doppler Tracking Analysis
los = r_earth_sc ./ dist_earth_sc; % Line-of-sight unit vector
rel_vel = vel_sc - earth_vel;      % Relative velocity [km/s]
radial_vel = dot(rel_vel, los, 2); % Radial component [km/s]
frequency = f0 * (1 - radial_vel / c); % Received frequency [Hz]

% Task 2: Solar Conjunction Analysis
r_earth_sun = sun_pos - earth_pos; % Vector from Earth to Sun
sep_angle = acosd(dot(r_earth_sun, r_earth_sc, 2) ./ ...
         (vecnorm(r_earth_sun, 2, 2) .* dist_earth_sc)); % SEP angle [deg]
conjunction_binary = sep_angle < crit_angle; % 1 = conjunction, 0 = no conjunction

% Task 3: DSN Visibility Analysis

% Precompute Moon positions
moon_pos = zeros(N, 3);
jd = juliandate(time_vec);
for i = 1:N
    [moon_pos(i,:), ~] = planetEphemeris(jd(i), 'Sun', 'Moon', '430');
end

dsn_stations = [
    35.426667, -116.89, 0.8;   % Goldstone
    40.433333, -4.25, 0.7;      % Madrid
    -35.401389, 148.981667, 0.6 % Canberra
];
dsn_binary = true(N,1); % Initialize as not visible
earth_ellipsoid = wgs84Ellipsoid('km'); % Use WGS84 ellipsoid

for i = 1:N
    % Convert SC position to ECEF
    dcm = dcmeci2ecef('IAU-2000/2006', time_vec(i));
    sc_ecef = dcm * r_earth_sc(i,:)'; % SC ECEF position (km)
    
    for station = 1:size(dsn_stations,1)
        % Get station LLA
        station_lla = dsn_stations(station,:);
        
        % Compute elevation using ECEF coordinates
        [~, elev_deg] = ecef2aer(...
            sc_ecef(1), sc_ecef(2), sc_ecef(3), ... % SC ECEF position
            station_lla(1), station_lla(2), station_lla(3), ... % Station LLA
            earth_ellipsoid);
        
        if elev_deg >= min_elevation
            dsn_binary(i) = false; % Visible
            break; % No need to check other stations
        end
    end
end

% Task 4: Light-Time Delay
light_time = dist_earth_sc / c; % One-way signal travel time [s]

% Sense Check: Minimum distances to planets
dist_to_earth = vecnorm(sc_pos - earth_pos, 2, 2);
dist_to_jupiter = vecnorm(sc_pos - jupiter_pos, 2, 2);
dist_to_neptune = vecnorm(sc_pos - neptune_pos, 2, 2);
min_dist = [min(dist_to_earth), min(dist_to_jupiter), min(dist_to_neptune)];
fprintf('Minimum distance to Earth: %.2f km\n', min_dist(1));
fprintf('Minimum distance to Jupiter: %.2f km\n', min_dist(2));
fprintf('Minimum distance to Neptune: %.2f km\n', min_dist(3));

% Generate all figures
% Figure 1: Radial Velocity vs. Time
figure('Name', 'Radial Velocity', 'NumberTitle', 'off');
plot(time_days, radial_vel, 'b-', 'LineWidth', 1.5);
xline(time_flyby, '--r', 'Jupiter flyby');
xlabel('Time Since Departure (days)');
ylabel('Radial Velocity (km/s)');
title('Doppler Tracking: Radial Velocity');
grid on;

% Figure 2: Frequency vs. Time
figure('Name', 'Frequency', 'NumberTitle', 'off');
plot(time_days, frequency, 'r-', 'LineWidth', 1.5);
xline(time_flyby, '--b', 'Jupiter flyby');
xlabel('Time Since Departure (days)');
ylabel('Frequency (Hz)');
title('Doppler Tracking: Received Frequency');
grid on;

% Figure 3: SEP Angle vs. Time
figure('Name', 'SEP Angle', 'NumberTitle', 'off');
plot(time_days, sep_angle, 'm-', 'LineWidth', 1.5);
xlabel('Time Since Departure (days)');
ylabel('SEP Angle (degrees)');
title('Solar Conjunction: Sun-Earth-Probe Angle');
grid on;
hold on;
xline(time_flyby, '--r', 'Jupiter flyby');
yline(3, '--r', 'Critical Angle (3°)');
hold off;

% Figure 4: Solar Conjunction Binary
figure('Name', 'Solar Conjunction', 'NumberTitle', 'off');
stairs(time_days, conjunction_binary, 'LineWidth', 1.5, 'Color', [0, 0.5, 0]);
xline(time_flyby, '--r', 'Jupiter flyby');
ylim([-0.1, 1.1]);
xlabel('Time Since Departure (days)');
ylabel('Conjunction Status (1: Yes, 0: No)');
title('Solar Conjunction: SEP < 3°');
grid on;

% Figure 5: DSN Visibility Binary
figure('Name', 'DSN Visibility', 'NumberTitle', 'off');
stairs(time_days, dsn_binary, 'LineWidth', 1.5, 'Color', [0.7, 0, 0]);
xline(time_flyby, '--b', 'Jupiter flyby');
ylim([-0.1, 1.1]);
xlabel('Time Since Departure (days)');
ylabel('Visibility Status (1: Not Visible, 0: Visible)');
title('DSN Visibility: No Station Visible');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elevation analysis
elevation_history = zeros(N,1);
for i = 1:10:N  % Sample every 100th point
    dcm = dcmeci2ecef('IAU-2000/2006', time_vec(i));
    sc_ecef = dcm * r_earth_sc(i,:)';
    for station = 1:size(dsn_stations,1)
        station_lla = dsn_stations(station,:);
        [~, elev_deg] = ecef2aer(...
            sc_ecef(1), sc_ecef(2), sc_ecef(3), ...
            station_lla(1), station_lla(2), station_lla(3), ...
            earth_ellipsoid);
        elevation_history(i) = max(elevation_history(i), elev_deg);
    end
end

% Plot elevation
figure('Name', 'Elevation Analysis', 'NumberTitle', 'off');
plot(time_days, elevation_history, 'b-');
xlabel('Time Since Departure (days)');
ylabel('Maximum Elevation (degrees)');
title('DSN Visibility: Maximum Station Elevation');
yline(min_elevation, 'r--', 'Visibility Threshold');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 6: Light-Time Delay
figure('Name', 'Light-Time Delay', 'NumberTitle', 'off');
plot(time_days, light_time, 'k-', 'LineWidth', 1.5);
xline(time_flyby, '--r', 'Jupiter flyby');
xlabel('Time Since Departure (days)');
ylabel('Signal Travel Time (seconds)');
title('Light-Time Delay to Earth');
grid on;