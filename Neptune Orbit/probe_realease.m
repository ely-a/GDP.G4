clear; clc; close all;

if exist('ScienceOrbit.mat', 'file') == 2
    load("ScienceOrbit.mat")
    mu_N = 6.8351e6; % km^3/s^2
    r_N = 24622;             % km (mean radius)
else
    error('No initial conditions');
end

fprintf('%6s %10s %10s %12s\n', 'Theta', 'Altitude', 'Latitude', 'Longitude');

for theta = 0:5:360
    % Position and velocity from orbital elements
    [r, v] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, theta, mu_N);
    
    % Magnitude of position vector (from Neptune's center)
    r_mag = norm(r);
    
    % Altitude above Neptune's surface
    altitude = r_mag - r_N;

    % Compute planetocentric latitude and longitude (from ECI)
    x = r(1); y = r(2); z = r(3);
    lat = asind(z / r_mag);              % latitude in degrees
    lon = atan2d(y, x);                  % longitude in degrees

    % Normalize longitude to [-180, 180]
    if lon > 180
        lon = lon - 360;
    end

    % Output
    fprintf('%6.0f %10.2f km %9.2f° %10.2f°\n', theta, altitude, lat, lon);
end
