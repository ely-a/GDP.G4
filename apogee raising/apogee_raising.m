clc
clear
close all


% Constants
mu = 398600.4418;    % Earth's gravitational parameter [km^3/s^2]
R_earth = 6378.0;    % Earth's radius [km]

%Things
h0 = 300;                       % Initial altitude of apogee (circular orbit)
ra = R_earth + h0;              % Iniital radius of spacecraft
rp = 300 + R_earth;             % Perigee in km
v_c = sqrt(mu/r0);              % Velocity at initial ciruclar orbit
a = (rp + ra)/2;                % Semi-major axis in km   
e = (ra - rp) ./ (ra + rp);     % Eccentricity
P = 2*pi * sqrt(a^3 / mu);      % Period of orbit in seconds
P_days = P/(24*3600);           % Period of orbit in days


max_days = 21;      % Maximum amount of days orbiting Earth



