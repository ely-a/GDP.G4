clear
clc
close all

mu_sun = 1.32712440018e11; % km^3/s^2 (gravitational parameter of the Sun)
r1 = 1.496e8;              % km (1 AU = Earth's orbit)
r2 = 4.498e9;              % km (30.07 AU = Neptune's orbit)

% Semi-major axis of transfer ellipse
a_t = (r1 + r2)/2;

% Velocities
v1 = sqrt(mu_sun / r1);             % Earth orbital speed
v2 = sqrt(mu_sun / r2);             % Neptune orbital speed

v_trans1 = sqrt(mu_sun * (2/r1 - 1/a_t));  % velocity at perihelion (departure)
v_trans2 = sqrt(mu_sun * (2/r2 - 1/a_t));  % velocity at aphelion (arrival)

% Delta-Vs
deltaV1 = v_trans1 - v1;   % at Earth
deltaV2 = v2 - v_trans2;   % at Neptune

% Total Delta-V
deltaV_total = abs(deltaV1) + abs(deltaV2);

% Time of flight (half period of the ellipse)
TOF_sec = pi * sqrt(a_t^3 / mu_sun);  % seconds
TOF_yrs = TOF_sec / (86400 * 365.25); % years
