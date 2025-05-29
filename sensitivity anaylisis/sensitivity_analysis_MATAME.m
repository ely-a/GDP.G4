clc
clear
close all

% Sensitivity Analysis

% Define knwon values
R1 = 149.6e6;       % Distance to Earth from Sun (km)
R2 = 227.9e6;         % Distance to Jupiter from Sun (km)
mu_S = 1.327e11;    % Sun's gravitational parameter (km^3/s^2)
mu_E = 398600;      % Earth's gravitational parameter (km^3/s^2)
r_p = 6678;         % Distance of spacecraft from Earth (km)


% Define the transfer orbit parameters
theta1 = 0;
theta2 = 180;
e = (R2-R1)/(R1*cosd(theta1) - R2*cosd(theta2));   % Eccentricity of transfer orbit
h = sqrt(mu_S*R1*(1+e*cosd(theta1)));              % Specific angular momentum of transfer orbit (km^2/s)
a = (R1+R2)/2;

% Define velocities
v_E = sqrt(mu_S/R1);                    % Velocity of Earth (km/s)
v_r = (mu_S/h)*e*sind(theta1);          % Radial velocity of spacecraft at departure (km/s)
v_per = (mu_S/h)*(1+e*cosd(theta1));    % Perpendicular velocity of spacecraft at departure (km/s)
v_inf = 2.943;                          % Departure hyperbolic excess velocity (km/s)
v_A = sqrt((v_r^2) + (v_per^2));        % Velocity magnitude of spcecraft at departure (km/s)
v_p = sqrt((v_inf^2) + (2*mu_E/r_p));   % Burnout velocity 
v_A1 = sqrt(mu_S * ((2/R1)-(1/a)));

% Define variation parameters
dvA_drp = mu_E/(v_inf*(r_p^2));     % Change in v_A due to variations in r_p
dvA_dvp = v_p/v_inf;                % Change in v_A due to variations in v_p

% Write factors
factor1 = (2/(1-((R1*v_A^2)/(2*mu_S))));
factor2 = (2*mu_S*(cosd(theta1) - cosd(theta2)))/(((mu_S*(cosd(theta1) - cosd(theta2)) + (v_per^2)*R1*cosd(theta2))));
k_rp = mu_E/(v_inf*r_p*v_A);
k_vp = ((v_inf)+(2*mu_E/(r_p*v_inf)))/(v_A);

%substitue into expression for dR2/R2
rp_variation1 = factor1*k_rp;
vp_variation1 = factor1*k_vp;

rp_variation2 = factor2*k_rp;
vp_variation2 = factor2*k_vp;
