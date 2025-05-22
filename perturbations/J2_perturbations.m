clc
clear
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   VARIATION OF PARAMETERS (J2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % rates_of_change.m
% % Compute nodal, perigee and mean motion rates due to J2
% 
% %---- Define parameters (all initialized to zero) ----
% J2 = 0;                 % Planets second zonal harmonic
% mu = 0;                 % Gravitational parameter [km^3/s^2]
% R  = 0;                 % Earth's equatorial radius [km]
% a  = 0;                 % Semi-major axis [km]
% e  = 0;                 % Eccentricity [–]
% i  = 0*pi/180;          % Inclination [deg]
% n  = sqrt(mu/a^3);      % Mean motion [rad/s]
% 
% %---- Precompute common factor ----
% factor = (3*J2*sqrt(mu)*R^2) / (2 * a^(7/2) * (1 - e^2)^2);
% 
% %---- Compute the rates ----
% Omega_dot   = - factor * cos(i);
% omega_dot   = - factor * ( (5/2)*sin(i)^2 - 2 );
% theta_dot   = n - factor * ( 1 - (3/2)*sin(i)^2 );
% 
