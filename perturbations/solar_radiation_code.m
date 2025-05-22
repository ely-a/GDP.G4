clc 
clear
close all

% Perturbations considered
% 1) All bodies
% 2) Solar radiation
% 3) J2 Jupiter
% 4) J2 Neptune

%venus, earth, mars, jupiter, saturn, uranus, neptune

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        SOLAR RADIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trajectories = load("trajectory.mat");
r_planets = trajectories.r;
r_sc = trajectories.r_out;
n = size(r_sc,1);
r_planets = r_planets(1:n,:);
radius = [6051.8, 6378.0, 3389.5, 69911, 58232, 25362, 24622]; %in km

S0 = 63.15e6;               % surface radiated power intensity
R_S = 696340;               % radius of the sun in km
R_sc = vecnorm(r_sc,2,2);   % distnace of probe from the sun in km
S = S0*(R_S./R_sc).^2;      % solar constant in W/m^2
c = 2.998e8;                % speed of light m/s^
P = S/c;                    % solar radiation pressure N/m^2

%Cannonball model - assume the spacecraft is a sphere
R = 3;                  % radius of spacecraft in m
v = zeros(n,7);         % shadow function
C_R = 1.25;             % radiation pressure coefficient (1-2)
A_s = pi*R^2;           % surface area of sc (m^2)
m = 42500;              % mass of sc (kg)

%Work out v for different positions
for i = 1:7
    r_b = r_planets(:,((3*i)-2):(3*i)); %in km
    R_b = vecnorm( r_b , 2 , 2 );       %in km

    theta = acosd((dot(r_sc, r_b , 2 ))./(R_sc.*R_b));
    theta1 = acosd(radius(i)./R_sc);
    theta2 = acosd(radius(i)./R_b);
    
    sight = (theta1 + theta2) >= theta;  % n×1 logical
    v(sight, i) = 1;
end

v = double(any(v, 2));

F = -v.*P.*C_R.*A_s;                    %perturbing force in N
a_vec = -(F/m).*((r_sc./R_sc)*0.001);   %perturbing acceleration vector in km/s^2
a_mag = vecnorm(a_vec, 2, 2);           %perturbing acceleration magnitude in km/s^2

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
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   REFERENCE TRAJECTORY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % import from Ben's code
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   PERTURBED TRAJECTORY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Add all of the perturbations to the reference trajectory
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
