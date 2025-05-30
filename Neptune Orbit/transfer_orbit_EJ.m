clc
clear
close all

% Define orbit around Earth 
a = 1008677538.4601;    % Semi-major axis (km)
e = 0.852087899;        % Eccentricity
i = 1.5183;             % Inclination (º)
RA = 5.3022;            % Right Ascension of Ascending Node (º)
w = 183.5846;           % Argument of Periapse (º)

% Define parameters
mu = 398600;                    % Gravitational Parameter
h = sqrt(mu * a * (1 - e^2));   % Specific angular momentum

% Position and velocity in perifocal frame
theta1 = 356.4033-360;                  % Earth's true anomaly (º) (heliocentric)
theta2 = 137.9554;                      % Jupiter's true anomaly (º) (heliocentric)
theta = theta1:0.01:theta2;
r = ((h^2)/mu)./(1+e*cosd(theta));

r_vec_peri = r.*[cosd(theta); sind(theta); zeros(size(theta))];               % Unperturbed position vector
v_vec_peri = (mu/h).*[-sind(theta); (e + cosd(theta)); zeros(size(theta))];   % Unperturbed vector transfer

% Position and velocity in inertial frame
A = [cosd(RA)*cosd(w) - sind(RA)*sind(w)*cosd(i), -cosd(RA)*sind(w) - sind(RA)*cosd(w)*cosd(i),  sind(RA)*sind(i);
            sind(RA)*cosd(w) + cosd(RA)*sind(w)*cosd(i), -sind(RA)*sind(w) + cosd(RA)*cosd(w)*cosd(i), -cosd(RA)*sind(i);
            sind(w)*sind(i), cosd(w)*sind(i), cosd(i)];

r_vec_inertial = A*r_vec_peri;
v_vec_inertial = A*v_vec_peri;


figure
plot3(r_vec_inertial(1,:),r_vec_inertial(2,:),r_vec_inertial(3,:));
hold on
scatter3(r_vec_inertial(1,1),r_vec_inertial(2,1),r_vec_inertial(3,1),'filled', 'MarkerFaceColor', [0,255,0]/255);
scatter3(r_vec_inertial(1,end),r_vec_inertial(2,end),r_vec_inertial(3,end),'filled', 'MarkerFaceColor', [255,0,0]/255);
hold off

% Get initial position ad velocity vectors and put it in the format that
% Ishaan's code likes
r0_sc = r_vec_inertial(:,1);
v0_sc = v_vec_inertial(:,1);