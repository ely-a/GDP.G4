clear; clc; close all;

%% === CONSTANTS ===
mu_N = 6.836529e6;      % Neptune gravitational parameter [km^3/s^2]
mu_T = 1.42e3;          % Triton's gravitational parameter [km^3/s^2]
R_N = 24622;            % Neptune radius [km]
R_T = 1353.4;           % Triton radius [km]
inclination_deg = 130;  % Triton's retrograde inclination [deg]
incl = deg2rad(inclination_deg);

%% === SPACECRAFT INITIAL ORBIT ===
ra_cap = 563400;                % Apoapsis [km]
rp_cap = R_N + 200;             % Periapsis [km]
a_cap = (rp_cap + ra_cap)/2;
e_cap = (ra_cap - rp_cap)/(ra_cap + rp_cap);
T_cap = 2 * pi * sqrt(a_cap^3 / mu_N);  % Period [s]

% Full true anomaly range (complete orbit)
theta_arr = linspace(0, 2*pi, 800);
r_arr = a_cap * (1 - e_cap^2) ./ (1 + e_cap * cos(theta_arr));

% Time of flight from apoapsis
E_arr = 2 * atan(sqrt((1 - e_cap)/(1 + e_cap)) .* tan(theta_arr/2));
M_arr = E_arr - e_cap .* sin(E_arr);
tof_arr = T_cap * M_arr / (2*pi);  % [s]

% Spacecraft trajectory in inclined plane
x_sc = r_arr .* cos(theta_arr);
y_sc = r_arr .* sin(theta_arr);
z_sc = r_arr .* sin(incl) .* sin(theta_arr);

%% === TRITON ORBIT: STARTING FROM ACTUAL POSITION ===
r_T0 = [-137198.419486464; 177212.813522969; 275007.819026313];  % km
r_T_mag = norm(r_T0);
v_T_mag = 4.39;                        % km/s
omega_T = v_T_mag / r_T_mag;          % rad/s

% Compute rotation axis from orbital angular momentum
z_hat = [0; 0; 1];
h_T = cross(r_T0, z_hat);             % Angular momentum direction
rot_axis = h_T / norm(h_T);           % Unit vector

% Preallocate Triton positions
x_T = zeros(1, length(tof_arr));
y_T = zeros(1, length(tof_arr));
z_T = zeros(1, length(tof_arr));

% Propagate using Rodrigues' rotation formula
for k = 1:length(tof_arr)
    theta_T = omega_T * tof_arr(k);
    r_rot = rodrigues_rotation(r_T0, rot_axis, theta_T);
    x_T(k) = r_rot(1);
    y_T(k) = r_rot(2);
    z_T(k) = r_rot(3);
end

%% === FINAL ORBIT AFTER FLYBY ===
rp_final = R_N + 2000;  % km
ra_final = 2.5e5;       % km
a_final = (rp_final + ra_final)/2;
e_final = (ra_final - rp_final)/(ra_final + rp_final);
theta_post = linspace(0, 2*pi, 500);
r_post = a_final * (1 - e_final^2) ./ (1 + e_final * cos(theta_post));
x_post = r_post .* cos(theta_post);
y_post = r_post .* sin(theta_post);
z_post = r_post .* sin(incl) .* sin(theta_post);

%% === PLOT AND ANIMATE EVERYTHING TOGETHER ===
figure('Color','w');
axis equal
hold on
grid on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Neptune Flyby and Final Orbit')
xlim([-6e5 6e5])
ylim([-6e5 6e5])
zlim([-6e5 6e5])
view(3)

% Plot Neptune
[XS, YS, ZS] = sphere(30);
surf(R_N*XS, R_N*YS, R_N*ZS, 'EdgeColor','none', ...
    'FaceAlpha', 0.3, 'FaceColor','b');

% Pre-allocate plot handles
sc_traj = plot3(nan, nan, nan, 'r', 'LineWidth', 1.4);
sc_dot = plot3(nan, nan, nan, 'ro', 'MarkerFaceColor','r');
tr_dot = plot3(nan, nan, nan, 'co', 'MarkerFaceColor','c');
final_traj = plot3(nan, nan, nan, 'g', 'LineWidth', 1.8);

% Animate spacecraft and Triton
for i = 2:length(x_sc)
    set(sc_traj, 'XData', x_sc(1:i), 'YData', y_sc(1:i), 'ZData', z_sc(1:i));
    set(sc_dot, 'XData', x_sc(i), 'YData', y_sc(i), 'ZData', z_sc(i));
    set(tr_dot, 'XData', x_T(i), 'YData', y_T(i), 'ZData', z_T(i));
    drawnow;
    pause(0.005)
end

% Continue animation with final orbit
for j = 1:length(x_post)
    set(final_traj, 'XData', x_post(1:j), 'YData', y_post(1:j), 'ZData', z_post(1:j));
    drawnow;
    pause(0.005)
end

legend({'Neptune','Spacecraft Entry Orbit','Final Orbit','Triton'}, ...
    'Location','northeastoutside')

%% === ROTATION FUNCTION ===
function r_rot = rodrigues_rotation(r, k, theta)
% Rotate vector r about axis k by angle theta (in radians)
    r = r(:); k = k(:);  % Ensure column vectors
    r_rot = r * cos(theta) + cross(k, r) * sin(theta) + k * dot(k, r) * (1 - cos(theta));
end
