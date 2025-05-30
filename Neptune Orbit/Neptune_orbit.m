clear
clc
close all

%% SET GRAPHING FONTS AND SIZES
set(groot,'defaultLineLineWidth',4) % Set all line widths to 2 unless stated otherwise
set(groot,'defaultAxesFontSize',34) % Set all axes font size to 16 unless stated otherwise
set(groot,'defaulttextfontsize',36) % Set all text font sizes to 16 unless stated otherwise
set(groot,'defaultLineMarkerSize',12) % Set all marker sizes to 16 unless stated otherwise
set(groot,'defaultAxesXGrid','on') % Turn x grid lines on
set(groot,'defaultAxesYGrid','on') % Turn y grid lines on

%% Orbit around Neptune 

% Neptune Orbital Trajectory Design
% Author: Team 04 - Trajectory Design
% Date: May 2025

%% === CONSTANTS ===
mu_N = 6.836529e15;              % Neptune gravitational parameter [m^3/s^2]
R_N = 24622e3;                   % Neptune radius [m]
r_triton = 354800e3;             % Triton's orbital radius [m]
incl_triton = deg2rad(156.8);    % Triton's orbital inclination [rad]

%% === FLYBY TRAJECTORY (Hyperbolic) ===
rp_flyby = r_triton - 1000e3;          % 1000 km below Triton orbit
e_flyby = 1.2;                         % hyperbolic eccentricity
a_flyby = -rp_flyby / (e_flyby - 1);   % semi-major axis

theta_fly = deg2rad(linspace(-100, 100, 1000));
r_fly = a_flyby * (e_flyby^2 - 1) ./ (1 + e_flyby * cos(theta_fly));
x_fly = r_fly .* cos(theta_fly);
y_fly = r_fly .* sin(theta_fly);
z_fly = zeros(size(x_fly));

% Rotate to Triton's inclination
x_fly_rot = x_fly;
y_fly_rot = y_fly * cos(incl_triton);
z_fly_rot = y_fly * sin(incl_triton);

r_flyby_end = [x_fly_rot(end), y_fly_rot(end), z_fly_rot(end)];

%% === POLAR ORBITS ===
polar_orbits = {
    struct('rp', R_N + 2000e3, 'ra', R_N + 2000e3, 'color', 'g'),
    struct('rp', R_N + 2000e3, 'ra', R_N + 5000e3, 'color', [1 0.5 0]),
    struct('rp', R_N + 2000e3, 'ra', R_N + 15000e3, 'color', [0.5 0 0.8])
};

figure('Color', 'w'); hold on;
plot3(x_fly_rot/1e3, y_fly_rot/1e3, z_fly_rot/1e3, 'b', 'LineWidth', 1.5);

for i = 1:length(polar_orbits)
    orbit = polar_orbits{i};
    rp = orbit.rp; ra = orbit.ra; color = orbit.color;
    a = (rp + ra)/2;
    e = (ra - rp)/(ra + rp);

    theta = linspace(0, 2*pi, 600);
    r = a * (1 - e^2) ./ (1 + e * cos(theta));
    x = r .* cos(theta);
    y = zeros(size(x));
    z = r .* sin(theta);

    plot3(x/1e3, y/1e3, z/1e3, 'Color', color, 'LineWidth', 1.5);

    % Transition arc (elliptical)
    r_end = [x(1), y(1), z(1)];
    r1 = norm(r_flyby_end);
    r2 = norm(r_end);
    rp_trans = min(r1, r2);
    ra_trans = max(r1, r2);
    a_trans = (rp_trans + ra_trans) / 2;
    e_trans = (ra_trans - rp_trans) / (ra_trans + rp_trans);

    theta_trans = linspace(0, pi, 300);
    r_trans = a_trans * (1 - e_trans^2) ./ (1 + e_trans * cos(theta_trans));

    v1 = r_flyby_end / norm(r_flyby_end);
    v2 = r_end / norm(r_end);
    normal = cross(v1, v2); normal = normal / norm(normal);
    basis_x = v1;
    basis_y = cross(normal, basis_x);

    arc_pts = zeros(3, length(theta_trans));
    for j = 1:length(theta_trans)
        arc_pts(:, j) = r_trans(j) * (cos(theta_trans(j)) * basis_x' + sin(theta_trans(j)) * basis_y');
    end
    plot3(arc_pts(1, :) / 1e3, arc_pts(2, :) / 1e3, arc_pts(3, :) / 1e3, '--', 'Color', color);
end

scatter3(0, 0, 0, 200, 'c', 'filled');
% Plot Triton orbit for reference
theta_triton = linspace(0, 2*pi, 500);
x_triton = r_triton * cos(-theta_triton);
y_triton = r_triton * sin(-theta_triton) * cos(incl_triton);
z_triton = r_triton * sin(-theta_triton) * sin(incl_triton);
plot3(x_triton/1e3, y_triton/1e3, z_triton/1e3, 'k--', 'LineWidth', 1.5);

legend('Flyby Trajectory','Polar Orbit 1','Transfer 1','Polar Orbit 2','Transfer 2','Polar Orbit 3','Transfer 3','Neptune','Triton Orbit');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Triton Flyby');
grid on; axis equal; view(30, 25);



%% MATLAB Script: Neptune Arrival → Triton Orbit → Return to Neptune Orbit
%% === CONSTANTS ===
mu_N = 6.836529e15;              % Neptune gravitational parameter [m^3/s^2]
mu_triton = 1.083e12;            % Triton gravitational parameter [m^3/s^2]
R_N = 24622e3;                   % Neptune radius [m]
r_triton = 354800e3;             % Triton's orbital radius [m]
incl_triton = deg2rad(156.8);    % Triton's orbital inclination [rad]

%% === PHASE 0: Parabolic Arrival Trajectory ===
rp_arrival = R_N + 2000e3; %2000e3
e_arrival = 1;  % parabolic
theta_arrival = linspace(-pi/3, pi/3, 500);  % tighter range around periapsis
r_arrival = rp_arrival * (1 + e_arrival) ./ (1 + e_arrival * cos(theta_arrival));
x_arrival = r_arrival .* cos(theta_arrival);
y_arrival = zeros(size(x_arrival));
z_arrival = r_arrival .* sin(theta_arrival);

%% === PHASE 1: Neptune Orbit Capture (Highly Elliptical) ===
rp_cap = rp_arrival;
ra_cap = r_triton + 2000e3; %can vary but lets say closest distance to triton bring 2000 km - this turns into vector x,y,z
a_cap = (rp_cap + ra_cap)/2;
e_cap = (ra_cap - rp_cap)/(ra_cap + rp_cap);
theta_cap = linspace(0, 2*pi, 600);
r_cap = a_cap * (1 - e_cap^2) ./ (1 + e_cap * cos(theta_cap));
x_cap = r_cap .* cos(theta_cap);
y_cap = zeros(size(x_cap));
z_cap = r_cap .* sin(theta_cap);

%% === PHASE 2: Transfer to Triton ===
% ra_trans = ra_cap;
% 
% rp_trans = -r_triton - 2000e3;
% 
% %rp_trans = ra_cap;
% %ra_trans = r_triton;
% 
% a_trans = (rp_trans + ra_trans)/2;
% e_trans = (ra_trans - rp_trans)/(ra_trans + rp_trans);
% theta_trans = linspace(0, pi, 300);
% r_trans = a_trans * (1 - e_trans^2) ./ (1 + e_trans * cos(theta_trans));
% x_trans = r_trans .* cos(theta_trans);
% y_trans = zeros(size(x_trans));
% z_trans = r_trans .* sin(theta_trans);

%% === PHASE 3: Orbit Around Triton (One Revolution) ===
% mu_T = 1.083e12;  % Triton gravity [m^3/s^2] 1.43405?? 
% rp_triton_orb = 1000e3;
% ra_triton_orb = 5000e3;
% a_triton_orb = (rp_triton_orb + ra_triton_orb)/2;
% e_triton_orb = (ra_triton_orb - rp_triton_orb)/(ra_triton_orb + rp_triton_orb);
% theta_triton_orb = linspace(0, 2*pi, 500);
% r_triton_orb = a_triton_orb * (1 - e_triton_orb^2) ./ (1 + e_triton_orb * cos(theta_triton_orb));
% x_triton_orb = r_triton_orb .* cos(theta_triton_orb);
% y_triton_orb = r_triton_orb .* sin(theta_triton_orb) * cos(incl_triton);
% z_triton_orb = r_triton_orb .* sin(theta_triton_orb) * sin(incl_triton);
% x_triton_orb = x_triton_orb + r_triton;

%% === PHASE 4: Return Transfer from Triton to Final Polar Orbit ===
% rp_return = r_triton;
% ra_return = R_N + 6000e3;
% a_return = (rp_return + ra_return)/2;
% e_return = (ra_return - rp_return)/(ra_return + rp_return);
% theta_return = linspace(0, pi, 300);
% r_return = a_return * (1 - e_return^2) ./ (1 + e_return * cos(theta_return));
% x_return = r_return .* cos(theta_return);
% y_return = zeros(size(x_return));
% z_return = r_return .* sin(theta_return);

%% === PHASE 5: Final Polar Elliptical Orbit Around Neptune ===
rp_final = R_N + 2000e3; %was 1000e3
ra_final = R_N + 200000e3; %was was 50,000km 
a_final = (rp_final + ra_final)/2;
e_final = (ra_final - rp_final)/(ra_final + rp_final);
theta_final = linspace(0, 2*pi, 600);
r_final = a_final * (1 - e_final^2) ./ (1 + e_final * cos(theta_final));
x_final = r_final .* cos(theta_final);
y_final = zeros(size(x_final));
z_final = r_final .* sin(theta_final);  % polar orbit (in X-Z plane)

%% === TRITON ORBIT ===
theta_triton = linspace(0, 2*pi, 500);
x_triton = r_triton * cos(-theta_triton);
y_triton = r_triton * sin(-theta_triton) * cos(incl_triton);
z_triton = r_triton * sin(-theta_triton) * sin(incl_triton);

%% === DELTA-V CALCULATIONS ===
% Gravitational parameter of Neptune
mu = mu_N;

% 1. From parabolic arrival to elliptical capture (at periapsis)
v_parabolic = sqrt(2 * mu_N / rp_arrival);
h_cap = sqrt(rp_cap * mu_N * (1 + e_cap));
vp_capture = h_cap / rp_cap;
%v_capture = sqrt(mu_N * (2/rp_cap - 1/a_cap));
dv_capture = abs(v_parabolic - vp_capture);

% Flyby at Triton 
v_triton = sqrt(mu_N/r_triton);
v_spc = h_cap / ra_cap;
vinf = v_triton - v_spc;  % Assumed flyby vinf
turn_angle = 2.3652;

v_after = sqrt(vinf^2 + v_triton^2);
h_flyby = v_after * r_triton;
v_p_fromtriton = h_flyby / rp_arrival;

vabeforeflyby = h_cap / (ra_cap);

% ΔV equivalent effect of gravity assist:
delta_v_flyby = 2 * vinf * sin(turn_angle / 2) ;%check this 
v_p_fromtriton = abs(vabeforeflyby - delta_v_flyby) ;%this is meant to be va before flyby not vp_capture, also - or +????

% 2. From capture orbit to Triton transfer (at apoapsis of capture)
%v1 = sqrt(mu_N * (2/ra_cap - 1/a_cap));
%v2 = sqrt(mu_N * (2/ra_cap - 1/a_trans));
%dv_transfer_to_triton = abs(v1 - v2);

% 3. From Triton transfer to return to Neptune (assumed symmetrical)
%v3 = sqrt(mu * (2/rp_return - 1/a_return));
%v4 = sqrt(mu * (2/rp_return - 1/a_final));
%dv_return_from_triton = abs(v3 - v4);

% 5. From capture orbit to polar orbit 
h_final = sqrt(rp_arrival * mu_N * (1 + e_final));
v4 = h_final/rp_arrival;
dv_return = abs(v_p_fromtriton - v4);


%% === PLOTTING ===
figure('Color','w'); hold on;
plot3(x_arrival/1e3, y_arrival/1e3, z_arrival/1e3, 'c--', 'LineWidth', 1.2);
plot3(x_cap/1e3, y_cap/1e3, z_cap/1e3, 'b', 'LineWidth', 1.5);
plot3(x_triton/1e3, y_triton/1e3, z_triton/1e3, 'k--', 'LineWidth', 1.2);
plot3(x_final/1e3, y_final/1e3, z_final/1e3, 'r--', 'LineWidth', 1.5);
scatter3(0, 0, 0, 200, 'c', 'filled');

legend('Parabolic Arrival','Neptune Capture Orbit','Triton Orbit','Final Polar Orbit','Neptune');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Neptune Capture → Triton Orbit → Return to Polar Neptune Orbit');
grid on; axis equal; view(30, 25);




%% Time Calculations

% Transfer time from Neptune to Triton 

P_phase1 = 2 * 3.14 * a_cap ^ (3/2) / sqrt(mu_N);
t_phase1 = P_phase1/2; 

% ADD THIS TO THE START OF LAUNCH DATE TO FIND THE VECTORS FOR TRITON on
% jpl horizons V

%%


r_triton = [-3.5481e8; 1.549e7; 1.243e6];     % meters
v_triton = [-52.3; -4379.2; -2.78];          % m/s


%% obtaining position and velocity vector of spacecraft at start of trajectory 

h = cross(r_triton, v_triton);        
e_vec = cross(v_triton, h)/mu_N - r_triton/norm(r_triton);  
theta_T = acos(dot(e_vec, r_triton) / (norm(e_vec) * norm(r_triton)));

% If dot(r, v) < 0, then angle is > 180°
if dot(r_triton, v_triton) < 0
    theta_T = 2*pi - theta_T;
end



i_triton = deg2rad(157);  
%theta_sc = deg2rad(180 - theta_T); 

theta_sc = mod(theta_T + pi, 2*pi);  % spacecraft at periapsis, opposite Triton

rp_cap = R_N + 2000e3;  % periapsis
ra_cap = norm(r_triton) + 2000e3;
a_cap = (rp_cap + ra_cap)/2;
e_cap = (ra_cap - rp_cap)/(ra_cap + rp_cap);

R_N 
norm(r_triton)

r_mag = a_cap * (1 - e_cap^2) / (1 + e_cap * cos(theta_sc))  % radius at true anomaly



[i_triton, Omega, omega, theta] = orbital_elements(r_triton, v_triton, mu_N);


r_sc = r_mag * [cos(theta_sc); sin(theta_sc) * cos(i_triton); sin(theta_sc) * sin(i_triton)];

relative_position = r_sc - r_triton;
distance = norm(relative_position) 

angle_between = acosd(dot(r_triton, r_sc) / (norm(r_triton) * norm(r_sc)));
fprintf("3D angle between vectors: %.2f deg\n", angle_between);


% Constants
R_N = 24622e3;              % Neptune radius in meters
altitude_perigee = 2000e3;  % 2000 km above surface
rp = R_N + altitude_perigee;

% Triton position vector at some epoch (example vector)
r_triton = [-3.5481e8; 1.549e7; 1.243e6];  % meters 
ra = norm(r_triton);

% Orbital parameters
a = (rp + ra)/2;
e = (ra - rp) / (ra + rp);  % or e = (ra - rp)/(2*a)

% Perigee position vector (in opposite direction of Triton)
unit_direction = -r_triton / norm(r_triton);  % unit vector pointing toward perigee
r_perigee_vec = unit_direction * rp;

% Output
disp('Spacecraft position vector at perigee: ');
disp(r_perigee_vec);


relative_position = r_perigee_vec - r_triton;
distance = norm(relative_position) 

angle_between = acosd(dot(r_triton, r_perigee_vec) / (norm(r_triton) * norm(r_perigee_vec)));
fprintf("3D angle between vectors: %.2f deg\n", angle_between);




R_N = 24622e3;
r_triton_vec = [-3.5481e8; 1.549e7; 1.243e6];  % meters
v_mag = 22e3;  % 22 km/s

[r_p, v_p] = get_initial_state(r_triton_vec, R_N, v_mag)

norm(v_p)


% Step 1: Position direction (Neptune to perigee)
r_dir = -r_triton_vec / norm(r_triton_vec);  % unit vector

% Step 2: Arbitrary vector not parallel to r_dir
ref = [0; 0; 1];
if abs(dot(r_dir, ref)) > 0.99
    ref = [0; 1; 0];
end

% Step 3: Compute tangential direction (perpendicular to r_dir)
v_tangent = cross(ref, r_dir);
v_tangent = v_tangent / norm(v_tangent);  % unit

% Step 4: Compute binormal (normal to orbital plane)
v_binormal = cross(r_dir, v_tangent);  % could be useful if needed

% Step 5: Assign velocity
v_perigee_vec = v_tangent * v_mag   % this just gave positive values of other one 



[r, v] = get_state_in_orbital_plane(r_triton_vec, v_triton, R_N, v_mag)

norm(v)


%% Calculating Delta Vs

r_triton = [-3.5481e5; 1.549e4; 1.243e3];     
v_triton = [-52.3e-3; -4379.2e-3; -2.78e-3];    

[total_dv, a_return, e_return, h_return] = DV_calc(r_triton, v_triton);


% % === CONSTANTS ===
% mu_N = 6.8365e6;           % Neptune gravitational parameter [km^3/s^2]
% R_N = 24760;               % Neptune radius [km]
% R_Triton = 1353.4;         % Triton radius [km]
% r_Triton_orbit = 354800;   % Circular orbital radius of Triton [km]
% incl_elliptical = 30;      % Elliptical orbit inclination [deg]
% incl_polar = 90;           % Final polar orbit inclination [deg]
% incl_triton = 23;          % Triton orbital inclination [deg]
% 
% % === RESOLUTION ===
% n_elliptical = 300;
% n_flyby = 50;
% n_polar = 300;
% n_total = n_elliptical + n_flyby + n_polar;
% 
% % === ELLIPTICAL ORBIT ===
% rp_ell = R_N + 2000;
% ra_ell = r_Triton_orbit + 10000;
% a_ell = (rp_ell + ra_ell)/2;
% e_ell = (ra_ell - rp_ell)/(ra_ell + rp_ell);
% 
% % === POLAR ORBIT ===
% rp_polar = R_N + 2000;
% ra_polar = R_N + 120000;
% a_polar = (rp_polar + ra_polar)/2;
% e_polar = (ra_polar - rp_polar)/(ra_polar + rp_polar);
% 
% % === Vis-Viva Spacing ===
% visviva_spacing = @(a,e,mu,N) ...
%     interp1(cumsum(1 ./ sqrt(mu * (2 ./ (a * (1 - e^2) ./ (1 + e * cos(linspace(0,2*pi,1000))) - 1 ./ a))), 'linear', 'extrap'), ...
%             linspace(0,1,N), ...
%             linspace(0,2*pi,1000));
% 
% theta_ell = visviva_spacing(a_ell, e_ell, mu_N, n_elliptical);
% theta_polar = visviva_spacing(a_polar, e_polar, mu_N, n_polar);
% theta_Triton = linspace(0, 2*pi, n_total);
% 
% % === Orbit Coordinates ===
% orbit_coords = @(a,e,theta,inc) [ ...
%     a*(1 - e^2)./(1 + e*cos(theta)) .* cos(theta); ...
%     a*(1 - e^2)./(1 + e*cos(theta)) .* sin(theta); ...
%     a*(1 - e^2)./(1 + e*cos(theta)) .* sin(theta)*sind(inc)];
% 
% orbit_ell = orbit_coords(a_ell, e_ell, theta_ell, incl_elliptical);
% orbit_polar = orbit_coords(a_polar, e_polar, theta_polar, incl_polar);
% triton_path = orbit_coords(r_Triton_orbit, 0, theta_Triton, incl_triton);
% 
% % === FLYBY ARC ===
% flyby_start = orbit_ell(:, end);
% flyby_end = orbit_polar(:, 1);
% t_fly = linspace(0, 1, n_flyby);
% flyby_mid = (flyby_start + flyby_end)/2 + [0; 0; 8000];  % bump up
% flyby_arc = (1 - t_fly).^2 .* flyby_start + ...
%             2 * (1 - t_fly) .* t_fly .* flyby_mid + ...
%             t_fly.^2 .* flyby_end;
% 
% % === FULL TRAJECTORY ===
% trajectory = [orbit_ell, flyby_arc, orbit_polar];
% 
% % === SETUP PLOT ===
% figure;
% hold on; grid on; axis equal;
% xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
% title('3D Neptune Flyby with Polar Capture');
% 
% % === Neptune ===
% [sx, sy, sz] = sphere(40);
% surf(R_N*sx, R_N*sy, R_N*sz, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% plot3(0, 0, 0, 'bo', 'MarkerFaceColor', 'blue', 'DisplayName', 'Neptune');
% 
% % === Orbits ===
% plot3(orbit_ell(1,:), orbit_ell(2,:), orbit_ell(3,:), 'b--', 'DisplayName', 'Elliptical Orbit');
% plot3(flyby_arc(1,:), flyby_arc(2,:), flyby_arc(3,:), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Flyby Arc');
% plot3(orbit_polar(1,:), orbit_polar(2,:), orbit_polar(3,:), 'g--', 'DisplayName', 'Polar Orbit');
% plot3(triton_path(1,:), triton_path(2,:), triton_path(3,:), 'k--', 'LineWidth', 1.2, 'DisplayName', 'Triton Orbit');
% 
% % === Spacecraft marker and trail ===
% sc = plot3(trajectory(1,1), trajectory(2,1), trajectory(3,1), ...
%     'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Spacecraft');
% trail = plot3(NaN, NaN, NaN, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Trajectory Trail');
% 
% % === Triton Sphere (scaled up) ===
% [u, v] = meshgrid(linspace(0, 2*pi, 30), linspace(0, pi, 20));
% xT = R_Triton*2.5 * cos(u) .* sin(v);
% yT = R_Triton*2.5 * sin(u) .* sin(v);
% zT = R_Triton*2.5 * cos(v);
% Triton = surf(xT + triton_path(1,1), ...
%               yT + triton_path(2,1), ...
%               zT + triton_path(3,1), ...
%               'FaceColor', 'orange', 'EdgeColor', 'black', ...
%               'LineWidth', 0.5, 'FaceAlpha', 1, 'DisplayName', 'Triton');
% 
% legend;
% view(45, 25);
% xlim([-4e5, 2e5]); ylim([-3e5, 3e5]); zlim([-2.5e5, 2.5e5]);
% 
% % === ANIMATION LOOP ===
% for i = 1:n_total
%     % Spacecraft update
%     sc.XData = trajectory(1,i);
%     sc.YData = trajectory(2,i);
%     sc.ZData = trajectory(3,i);
% 
%     trail.XData = trajectory(1,1:i);
%     trail.YData = trajectory(2,1:i);
%     trail.ZData = trajectory(3,1:i);
% 
%     % Triton update
%     Triton.XData = xT + triton_path(1,i);
%     Triton.YData = yT + triton_path(2,i);
%     Triton.ZData = zT + triton_path(3,i);
% 
%     drawnow;
% end





%% Functions

function [total_dv, a_return, e_return, h_return] = DV_calc(r_Triton_vec, v_Triton_vec)
% r_Triton_vec: 3x1 vector [km] of Triton position relative to Neptune
% v_Triton_vec: 3x1 vector [km/s] of Triton velocity relative to Neptune

% === CONSTANTS ===
    mu_N = 6.8365e6;           % Neptune gravitational parameter [km^3/s^2]
    R_N = 24760;               % Neptune radius [km]

% === Initial parabolic arrival === 
    rp_arrival = R_N + 2000; % Adjust for altitude above Neptune - km

% === PHASE 1: CAPTURE DELTA-V ===
    % Post-capture orbit
    ra_cap = norm(r_Triton_vec) + 2000; % Adjust for altitude above Triton
    rp_cap = rp_arrival;                  
    a_cap = (rp_cap + ra_cap) / 2;
    e_cap = (ra_cap - rp_cap) / (ra_cap + rp_cap);

    v_parabolic = sqrt(2 * mu_N / rp_cap); % Formula from vis-viva equation, setting 1/a = 0 
    h_cap = sqrt(mu_N * rp_cap * (1 + e_cap));
    vp_cap = h_cap / rp_cap;
    dv_capture = abs(v_parabolic - vp_cap);

% === PHASE 2: FLYBY AT TRITON ===
    r_hat = r_Triton_vec / norm(r_Triton_vec);
    v_sc_mag = sqrt(mu_N * (2 / norm(r_Triton_vec) - 1 / a_cap));
    v_sc_vec = v_sc_mag * r_hat;       

    % v_hat = cross(cross(r_Triton_vec, v_Triton_vec), r_Triton_vec);
    % v_hat = v_hat / norm(v_hat);
    % v_sc_vec = v_sc_mag * v_hat; 


    % Relative hyperbolic excess velocity
    vinf_vec = v_sc_vec - v_Triton_vec;
    vinf_mag = norm(vinf_vec);
    vinf_hat = vinf_vec / vinf_mag;

    % Calculating the turn angle flyby 
    rp_final = R_N + 2000;                 % 2000 km above Neptune
    rp_flyby = norm(r_Triton_vec) + 1000;  % 1000 km above Triton
    
    e_hyper = 1 + (rp_flyby * vinf_mag^2) / mu_N;
    ra_finaltrans = norm(r_Triton_vec);
    e_finaltrans = (ra_finaltrans - rp_final)/(ra_finaltrans + rp_final)
    delta_fin = 2 * asin(1/e_hyper);
    %deltadeg_fin = rad2deg(delta_fin);
    %dv_flyby = 2 * vinf_mag * sin(delta_fin / 2);

    % Rotation axis - normal to the plane of flyby (r × vinf)
    n_hat = cross(r_hat, vinf_hat);
    n_hat = n_hat / norm(n_hat);

    % Rotate vinf by turn angle δ - using Rodrigues' formula
    v_rotated = vinf_mag * ( ...
        vinf_hat * cos(delta_fin) + ...
        cross(n_hat, vinf_hat) * sin(delta_fin) ...
    );
    
    % Add back Triton's velocity (Neptune frame)
    v_after_vec = v_Triton_vec + v_rotated;
    
    % Return orbit (Neptune frame)
    r_flyby = norm(r_Triton_vec);
    v_after_mag = norm(v_after_vec);
    epsilon_return = v_after_mag^2 / 2 - mu_N / r_flyby;

    % Catch if return is hyperbolic
    if epsilon_return > 0
       error('Return orbit is hyperbolic (unbound)');
    end

    a_return = -mu_N / (2 * epsilon_return);
    rp_return = rp_arrival;

    % Final speed at periapsis of return orbit (based on energy)
    vp_return = sqrt(mu_N * (2 / rp_return - 1 / a_return));

    % Angular momentum and eccentricity calculation
     v_r = dot(v_after_vec, r_hat);
     v_perp_vec = v_after_vec - v_r * r_hat;
     v_perp = norm(v_perp_vec);
     
     h_return = v_perp * r_flyby;
     a_return = -mu_N / (2 * epsilon_return);
     e_return = sqrt(1 - h_return^2 / (mu_N * a_return))
    
    % Check turn angle 
    vinf_in = v_sc_vec - v_Triton_vec;
    vinf_out = v_after_vec - v_Triton_vec;
    
    delta_true = acosd(dot(vinf_in, vinf_out) / (norm(vinf_in) * norm(vinf_out)))

 % === PHASE 3: FINAL POLAR ORBIT AROUND NEPTUNE ===
    ra_final = R_N + 200000;                              % Larger than the apoasis of closest moon ~ 200,000 km (apo: 120,000)
    %a_final = (rp_final + ra_final) / 2;
    e_final = (ra_final - rp_final) / (ra_final + rp_final)
    h_final = sqrt(mu_N * rp_final * (1 + e_final));
    vp_final = h_final / rp_final;

    dv_return = abs(vp_return - vp_final);

 % === TOTAL DELTA-V ===
    total_dv = dv_capture + dv_return;

 % === OUTPUT ===
    fprintf('\n--- DELTA-V ---\n');
    fprintf('Capture ΔV:       %.2f km/s\n', dv_capture);
    fprintf('Return ΔV:        %.2f km/s\n', dv_return);
    fprintf('TOTAL ΔV:         %.2f km/s\n', total_dv);
    fprintf('-----------------------\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_perigee_vec, v_perigee_vec] = get_initial_state(r_triton_vec, R_N, v_mag)
    % Inputs:
    % - r_triton_vec: Triton's position vector from Neptune [3x1]
    % - R_N: Neptune radius [m]
    % - v_mag: magnitude of spacecraft velocity at perigee [m/s]

    altitude_perigee = 2000e3;  % 2000 km
    rp = R_N + altitude_perigee;

    % Perigee position vector: opposite direction of Triton
    r_dir = -r_triton_vec / norm(r_triton_vec);
    r_perigee_vec = r_dir * rp;

    % Construct velocity vector perpendicular to position vector
    % Choose an arbitrary reference vector not parallel to r_dir
    ref = [0; 0; 1];
    if abs(dot(ref, r_dir)) > 0.99
        ref = [0; 1; 0];  % fallback if near parallel
    end

    % Get a vector perpendicular to r_dir
    v_dir = cross(r_dir, ref);
    v_dir = v_dir / norm(v_dir);  % unit vector

    % Scale to desired speed
    v_perigee_vec = v_dir * v_mag;

    % Output
    disp('Position vector at perigee:');
    disp(r_perigee_vec);
    disp('Velocity vector at perigee:');
    disp(v_perigee_vec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [i, Omega, omega, theta] = orbital_elements(r, v, mu)
    % Angular momentum vector
    h = cross(r, v);
    h_mag = norm(h);
    
    % Inclination
    i = acos(h(3) / h_mag);

    % Node vector (points toward ascending node)
    K = [0; 0; 1];
    n = cross(K, h);
    n_mag = norm(n);
    
    % Eccentricity vector
    e_vec = (cross(v, h) / mu) - r / norm(r);
    e = norm(e_vec);
    
    % Longitude of ascending node (Ω)
    if n_mag ~= 0
        Omega = acos(n(1) / n_mag);
        if n(2) < 0
            Omega = 2*pi - Omega;
        end
    else
        Omega = 0;  % Circular equatorial orbit
    end

    % Argument of periapsis (ω)
    if e > 1e-8  % non-circular
        omega = acos(dot(n, e_vec) / (n_mag * e));
        if e_vec(3) < 0
            omega = 2*pi - omega;
        end
    else
        omega = 0;  % circular orbit
    end

    % True anomaly (θ)
    theta = acos(dot(e_vec, r) / (e * norm(r)));
    if dot(r, v) < 0
        theta = 2*pi - theta;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_p_vec, v_p_vec] = get_state_in_orbital_plane(r_triton_vec, v_triton_vec, R_N, v_mag)

    % Step 1: Define perigee position vector (opposite to Triton)
    altitude_perigee = 2000e3;
    rp = R_N + altitude_perigee;
    r_dir = -r_triton_vec / norm(r_triton_vec);
    r_p_vec = rp * r_dir;

    % Step 2: Compute orbital plane normal (angular momentum direction)
    h_vec = cross(r_triton_vec, v_triton_vec);
    h_unit = h_vec / norm(h_vec);  % unit vector normal to orbital plane

    % Step 3: Tangential direction in the orbital plane
    v_dir = cross(h_unit, r_dir);  % orthogonal to r_dir and in the plane
    v_dir = v_dir / norm(v_dir);

    % Step 4: Assign velocity vector
    v_p_vec = v_mag * v_dir;

    % Output
    disp('Position vector at perigee:');
    disp(r_p_vec);
    disp('Velocity vector at perigee:');
    disp(v_p_vec);
end

%% 

save('vectors','r','v')