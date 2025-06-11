clear
clc
close all

%% SET GRAPHING FONTS AND SIZES
set(groot,'defaultLineLineWidth',2) % Set all line widths to 2 unless stated otherwise
set(groot,'defaultAxesFontSize',16) % Set all axes font size to 16 unless stated otherwise
set(groot,'defaulttextfontsize',16) % Set all text font sizes to 16 unless stated otherwise
set(groot,'defaultLineMarkerSize',12) % Set all marker sizes to 16 unless stated otherwise
set(groot,'defaultAxesXGrid','on') % Turn x grid lines on
set(groot,'defaultAxesYGrid','on') % Turn y grid lines on

%% Neptune Constants

mu_N = 6.8351e6; % km^3/s^2
mass_neptune = 102.409e24; % kg
r_N = 24622;             % km (mean radius)
r_N_equ = 24764;         % km (equatorial radius)

%% Load drag data

A = pi * 3^2;          % m^2
m = (1740 - 380) * 1e6;              % kg
C_D = 1.22;
beta_nom = m / (C_D * A);

% Build interpolant functions for nominal run
load("gram_profiles.mat"); % rho_mean, ew_mean, ns_mean, alt_steps
atmos_nom = @(alt_km) interp1(alt_steps, rho_mean * 1e9, alt_km, 'linear', 0); % kg/km^3
zonal_nom = @(alt_km) interp1(alt_steps, ew_mean/1000, alt_km, 'linear', 0);   % km/s
merid_nom = @(alt_km) interp1(alt_steps, ns_mean/1000, alt_km, 'linear', 0);   % km/s

%% Load moon data

if exist('moon_orbits.mat', 'file') ~= 2 
    load_moon_data()
end 
load('moon_orbits.mat','orbits_mat','moon_names', 'max_len', 'times_all', 'velocities_mat');

%% Unperturbed Trajectory

if exist('ScienceOrbit.mat', 'file') == 2
    load("ScienceOrbit.mat")
else
    error('no intitial conditions')
end

% [r0, v0] = rv_from_oe(3e5, 0.9, 90, 90, 250, 0, mu_neptune);
[r0, v0] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, 180, mu_N);
[a, e, h, i, Omega, omega, theta] = oe_from_rv(r0, v0, mu_N);
P = 2 * pi * sqrt(a^3 / mu_N);  % seconds
N_orbits = 1;
T = P * N_orbits;

delta_theta = 1;
thetas = mod(theta:delta_theta:theta+360 * T/P, 360);

r_unperturbed = zeros(3, length(thetas));
v_unperturbed = zeros(3, length(thetas));

r_unperturbed(:,1) = r0;
v_unperturbed(:,1) = v0;

for j = 2:length(thetas)
    [r_j, v_j] = rv_from_oe(a, e, i, Omega, omega, thetas(j), mu_N);
    r_unperturbed(:,j) = r_j;
    v_unperturbed(:,j) = v_j;
end

%% N-body perturbations

jd0 = 2.466559716824276e+06;

idx_triton = find(strcmpi("Triton", moon_names));  % case-insensitive
idx_sun = find(strcmpi("Sun", moon_names));  % case-insensitive

%% Perturbed trajectory

% Initial state vector (km and km/s)
Y0 = [r_unperturbed(:,1); v_unperturbed(:,1)];

% Orbital period (estimate from semi-major axis)
tspan = [0, T];         % time steps

% Integrate perturbed motion
opts = odeset('RelTol',1e-9, 'AbsTol',1e-9);
perturbs = ["J2", "drag"];

orbits_sun_triton = orbits_mat(:,:,[idx_sun, idx_triton]);
times_sun_triton = {times_all{idx_sun}, times_all{idx_triton}};
mu_sun = 1.327e11;
mu_triton = 1428.495;
mu_sun_triton = [mu_sun, mu_triton];

[t_out, Y_out] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, perturbs, jd0, ...
                                             orbits_sun_triton, times_sun_triton, mu_sun_triton, ...
                                             atmos_nom, zonal_nom, merid_nom, beta_nom) ...
                                             , tspan, Y0, opts);

r_perturbed = Y_out(:,1:3)';
v_perturbed = Y_out(:,4:6)';

%% Plot trajectories

figure;
plot3(r_unperturbed(1,:), r_unperturbed(2,:), r_unperturbed(3,:), 'b-');  % unperturbed
hold on;
plot3(r_perturbed(1,:), r_perturbed(2,:), r_perturbed(3,:), 'r--', 'LineWidth', 1.5);  % perturbed

% Mark start and end points
plot3(r_unperturbed(1,1), r_unperturbed(2,1), r_unperturbed(3,1), 'bo', 'MarkerFaceColor', 'b');
text(r_unperturbed(1,1), r_unperturbed(2,1), r_unperturbed(3,1), '  Start', 'Color', 'b');

plot3(r_unperturbed(1,end), r_unperturbed(2,end), r_unperturbed(3,end), 'bs');
text(r_unperturbed(1,end), r_unperturbed(2,end), r_unperturbed(3,end), '  End (Unperturbed)', 'Color', 'b');

plot3(r_perturbed(1,end), r_perturbed(2,end), r_perturbed(3,end), 'rs');
text(r_perturbed(1,end), r_perturbed(2,end), r_perturbed(3,end), '  End (Perturbed)', 'Color', 'r');

legend('Unperturbed', 'Perturbed (J2)', 'Location', 'best');
grid on;
axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('3D Trajectory with J2 Perturbation');

%% Orbital element variation
% Preallocate orbital elements over time
N = length(t_out);
a_vals      = zeros(1, N);
e_vals      = zeros(1, N);
i_vals      = zeros(1, N);
Omega_vals  = zeros(1, N);
omega_vals  = zeros(1, N);
theta_vals  = zeros(1, N);

for k = 1:N
    r_k = r_perturbed(:,k);
    v_k = v_perturbed(:,k);
    [a_k, e_k, h_k, i_k, Omega_k, omega_k, theta_k] = oe_from_rv(r_k, v_k, mu_N);
    a_vals(k) = a_k;
    e_vals(k) = e_k;
    i_vals(k) = i_k;
    Omega_vals(k) = Omega_k;
    omega_vals(k) = omega_k;
    theta_vals(k) = theta_k;
end

% Convert time to days
t_days = t_out / 86400;

% Plot
figure;
subplot(3,2,1); plot(t_days, a_vals); hold on; plot([0, t_days(end)], [a, a]);      ylabel('a [km]'); title('Semi-major Axis');
subplot(3,2,2); plot(t_days, e_vals);       ylabel('e'); title('Eccentricity');
subplot(3,2,3); plot(t_days, (i_vals));      ylabel('i [deg]');
subplot(3,2,4); plot(t_days, (Omega_vals));  ylabel('\Omega [deg]');
subplot(3,2,5); plot(t_days, (omega_vals));  ylabel('\omega [deg]');
subplot(3,2,6); plot(t_days, (theta_vals));  ylabel('\theta [deg]'); xlabel('Time [days]');
sgtitle('Orbital Elements Over Time (Perturbed)');

% %% J2 Omega and omega rate of change comparison
% % === Analytical J2 Drift Rates ===
% 
% J2 = 3.411e-3;              % Neptune's J2 coefficient
% mu = mu_neptune;           % km^3/s^2
% R = radius_neptune;        % km
% 
% % Use initial orbital elements (from earlier in script)
% n = sqrt(mu / a^3);        % mean motion [rad/s]
% 
% % RAAN drift (rad/s)
% Omega_dot_rad = -1.5 * J2 * n * (R/a)^2 * cosd(i) / (1 - e^2)^2;
% 
% % Argument of perigee drift (rad/s)
% omega_dot_rad = 0.75 * J2 * n * (R/a)^2 * (5 * cosd(i)^2 - 1) / (1 - e^2)^2;
% 
% % Convert to deg/day
% Omega_dot_deg_day = rad2deg(Omega_dot_rad) * 86400;
% omega_dot_deg_day = rad2deg(omega_dot_rad) * 86400;
% 
% fprintf('Analytical RAAN drift:       %.6f deg/day\n', Omega_dot_deg_day);
% fprintf('Analytical perigee drift:    %.6f deg/day\n', omega_dot_deg_day);
% 
% RAAN_drift_sim = (Omega_vals(end) - Omega_vals(1)) / (t_days(end) - t_days(1));
% omega_drift_sim = (omega_vals(end) - omega_vals(1)) / (t_days(end) - t_days(1));
% 
% fprintf('Simulated RAAN drift:        %.6f deg/day\n', RAAN_drift_sim);
% fprintf('Simulated perigee drift:     %.6f deg/day\n', omega_drift_sim);

%% J2 Perturbations

function a_J2 = J2_acc(r_vec, mu, R)
    x = r_vec(1);
    y = r_vec(2);
    z = r_vec(3);
    r = norm(r_vec);

    J2 = 3.411e-3;
    J2_const = (1.5 * J2 * mu * R^2)/(r^5);

    a_J2 = J2_const * [
        x * (5 * (z/r)^2 - 1);
        y * (5 * (z/r)^2 - 1);
        z * (5 * (z/r)^2 - 3)
    ];
end

%% N-body Perturbations

function a_N_body = N_body_acc(jd, orbits, times, mus, r_sc)
    persistent pos_sun; persistent pos_triton;

    if isempty(pos_sun)

        % Create anonymous functions for each body
        pos_sun = @(tq)  interp1(times{1}, orbits(:,:,1), jd, 'spline', 0);

        pos_triton = @(tq) interp1(times{2}, orbits(:,:,2), jd, 'spline', 0);
    end
    
    r1 = pos_sun(jd)';
    r2 = pos_triton(jd)';
    r_bodies = {r1, r2};

    % Compute total acceleration
    a_N_body = zeros(3,1);
    for i = 1:2
        delta_r = r_bodies{i} - r_sc;
        d = norm(delta_r);
        if d > 0
            a_N_body = a_N_body + mus(i) * delta_r / d^3;
        end
    end
end

%% Drag Pertubration

function a_drag = drag_acc(r, v, R, atmos, zonal_model, meridional_model, beta)
    omega_N = [0; 0; 2 * pi/(16.11 * 3600)];
    v_rel = v - cross(omega_N, r) - v_wind_inertial(r, R, zonal_model, meridional_model);
    rho = atmos(norm(r) - R);
    a_drag = -0.5 * rho * norm(v_rel) * v_rel / beta;
end

%% Perturbed equations of motion

function dY = eom_perturbed(t, Y, mu_N, R, R_equ, perturbations, jd0, orbits, times, mus, atmos, zonal, merid, beta)
    r = Y(1:3);
    v = Y(4:6);
    r_norm = norm(r);
    
    a_grav = -mu_N * r / r_norm^3;
    a_total = a_grav;

    for k = 1:length(perturbations)
        switch perturbations{k}
            case "drag"
                if r_norm - R <= 4000
                    a_drag = drag_acc(r, v, R, atmos, zonal, merid, beta);
                    a_total = a_total + a_drag;
                end
            case "J2"
                a_J2 = J2_acc(r, mu_N, R_equ);
                a_total = a_total + a_J2;
            case "N-body"
                jd = jd0 + t/86400;
                a_N_body = N_body_acc(jd, orbits, times, mus, r);
                a_total = a_total + a_N_body;
        end
    end
    
    dY = [v; a_total];
end

%% === Wind function ===
function v_wind_inertial = v_wind_inertial(r_vec, R, zonal_model, meridional_model)
    x = r_vec(1); y = r_vec(2); z = r_vec(3);
    r_mag = norm(r_vec);
    alt = r_mag - R;
    lat = asind(z / r_mag);
    lon = atan2d(y, x);
    lat_rad = deg2rad(lat);
    lon_rad = deg2rad(lon);
    e_hat = [-sind(lon_rad); cosd(lon_rad); 0];
    n_hat = [-sind(lat_rad)*cosd(lon_rad);
             -sind(lat_rad)*sind(lon_rad);
              cosd(lat_rad)];
    v_zonal = zonal_model(alt);
    v_meridional = meridional_model(alt);
    v_wind_inertial = v_zonal * e_hat + v_meridional * n_hat;
end