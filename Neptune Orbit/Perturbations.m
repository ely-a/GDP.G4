clear
clc
close all

%% SET GRAPHING FONTS AND SIZES
set(groot,'defaultLineLineWidth',3) % Set all line widths to 2 unless stated otherwise
set(groot,'defaultAxesFontSize',20) % Set all axes font size to 16 unless stated otherwise
set(groot,'defaulttextfontsize',20) % Set all text font sizes to 16 unless stated otherwise
set(groot,'defaultLineMarkerSize',12) % Set all marker sizes to 16 unless stated otherwise
set(groot,'defaultAxesXGrid','on') % Turn x grid lines on
set(groot,'defaultAxesYGrid','on') % Turn y grid lines on

%% Neptune Constants

mu_N = 6.8351e6; % km^3/s^2
mass_neptune = 102.409e24; % kg
r_N = 24622;             % km (mean radius)
r_N_equ = 24764;         % km (equatorial radius)

%% Initialise glocal variables
global acc_log
acc_log = struct('time', [], 'a_grav', [], 'a_drag', [], 'a_J2', [], 'a_Nbody', []);

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

%% ULoad Science Orbit

if exist('ScienceOrbit.mat', 'file') == 2
    load("ScienceOrbit.mat")
else
    error('no intitial conditions')
end

%% Uneprturbed Trajectory

% [r0, v0] = rv_from_oe(3e5, 0.9, 90, 90, 250, 0, mu_neptune);
[r0, v0] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, 180, mu_N);
[a, e, h, i, Omega, omega, theta] = oe_from_rv(r0, v0, mu_N);
P = 2 * pi * sqrt(a^3 / mu_N);  % seconds
N_orbits = 1;
N_orbits_total = (5 * 365.25 - t_days) / (P / 86400);  % number of orbits in 5 years
% T = P/100 * N_orbits_total;
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
perturbations = ["J2", "drag", "N-body"];

orbits_sun_triton = orbits_mat(:,:,[idx_sun, idx_triton]);
times_sun_triton = {times_all{idx_sun}, times_all{idx_triton}};
mu_sun = 1.327e11;
mu_triton = 1428.495;
mu_sun_triton = [mu_sun, mu_triton];

opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
              'Events', @(t,Y) stop_at_apoapsis(t, Y, N_orbits));

[t_out, Y_out] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, ...
                                perturbations, jd0, orbits_sun_triton, ...
                                times_sun_triton, mu_sun_triton, ...
                                atmos_nom, zonal_nom, merid_nom, beta_nom), ...
                                tspan, Y0, opts);


t_plot = acc_log.time / 3600;  % convert to days

figure;
semilogy(t_plot, acc_log.a_grav * 1000, 'k-', 'DisplayName', "Neptune's Gravity");
hold on;
semilogy(t_plot, acc_log.a_drag * 1000, 'b:', 'DisplayName', 'Atmospheric Drag');
semilogy(t_plot, acc_log.a_J2 * 1000,   'r-.', 'DisplayName', '$J_2$ ');
semilogy(t_plot, acc_log.a_Nbody * 1000,'g--', 'DisplayName', 'N-body:Sun+Triton');

xlabel('Time [hours]', Interpreter='latex');
ylabel('Acceleration [m/s$^2$]', Interpreter='latex');
legend('Location', 'best', Interpreter='latex');
grid on;

axis([0 7 1e-8 1e2])

r_perturbed = Y_out(:,1:3)';
v_perturbed = Y_out(:,4:6)';

%% Plot trajectories

figure;
plot3(r_unperturbed(1,:), r_unperturbed(2,:), r_unperturbed(3,:), 'b-');  % unperturbed
hold on;
plot3(r_perturbed(1,:), r_perturbed(2,:), r_perturbed(3,:), 'r--', 'LineWidth', 1.5);  % perturbed

% Mark start and end points
% plot3(r_unperturbed(1,1), r_unperturbed(2,1), r_unperturbed(3,1), 'bo', 'MarkerFaceColor', 'b');
% text(r_unperturbed(1,1), r_unperturbed(2,1), r_unperturbed(3,1), '  Start', 'Color', 'b');
% 
% plot3(r_unperturbed(1,end), r_unperturbed(2,end), r_unperturbed(3,end), 'bs');
% text(r_unperturbed(1,end), r_unperturbed(2,end), r_unperturbed(3,end), '  End (Unperturbed)', 'Color', 'b');
% 
% plot3(r_perturbed(1,end), r_perturbed(2,end), r_perturbed(3,end), 'rs');
% text(r_perturbed(1,end), r_perturbed(2,end), r_perturbed(3,end), '  End (Perturbed)', 'Color', 'r');

legend('Unperturbed', 'Perturbed (J2 + Drag)', 'Location', 'best');
grid on;
axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

r0_perturb = r_perturbed(:,1);
v0_perturb = v_perturbed(:,1);

r_end_perturb = r_perturbed(:,end);
v_end_perturb = v_perturbed(:,end);

[a_0, e_0, ~, i_0, ~, ~ , ~] = oe_from_rv(r0_perturb, v0_perturb, mu_N);
[a_end, e_end, ~, i_end, ~, ~ ,theta_end] = oe_from_rv(r_end_perturb, v_end_perturb, mu_N);

theta_end
fprintf('variation in a:        %.6f km/day\n', (a_end - a_0)/(T/86400));
fprintf('variation in e:        %.6f /day\n', (e_end - e_0)/(T/86400));
fprintf('variation in i:        %.6f deg/day\n', (i_end - i_0)/(T/86400));


%% Orbital element variation
% Preallocate orbital elements over time
N = length(t_out);
a_vals      = zeros(1, N);
e_vals      = zeros(1, N);
i_vals      = zeros(1, N);
Omega_vals  = zeros(1, N);
omega_vals  = zeros(1, N);
theta_vals  = zeros(1, N);
h_vals = zeros(1, N);

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
    h_vals(k) = h_k;
end

% Convert time to days
t_days_plot = t_out / 86400;

% Plot
figure; 
subplot(3,2,1); plot(t_days_plot, a_vals);     ylabel('a [km]');  xlabel('Time [days]');
subplot(3,2,2); plot(t_days_plot, e_vals);       ylabel('e'); xlabel('Time [days]');
subplot(3,2,3); plot(t_days_plot, (i_vals));      ylabel('i [deg]'); xlabel('Time [days]');
subplot(3,2,6); plot(t_days_plot, (h_vals));  ylabel('h [km^2/s]'); xlabel('Time [days]');
subplot(3,2,4); plot(t_days_plot, (Omega_vals));  ylabel('\Omega [deg]'); xlabel('Time [days]');
subplot(3,2,5); plot(t_days_plot, (omega_vals));  ylabel('\omega [deg]'); xlabel('Time [days]');

%% J2 Omega and omega rate of change comparison
% === Analytical J2 Drift Rates ===

J2 = 3.411e-3;              % Neptune's J2 coefficient
mu = mu_N;           % km^3/s^2
R = r_N_equ;        % km

% Use initial orbital elements (from earlier in script)
n = sqrt(mu / a^3);        % mean motion [rad/s]

% RAAN drift (rad/s)
Omega_dot_rad = -1.5 * J2 * n * (R/a)^2 * cosd(i) / (1 - e^2)^2;

% Argument of perigee drift (rad/s)
omega_dot_rad = 0.75 * J2 * n * (R/a)^2 * (5 * cosd(i)^2 - 1) / (1 - e^2)^2;

% Convert to deg/day
Omega_dot_deg_day = rad2deg(Omega_dot_rad) * 86400;
omega_dot_deg_day = rad2deg(omega_dot_rad) * 86400;

fprintf('Analytical RAAN drift:       %.6f deg/day\n', Omega_dot_deg_day);
fprintf('Analytical perigee drift:    %.6f deg/day\n', omega_dot_deg_day);

RAAN_drift_sim = (Omega_vals(end) - Omega_vals(1)) / (t_days(end) - t_days(1));
omega_drift_sim = (omega_vals(end) - omega_vals(1)) / (t_days(end) - t_days(1));

fprintf('Simulated RAAN drift:        %.6f deg/day\n', RAAN_drift_sim);
fprintf('Simulated perigee drift:     %.6f deg/day\n', omega_drift_sim);

%% J2 reference trajectory

opts_J2 = odeset('RelTol',1e-9, 'AbsTol',1e-9);
[t_J2, Y_J2] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, ...
                                ["J2"], jd0, orbits_sun_triton, ...
                                times_sun_triton, mu_sun_triton, ...
                                atmos_nom, zonal_nom, merid_nom, beta_nom), ...
                                [0, P*10], Y0, opts_J2);

% Extract position and velocity from simulation
r_J2 = Y_J2(:,1:3)';
v_J2 = Y_J2(:,4:6)';

% Preallocate
N = size(Y_J2,1);
a_vals = zeros(1,N);
e_vals = zeros(1,N);
theta_vals = zeros(1,N);

for k = 1:N
    r_k = r_J2(:,k);
    v_k = v_J2(:,k);
    [a_k, e_k, ~, ~, ~, ~, theta_k] = oe_from_rv(r_k, v_k, mu_N);
    a_vals(k) = a_k;
    e_vals(k) = e_k;
    theta_vals(k) = mod(theta_k, 360);  % wrap to [0, 360)
end

% Sort by theta for smooth plotting
[theta_sorted, sort_idx] = sort(theta_vals);
a_sorted = a_vals(sort_idx);
e_sorted = e_vals(sort_idx);

% Plot a vs theta
figure;
subplot(2,1,1)
plot(theta_sorted, a_sorted, 'b-');
xlabel('\theta [deg]');
ylabel('a [km]');
title('Semi-Major Axis Variation with True Anomaly');
grid on;

% Plot e vs theta
subplot(2,1,2)
plot(theta_sorted, e_sorted, 'r-');
xlabel('\theta [deg]');
ylabel('e');
title('Eccentricity Variation with True Anomaly');
grid on;

a_a_target = a_final;
e_a_target = e_final;
rp_a_target = a_a_target * (1 - e_a_target);
ra_a_target = a_a_target * (1 + e_a_target);


a_p_target = interp1([theta_sorted(end), theta_sorted(1) + 360], [a_sorted(end), a_sorted(1)], 360, 'spline');
e_p_target = interp1([theta_sorted(end), theta_sorted(1) + 360], [e_sorted(end), e_sorted(1)], 360, 'spline');
rp_p_target = a_p_target * (1 - e_p_target);
ra_p_target = a_p_target * (1 + e_p_target);

%% Burns

burns = [1.446, 3.100, 1.440, 3.080, 1.424, 3.103, 1.408, 3.122, 1.447, 3.098, ...
         1.439, 3.095, 1.434, 3.092, 1.432, 3.115, 1.405, 3.137, 1.442, 3.107, ...
         1.442, 3.084, 1.434, 3.092, 1.431, 3.116, 1.425, 3.098, 1.432, 3.100, ...
         1.441, 3.079, 1.439, 3.096, 1.440, 3.101, 1.447, 3.097, 1.418, 3.115, ...
         1.443, 3.097, 1.406, 3.128, 1.411, 3.111, 1.431, 3.110, 1.442, 3.092];
length(burns)
sum(burns)
mean(burns(1:2:end))
mean(burns(2:2:end))

%% === TRAJECTORY CORRECTION MANOEUVRES LOOP ===
% a_thresh = 10;         % km
% % N_orbits_corr = N_orbits_total;   % Number of corrected orbits
% N_orbits_corr = N_orbits_total;
% T_total_corr = N_orbits_corr * P;
% 
% 
% % Initial state
% Y_corr = Y0;
% t_corr = 0;
% t_log_corr = 0;
% dV_log_corr = [];
% 
% opts_event_apo = odeset('RelTol',1e-9, 'AbsTol',1e-9, 'Events', @(t,Y) detect_apoapsis(t,Y, mu_N));
% opts_event_peri = odeset('RelTol',1e-9, 'AbsTol',1e-9, 'Events', @(t,Y) detect_periapsis(t,Y, mu_N));
% opts_J2_apo = odeset('RelTol',1e-9, 'AbsTol',1e-9, 'Events', @(t,Y) detect_apoapsis(t, Y, mu_N));
% opts_J2_peri = odeset('RelTol',1e-9, 'AbsTol',1e-9, 'Events', @(t,Y) detect_periapsis(t, Y, mu_N));
% 
% a_vary_a = [a_final];
% t_vary_a = [0];
% 
% while t_corr < T_total_corr
% 
%     % Integrate until next apoapsis
%     [t_seg, Y_seg] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, ...
%                                 perturbations, jd0+ t_corr/86400, orbits_sun_triton, ...
%                                 times_sun_triton, mu_sun_triton, ...
%                                 atmos_nom, zonal_nom, merid_nom, beta_nom), ...
%                                 [0, P * 1.1], Y_corr, opts_event_apo);
% 
%     % Final state at apoapsis
%     Y_end =  Y_seg(end,:)';
%     r_curr = Y_end(1:3);
%     v_curr = Y_end(4:6);
%     v_norm = v_curr/norm(v_curr);
%     % a_curr = norm(r_curr) / (2 - norm(r_curr)*norm(v_curr)^2 / mu_N);  % Vis-viva
%     [a_curr, e_curr, ~, ~, ~, ~, ~] = oe_from_rv(r_curr, v_curr, mu_N);
% 
%     a_vary_a(end+1) = a_curr;
%     t_corr = t_corr + t_seg(end);
% 
%     t_vary_a(end+1) = t_corr;
% 
%     % Check if correction needed
%     fprintf('a difference = %.4f km \n', abs(a_curr - a_a_target));
% 
%     if abs(a_curr - a_a_target) > a_thresh
%         % Correct periapsis at apoapsis
% 
%         dV_apo_high = 20;
%         dV_apo_low = -20;
%         peri_found = false;
% 
%         while peri_found == false
%             dV_apo = 0.5*(dV_apo_low + dV_apo_high);
%             fprintf('Apoapsis dv test: %.5f m/s\n', dV_apo * 1000);
% 
%             v_J2_apo = v_curr +  dV_apo * v_norm;
%             Y_J2_apo_0 = [r_curr; v_J2_apo];
%             [t_J2_apo, Y_J2_apo] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, ...
%                                 ["J2"], jd0, orbits_sun_triton, ...
%                                 times_sun_triton, mu_sun_triton, ...
%                                 atmos_nom, zonal_nom, merid_nom, beta_nom), ...
%                                 [0, P*0.7], Y_J2_apo_0, opts_J2_peri);
%             Y_J2_apo_rp = Y_J2_apo(end, :)';
%             r_p_apo_targeting = Y_J2_apo_rp(1:3);
%             r_p_mag = norm(r_p_apo_targeting);
% 
%             if abs(r_p_mag - rp_a_target) < 0.1
%                 peri_found = true;
%             elseif r_p_mag > rp_a_target
%                 dV_apo_high = dV_apo;
%             elseif r_p_mag < rp_a_target
%                 dV_apo_low = dV_apo;
%             end
%         end
% 
%         Y_corr = Y_J2_apo_rp;
% 
%         % Log correction
%         dV_log_corr(end+1,1) = dV_apo;
%         t_log_corr(end+1,1) = t_corr;
% 
%         fprintf('Apoapsis Manouvre: %.3f m/s\n', dV_apo * 1000);
% 
%         [t_seg, Y_seg] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, ...
%                                 perturbations, jd0 + t_corr/86400, orbits_sun_triton, ...
%                                 times_sun_triton, mu_sun_triton, ...
%                                 atmos_nom, zonal_nom, merid_nom, beta_nom), ...
%                                 [0, P * 0.7], Y_corr, opts_event_peri);
% 
%         Y_end =  Y_seg(end,:)';
%         r_curr = Y_end(1:3);
%         v_curr = Y_end(4:6);
%         v_norm = v_curr/norm(v_curr);
% 
%         t_corr = t_corr + t_seg(end);
% 
%         dV_peri_high = 20;
%         dV_peri_low = -20;
%         apo_found = false;
% 
%         while apo_found == false
%             dV_peri = 0.5*(dV_peri_low + dV_peri_high);
%             fprintf('Periapsis dv test: %.5f m/s\n', dV_peri * 1000);
% 
%             v_J2_peri = v_curr +  dV_peri * v_norm;
%             Y_J2_peri_0 = [r_curr; v_J2_peri];
%             [t_J2_apo, Y_J2_peri] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, ...
%                                 ["J2"], jd0, orbits_sun_triton, ...
%                                 times_sun_triton, mu_sun_triton, ...
%                                 atmos_nom, zonal_nom, merid_nom, beta_nom), ...
%                                 [0, P*0.7], Y_J2_peri_0, opts_J2_apo);
%             Y_J2_peri_ra = Y_J2_peri(end, :)';
%             r_a_peri_targeting = Y_J2_peri_ra(1:3);
%             r_a_mag = norm(r_a_peri_targeting);
% 
%             if abs(r_a_mag - ra_p_target) < 0.1
%                 apo_found = true;
%             elseif r_a_mag > ra_p_target
%                 dV_peri_high = dV_peri;
%             elseif r_a_mag < ra_p_target
%                 dV_peri_low = dV_peri;
%             end
%         end
% 
%         Y_corr = Y_J2_peri_ra;
% 
%         % Log correction
%         dV_log_corr(end+1,1) = dV_peri;
%         t_log_corr(end+1,1) = t_corr;
% 
%         fprintf('Periapsis Manouvre: %.3f m/s\n', dV_peri * 1000);
% 
%     else
%         Y_corr = Y_end;
%     end
% end
% 
% fprintf('\n=== CORRECTION SUMMARY ===\n');
% fprintf('Number of burns: %d\n', numel(dV_log_corr));
% fprintf('Total ΔV spent : %.4f m/s\n', sum(abs(dV_log_corr)) * 1000);

figure
plot(t_vary_a, a_vary_a)

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
    persistent pos_sun; 
    persistent pos_triton;

    if isempty(pos_sun)
        pos_sun = @(tq) interp1(times{1}, orbits(:,:,1), tq, 'spline', 0);
        pos_triton = @(tq) interp1(times{2}, orbits(:,:,2), tq, 'spline', 0);
    end

    r1 = pos_sun(jd)';
    r2 = pos_triton(jd)';
    r_bodies = {r1, r2};

    a_N_body = zeros(3,1);
    for i = 1:2
        r_i = r_bodies{i};              % position of body i relative to Neptune
        delta_r = r_i - r_sc;           % vector from spacecraft to body i
        d_sc = norm(delta_r);
        d_nep = norm(r_i);
        if d_sc > 0 && d_nep > 0
            term1 = mus(i) * delta_r / d_sc^3;
            term2 = mus(i) * r_i / d_nep^3;
            a_N_body = a_N_body + (term1 - term2);
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
    % global acc_log
    % 
    % a_drag = zeros(3,1);
    % a_J2   = zeros(3,1);
    % a_Nbody = zeros(3,1);

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
                a_Nbody = N_body_acc(jd, orbits, times, mus, r);
                a_total = a_total + a_Nbody;
        end
    end

    % acc_log.time(end+1,1)     = t;
    % acc_log.a_grav(end+1,1)   = norm(a_grav);
    % acc_log.a_drag(end+1,1)   = norm(a_drag);
    % acc_log.a_J2(end+1,1)     = norm(a_J2);
    % acc_log.a_Nbody(end+1,1)  = norm(a_Nbody);

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

function [value, isterminal, direction] = stop_at_apoapsis(t, Y, N_max)
    persistent last_r norm_r count
    
    r = Y(1:3);
    v = Y(4:6);
    norm_r = norm(r);
    
    % Initialise memory on first call
    if isempty(last_r)
        last_r = norm_r;
        count = 0;
    end

    % Detect local maxima (apoapsis)
    drdt = dot(r, v) / norm_r;
    value = drdt;            % Crosses zero at periapsis/apoapsis
    direction = -1;          % Detect only negative-going zero crossings (apoapsis)
    
    % Stop only if apoapsis has occurred N_max times
    if direction < 0 && value == 0
        count = count + 1;
    end
    isterminal = (count >= N_max); % Stop integration when true
end

function [value, isterminal, direction] = detect_apoapsis(~, Y, mu)
    % Y: [r; v] state vector
    % mu: gravitational parameter
    r = Y(1:3);
    v = Y(4:6);
    [~, ~, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu); % theta in degrees

    % Detect when theta crosses 180 deg (apoapsis)
    value = mod(theta, 360) - 179.99; % Crosses zero at apoapsis (180 degrees)
    isterminal = 1;    % Stop the integration
    direction = 1;     % Only detect crossing in positive direction
end

function [value, isterminal, direction] = detect_periapsis(~, Y, mu)
    r = Y(1:3); v = Y(4:6);
    [~, ~, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu); % theta in degrees
    value = sin(deg2rad(theta)); % triggers at theta = 0°
    isterminal = 1;
    direction = 1; % Only detect when theta increases through 0
end