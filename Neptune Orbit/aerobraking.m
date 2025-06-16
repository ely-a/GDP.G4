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
m_N = 102.409e24; % kg
r_N = 24622; % km
r_N_equ = 24764; %km

%% Neptune atmosphere model

density_data = readmatrix("densitymodel.txt");
rhos = density_data(:, 1) * 1e9;  % kg/km^3
alts = density_data(:, 2) / 1e3;  % km

atmos = @(alt_km) interp1(alts, rhos, alt_km, 'linear', 0);

%% Initial orbit parameters

alt_cap = 10000;
rp_cap = alt_cap + r_N;
e_cap = 1;
vp_cap = sqrt(2 * mu_N/rp_cap);
i_cap = 10;
RAAN_cap = 40;
omega_cap = 50;

[rp_cap_vec, vp_cap_vec] = rv_parabolic(rp_cap, vp_cap, i_cap, RAAN_cap, omega_cap);

%% Load Science Orbit

load("ScienceOrbit.mat")

%% Finding apogee of capture orbit

dv_init = 0.01:0.01:5;
vp_init = vp_cap - dv_init;
rp_init = rp_cap;
eps_init = (vp_init.^2)/2 - mu_N/rp_init;
a_init = -mu_N./(2*eps_init);
e_init = 1 - rp_init./a_init;
ra_init = a_init.*(1+e_init);

figure
plot(dv_init, ra_init)
xlabel("$\Delta V$ (km s$^{-1})$", "Interpreter", "latex")
ylabel("r$_a$ (km)", "Interpreter","latex")
title("r$_a$ vs initial $\Delta V$ post capture", "Interpreter", "latex")

figure
plot(dv_init, e_init)
xlabel("$\Delta V$ (km s$^{-1})$", "Interpreter", "latex")
ylabel("r$_a$ (km)", "Interpreter","latex")
title("r$_a$ vs initial $\Delta V$ post capture", "Interpreter", "latex")

%% Chemical to rach science

ra_final = a_final * (1+e_final);
rp_final = a_final * (1-e_final);
vp_final = h_final/rp_final;
va_final = h_final/ra_final;

e_inter = (ra_final - rp_cap)/(ra_final + rp_cap); 
a_inter = 0.5 * (ra_final + rp_cap);
h_inter = sqrt(mu_N * a_inter * (1 - e_inter^2));
va_inter = h_inter/ra_final;
vp_inter = h_inter/rp_cap;
dv1 = vp_inter - vp_cap;
dv2 = va_final - va_inter;


%% Select initial dv based on trade-off

dv_init = 0.22;
vp_init = vp_cap - dv_init;
rp_init_vec = rp_cap_vec;
vp_init_vec = vp_cap_vec * (1 - dv_init/vp_cap);

[a_init, e_init, h_init, i_init, RAAN_init, omega_init, ~] = oe_from_rv(rp_init_vec, vp_init_vec, mu_N);
[ra_init_vec, va_init_vec] = rv_from_oe(a_init, e_init, i_init, RAAN_init, omega_init, 180, mu_N);

ra_init = norm(ra_init_vec);
va_init = norm(va_init_vec);
%% Delta v at apoapsis to reach atmosphere

ra_brake1 = ra_init;
alt_brake1 = 460:1:460;
rp_brake1 = r_N + alt_brake1;

e_brake1 = (ra_brake1 - rp_brake1)./(ra_brake1 + rp_brake1);
a_brake1 = 0.5 * (ra_brake1 + rp_brake1);
h_brake1 = sqrt(mu_N * a_brake1 .* (1 - e_brake1.^2));
va_brake1 = h_brake1./ra_brake1;

dv_brake1 = va_brake1 - va_init;

%% Aerobrake propogation for densty uncertatinities

max_atm_alt = alts(end);
max_atm_r = r_N + max_atm_alt;
density_factors = [0.9, 1, 1.1];
density_labels = {'-10% density', 'Nominal', '+10% density'};
colors = {'r', 'k', 'b'};

results = struct();

for dcase = 1:3
    fprintf('\n--- Examining %s ---\n', density_labels{dcase});
    % Adjust density model
    rhos_mod = rhos * density_factors(dcase);
    atmos_mod = @(alt_km) interp1(alts, rhos_mod, alt_km, 'linear', 0);

    % --- Aerobrake propagation ---
    state_exit = zeros(6, length(alt_brake1));
    event_type = zeros(1, length(alt_brake1));
    for idx = 1:length(alt_brake1)
        fprintf('  Periapsis altitude: %d km\n', alt_brake1(idx));
        theta_entry = acosd(((h_brake1(idx)^2 / (mu_N*max_atm_r)) - 1)/e_brake1(idx));
        theta_entry = 360 - theta_entry;
        [r_entry, v_entry] = rv_from_oe(a_brake1(idx), e_brake1(idx), i_init, RAAN_init, omega_init, theta_entry, mu_N);
        Y0 = [r_entry; v_entry];
        opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
                  'Events', @(t, Y) exit_atmosphere_event(t, Y, max_atm_r));
        [t_out, Y_out, te, ye, ie] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, atmos_mod), ...
                            [0, 0.5 * 365.25 * 86400], Y0, opts);
        state_exit(:,idx) = Y_out(end, :)';
        if ~isempty(ie)
            event_type(idx) = ie(end);
        else
            event_type(idx) = 0;
        end
    end

    % Filter successful escapes
    success_idx = find(event_type == 1);
    alt_success = alt_brake1(success_idx);
    state_success = state_exit(:, success_idx);

    % --- Post-processing: orbital elements ---
    N_success = size(state_success, 2);
    apoapsis = zeros(1, N_success);
    eccentricity = zeros(1, N_success);
    inclination = zeros(1, N_success);
    RAAN = zeros(1, N_success);
    omega = zeros(1, N_success);

    for k = 1:N_success
        r_vec = state_success(1:3, k);
        v_vec = state_success(4:6, k);
        [a, e, ~, i, RAAN_k, omega_k, ~] = oe_from_rv(r_vec, v_vec, mu_N);
        apoapsis(k) = a * (1 + e); % km
        eccentricity(k) = e;
        inclination(k) = i;
        RAAN(k) = RAAN_k;
        omega(k) = omega_k;
    end

    % Store results
    results(dcase).alt_success = alt_success;
    results(dcase).apoapsis = apoapsis;
    results(dcase).eccentricity = eccentricity;
    results(dcase).inclination = inclination;
    results(dcase).RAAN = RAAN;
    results(dcase).omega = omega;
end

%% Plot orbital elements against aiminig periapsis
figure;
hold on;
for dcase = 1:3
    plot(results(dcase).alt_success, results(dcase).apoapsis, 'Color', colors{dcase}, 'DisplayName', density_labels{dcase});
end
xlabel('Periapsis Altitude (km)');
ylabel('Apoapsis Radius (km)');
title('Apoapsis Radius vs Periapsis Altitude');
legend;
hold off;

figure;
subplot(2,2,1); hold on;
for dcase = 1:3
    plot(results(dcase).alt_success, results(dcase).eccentricity, 'Color', colors{dcase}, 'DisplayName', density_labels{dcase});
end
xlabel('Periapsis Altitude (km)'); ylabel('Eccentricity'); title('Eccentricity');
hold off;

subplot(2,2,2); hold on;
for dcase = 1:3
    plot(results(dcase).alt_success, results(dcase).inclination, 'Color', colors{dcase});
end
xlabel('Periapsis Altitude (km)'); ylabel('Inclination (deg)'); title('Inclination');
hold off;

subplot(2,2,3); hold on;
for dcase = 1:3
    plot(results(dcase).alt_success, results(dcase).RAAN, 'Color', colors{dcase});
end
xlabel('Periapsis Altitude (km)'); ylabel('RAAN (deg)'); title('RAAN');
hold off;

subplot(2,2,4); hold on;
for dcase = 1:3
    plot(results(dcase).alt_success, results(dcase).omega, 'Color', colors{dcase});
end
xlabel('Periapsis Altitude (km)'); ylabel('Argument of Periapsis (deg)'); title('\omega');
hold off;

sgtitle('Post-Aerobrake Orbital Elements vs Periapsis Altitude');

%% Plot chosen altitude

% --- User input: choose periapsis altitude (must be in alt_brake1) ---
chosen_alt = 460; % <-- Set this to your desired periapsis altitude (km)
[~, idx_chosen] = min(abs(alt_brake1 - chosen_alt)); % Find closest index

% Choose which density case to use (1: -10%, 2: nominal, 3: +10%)
chosen_case = 2; % 2 = nominal

% Re-run propagation for this case and periapsis
rhos_mod = rhos * density_factors(chosen_case);
atmos_mod = @(alt_km) interp1(alts, rhos_mod, alt_km, 'linear', 0);

theta_entry = acosd(((h_brake1(idx_chosen)^2 / (mu_N*max_atm_r)) - 1)/e_brake1(idx_chosen));
theta_entry = 360 - theta_entry;
[r_entry, v_entry] = rv_from_oe(a_brake1(idx_chosen), e_brake1(idx_chosen), i_init, RAAN_init, omega_init, theta_entry, mu_N);
Y0 = [r_entry; v_entry];
opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
          'Events', @(t, Y) exit_atmosphere_event(t, Y, max_atm_r));
[t_out, Y_out] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, atmos_mod), ...
                        [0, 2e4], Y0, opts);

% Calculate acceleration and altitude over time
r_vecs = Y_out(:,1:3);
v_vecs = Y_out(:,4:6);
N_pts = size(Y_out,1);
accel = zeros(N_pts,1);
altitude = zeros(N_pts,1);
for n = 1:N_pts
    a_vec = eom_perturbed(0, Y_out(n,:)', mu_N, r_N, atmos_mod);
    accel(n) = norm(a_vec(4:6));
    altitude(n) = norm(r_vecs(n,:)) - r_N;
end

% --- Plot acceleration vs time ---
figure;
plot(t_out/60, accel);
xlabel('Time (min)');
ylabel('Acceleration (km/s^2)');
title(['Total Acceleration vs Time (Periapsis Altitude = ', num2str(alt_brake1(idx_chosen)), ' km, ', density_labels{chosen_case}, ')']);

% --- Plot velocity vs time ---
velocity = vecnorm(v_vecs, 2, 2); % km/s
figure;
plot(t_out/60, velocity);
xlabel('Time (min)');
ylabel('Velocity (km/s)');
title(['Velocity vs Time (Periapsis Altitude = ', num2str(alt_brake1(idx_chosen)), ' km, ', density_labels{chosen_case}, ')']);

% --- Plot altitude vs time ---
figure;
plot(t_out/60, altitude);
xlabel('Time (min)');
ylabel('Altitude above Neptune (km)');
title(['Altitude vs Time (Periapsis Altitude = ', num2str(alt_brake1(idx_chosen)), ' km, ', density_labels{chosen_case}, ')']);

% --- Plot 3D trajectory around Neptune ---
figure;
hold on;
rs = r_vecs';
plot3(rs(1,:), rs(2,:), rs(3,:), 'b', 'LineWidth', 2);
[X, Y, Z] = sphere(100);  % Neptune for reference
surf(r_N*X, r_N*Y, r_N*Z, ...
    'FaceColor', [0.2 0.4 1], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.4);
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
title(['3D Trajectory (Periapsis Altitude = ', num2str(alt_brake1(idx_chosen)), ' km, ', density_labels{chosen_case}, ')']);
axis equal;
grid on;
hold off;
view(3);

aero_data = [t_out/60, velocity, altitude];
writematrix(aero_data, 'aerobraking_data.txt');
%% Alternative aerobrake 1 propogation (commented out)
% state_exit = zeros(6, length(alt_brake1));
% event_type = zeros(1, length(alt_brake1)); % 1: exited atmosphere, 2: impacted Neptune, 0: no event

% for idx = 1:length(alt_brake1)
%     disp(strcat("Aiming Periapsis: ", num2str(alt_brake1(idx)), " km"))
%     theta_entry = acosd(((h_brake1(idx)^2 / (mu_N*max_atm_r)) - 1)/e_brake1(idx));
%     theta_entry = 360 - theta_entry;
%     [r_entry, v_entry] = rv_from_oe(a_brake1(idx), e_brake1(idx), i_init, RAAN_init, omega_init, theta_entry, mu_N);
    
%     Y0 = [r_entry; v_entry];

%     opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
%               'Events', @(t, Y) exit_atmosphere_event(t, Y, max_atm_r));

%     [t_out, Y_out, te, ye, ie] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, atmos), ...
%                         [0, 2e4], Y0, opts);  % 1day as placeholder max time
    
%     state_exit(:,idx) = Y_out(end, :)';
    
%     if ~isempty(ie)
%         switch ie(end)
%             case 1
%                 fprintf('Exited atmosphere at t = %.2f s, altitude = %.2f km (event 1)\n', te(end), norm(ye(end,1:3)) - r_N);
%             case 2
%                 fprintf('Impacted Neptune at t = %.2f s, radius = %.2f km (event 2)\n', te(end), norm(ye(end,1:3)));
%             otherwise
%                 fprintf('Exited for unknown reason (event %d) at t = %.2f s\n', ie(end), te(end));
%         end
%     else
%         fprintf('Integration ended without triggering an event.\n');
%     end

%     if ~isempty(ie)
%         event_type(idx) = ie(end);
%     else
%         event_type(idx) = 0;
%     end
% end

% % figure
% % hold on
% % rs = Y_out(:, 1:3)';
% % plot3(rs(1,:), rs(2,:), rs(3, :))
% % [X, Y, Z] = sphere(100);  % resolution = 100
% % surf(r_N*X, r_N*Y, r_N*Z, ...
% %     'FaceColor', [0.2 0.4 1], ...
% %     'EdgeColor', 'none', ...
% %     'FaceAlpha', 0.4);

% success_idx = find(event_type == 1); % indices of successful escapes
% alt_success = alt_brake1(success_idx); % corresponding periapsis altitudes
% state_success = state_exit(:, success_idx); % corresponding final states

% % Example: print or plot successful periapsis altitudes
% disp('Successful escape periapsis altitudes (km):');
% disp(alt_success);

% %% Investigate apoapsis post aerobrake 1

% N_success = size(state_success, 2);
% apoapsis = zeros(1, N_success);
% eccentricity = zeros(1, N_success);
% inclination = zeros(1, N_success);
% RAAN = zeros(1, N_success);
% omega = zeros(1, N_success);

% for k = 1:N_success
%     r_vec = state_success(1:3, k);
%     v_vec = state_success(4:6, k);
%     [a, e, ~, i, RAAN_k, omega_k, ~] = oe_from_rv(r_vec, v_vec, mu_N);
%     apoapsis(k) = a * (1 + e); % km
%     eccentricity(k) = e;
%     inclination(k) = i;
%     RAAN(k) = RAAN_k;
%     omega(k) = omega_k;
% end

% % Plot apoapsis radius vs periapsis altitude
% figure;
% plot(alt_success, apoapsis, 'o-');
% xlabel('Periapsis Altitude (km)');
% ylabel('Apoapsis Radius (km)');
% title('Apoapsis Radius vs Periapsis Altitude');

% % Plot other elements vs periapsis altitude
% figure;
% subplot(2,2,1)
% plot(alt_success, eccentricity, 'o-');
% xlabel('Periapsis Altitude (km)');
% ylabel('Eccentricity');
% title('Eccentricity');

% subplot(2,2,2)
% plot(alt_success, inclination, 'o-');
% xlabel('Periapsis Altitude (km)');
% ylabel('Inclination (deg)');
% title('Inclination');

% subplot(2,2,3)
% plot(alt_success, RAAN, 'o-');
% xlabel('Periapsis Altitude (km)');
% ylabel('RAAN (deg)');
% title('RAAN');

% subplot(2,2,4)
% plot(alt_success, omega, 'o-');
% xlabel('Periapsis Altitude (km)');
% ylabel('Argument of Periapsis (deg)');
% title('\omega');

% sgtitle('Post-Aerobrake Orbital Elements vs Periapsis Altitude');

%% Exit atmosphere detetction

function [value, isterminal, direction] = exit_atmosphere_event(~, Y, r_exit)
    r = Y(1:3);
    r_mag = norm(r);

    % Event 1: Exit atmosphere
    value(1) = r_mag - r_exit;   % Zero when spacecraft exits the atmosphere
    isterminal(1) = 1;           % Stop integration
    direction(1) = 1;            % Only trigger when radius is increasing (exit)

    % Event 2: Impact Neptune
    r_N = 24622;                 % Neptune radius [km] (should match your main code)
    value(2) = r_mag - r_N;      % Zero when spacecraft hits Neptune
    isterminal(2) = 1;           % Stop integration
    direction(2) = -1;           % Only trigger when radius is decreasing (impact)
end

%% Drag Perturbations

function a_drag = drag_acc(r, v, R, atmos)
    beta = 1800e6/(1.15 * pi * 2.5^2);

    sidereal_period_N = 16.11 * 3600; % 16.11 hours
    omega_N = [0;0;2 * pi/sidereal_period_N]; % rad/s
    v_atm = cross(omega_N, r);
    v_rel = v - v_atm;

    alt = norm(r) - R;
    rho = atmos(alt);

    a_drag = -0.5 * rho * norm(v_rel) * v_rel / beta;
end


%% Pertubed propogation

function dY = eom_perturbed(~, Y, mu, R, atmos)
    r = Y(1:3);
    v = Y(4:6);
    r_norm = norm(r);

    a_grav = -mu * r / r_norm^3;

    a_drag = drag_acc(r, v, R, atmos);

    a_total = a_grav + a_drag;

    dY = [v; a_total];
end


