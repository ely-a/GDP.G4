% filepath: c:\Users\ishaa\OneDrive - Imperial College London\Uni\Year 3\GDP\GDP.G4\Neptune Orbit\monte_carlo_heatshield.m
% Monte Carlo simulation for Neptune aerobraking heat shield safety factor

clear; clc; close all

%% === SETTINGS ===
N_MC = 100; % Number of Monte Carlo runs

% Load atmospheric and wind statistics
load("gram_profiles.mat"); % rho_mean, rho_std, ew_mean, ew_std, ns_mean, ns_std, alt_steps

% Nominal values (edit as appropriate)
mu_N = 6.8351e6;         % km^3/s^2
r_N = 24622;             % km (mean radius)
r_N_equ = 24764;         % km (equatorial radius)
C_D_nom = 1.22; % or your nominal value
C_D_std = 0.1 * C_D_nom;
A = pi * 3^2; % m^2, or your reference area
m = 1740e6;     % kg, or your nominal mass
alt_entry_nominal = 581;                   % km, as in your script
alt_entry_std = 5;                         % km, 1-sigma

% New: Standard deviations for orbital elements
i_cap_nom = 65.0; i_cap_std = 0.5;             % deg
RAAN_cap_nom = 24.8; RAAN_cap_std = 0.05;       % deg
omega_cap_nom = 99.8; omega_cap_std = 0.1;    % deg

% Preallocate results
max_m_TPS = zeros(N_MC,1);
max_t_TPS = zeros(N_MC,1);
time_taken_days = zeros(N_MC,1);
num_passes = zeros(N_MC,1);
min_m_TPS = zeros(N_MC,1);
min_t_TPS = zeros(N_MC,1);

n_failed = 0;

for mc = 1:N_MC
    % --- Sample density and wind profiles ---
    rho_mc = rho_mean + randn(size(rho_mean)) .* rho_std;
    ew_mc  = ew_mean  + randn(size(ew_mean))  .* ew_std;
    ns_mc  = ns_mean  + randn(size(ns_mean))  .* ns_std;
    temp_mc = temps_mean + randn(size(temps_mean)) .* temps_std;

    % --- Sample beta and entry altitude ---
    C_D_mc = C_D_nom + randn * C_D_std;
    beta_mc = m / (C_D_mc * A);
    alt_entry_mc = alt_entry_nominal + randn * alt_entry_std;

    % --- Sample orbital elements ---
    i_cap_mc = i_cap_nom + randn * i_cap_std;
    RAAN_cap_mc = RAAN_cap_nom + randn * RAAN_cap_std;
    omega_cap_mc = omega_cap_nom + randn * omega_cap_std;

    % --- Build interpolant functions for this run ---
    atmos_mc = @(alt_km) interp1(alt_steps, rho_mc * 1e9, alt_km, 'linear', 0); % kg/km^3
    zonal_mc = @(alt_km) interp1(alt_steps, ew_mc/1000, alt_km, 'linear', 0);   % km/s
    merid_mc = @(alt_km) interp1(alt_steps, ns_mc/1000, alt_km, 'linear', 0);   % km/s
    temps_mc = @(alt_km) interp1(alt_steps, temp_mc, alt_km, 'linear', 0);

    try
        [m_TPS_vec, t_TPS_cm_vec, t_days, pass] = run_aerobrake_sim(atmos_mc, zonal_mc, merid_mc, temps_mc, beta_mc, alt_entry_mc, mu_N, r_N, r_N_equ, i_cap_mc, RAAN_cap_mc, omega_cap_mc);
        max_m_TPS(mc) = max(m_TPS_vec);
        max_t_TPS(mc) = max(t_TPS_cm_vec);
        time_taken_days(mc) = t_days;
        num_passes(mc) = pass;
        min_m_TPS(mc) = min(m_TPS_vec);
        min_t_TPS(mc) = min(t_TPS_cm_vec);
    catch
        max_m_TPS(mc) = NaN;
        max_t_TPS(mc) = NaN;
        time_taken_days(mc) = NaN;
        num_passes(mc) = NaN;
        n_failed = n_failed + 1;
        min_m_TPS(mc) = NaN;
        min_t_TPS(mc) = NaN;
        fprintf('MC run %d FAILED\n', mc);
    end
    fprintf('Completed MC run %d of %d\n', mc, N_MC);
end
fprintf('Total failed MC runs: %d out of %d\n', n_failed, N_MC);

% Remove failed runs
max_m_TPS = max_m_TPS(~isnan(max_m_TPS));
max_t_TPS = max_t_TPS(~isnan(max_t_TPS));
min_m_TPS = min_m_TPS(~isnan(min_m_TPS));
min_t_TPS = min_t_TPS(~isnan(min_t_TPS));
valid = ~isnan(max_m_TPS) & ~isnan(time_taken_days) & ~isnan(num_passes);
time_taken_days = time_taken_days(valid);
num_passes = num_passes(valid);


%% === Plot Results ===

% Plot histogram for max heat shield mass
figure;
histogram(max_m_TPS, 25, 'FaceColor', [0.2 0.6 1], 'Normalization', 'pdf');
xlabel('Maximum Heat Shield Mass per MC Run (kg)');
ylabel('Probability Density');
title('Monte Carlo: Maximum Heat Shield Mass Distribution');
grid on;
hold on;
pd1 = fitdist(max_m_TPS, 'Normal');
x1 = linspace(min(max_m_TPS), max(max_m_TPS), 100);
y1 = normpdf(x1, pd1.mu, pd1.sigma);
plot(x1, y1, 'r-', 'LineWidth', 2);
legend('Histogram', 'Normal Fit');
hold off;

% Plot histogram for max heat shield thickness
figure;
histogram(max_t_TPS, 25, 'FaceColor', [1 0.6 0.2], 'Normalization', 'pdf');
xlabel('Maximum Heat Shield Thickness per MC Run (cm)');
ylabel('Probability Density');
title('Monte Carlo: Maximum Heat Shield Thickness Distribution');
grid on;
hold on;
pd2 = fitdist(max_t_TPS, 'Normal');
x2 = linspace(min(max_t_TPS), max(max_t_TPS), 100);
y2 = normpdf(x2, pd2.mu, pd2.sigma);
plot(x2, y2, 'r-', 'LineWidth', 2);
legend('Histogram', 'Normal Fit');
hold off;

% Plot histogram for time taken (days)
figure;
histogram(time_taken_days, 25, 'FaceColor', [0.3 0.7 0.3], 'Normalization', 'pdf');
xlabel('Time Taken (days)');
ylabel('Probability Density');
title('Monte Carlo: Time Taken Distribution');
grid on;
hold on;
pd3 = fitdist(time_taken_days, 'Normal');
x3 = linspace(min(time_taken_days), max(time_taken_days), 100);
y3 = normpdf(x3, pd3.mu, pd3.sigma);
plot(x3, y3, 'r-', 'LineWidth', 2);
legend('Histogram', 'Normal Fit');
hold off;

% Plot histogram for number of passes
figure;
histogram(num_passes, 25, 'FaceColor', [0.7 0.3 0.7], 'Normalization', 'pdf');
xlabel('Number of Passes');
ylabel('Probability Density');
title('Monte Carlo: Number of Passes Distribution');
grid on;
hold on;
pd4 = fitdist(num_passes, 'Normal');
x4 = linspace(min(num_passes), max(num_passes), 100);
y4 = normpdf(x4, pd4.mu, pd4.sigma);
plot(x4, y4, 'r-', 'LineWidth', 2);
legend('Histogram', 'Normal Fit');
hold off;

% Plot histogram for min heat shield mass
figure;
histogram(min_m_TPS, 25, 'FaceColor', [0.2 0.3 0.8], 'Normalization', 'pdf');
xlabel('Minimum Heat Shield Mass per MC Run (kg)');
ylabel('Probability Density');
title('Monte Carlo: Minimum Heat Shield Mass Distribution');
grid on;
hold on;
pd_min1 = fitdist(min_m_TPS, 'Normal');
x_min1 = linspace(min(min_m_TPS), max(min_m_TPS), 100);
y_min1 = normpdf(x_min1, pd_min1.mu, pd_min1.sigma);
plot(x_min1, y_min1, 'r-', 'LineWidth', 2);
legend('Histogram', 'Normal Fit');
hold off;

% Plot histogram for min heat shield thickness
figure;
histogram(min_t_TPS, 25, 'FaceColor', [0.8 0.5 0.2], 'Normalization', 'pdf');
xlabel('Minimum Heat Shield Thickness per MC Run (cm)');
ylabel('Probability Density');
title('Monte Carlo: Minimum Heat Shield Thickness Distribution');
grid on;
hold on;
pd_min2 = fitdist(min_t_TPS, 'Normal');
x_min2 = linspace(min(min_t_TPS), max(min_t_TPS), 100);
y_min2 = normpdf(x_min2, pd_min2.mu, pd_min2.sigma);
plot(x_min2, y_min2, 'r-', 'LineWidth', 2);
legend('Histogram', 'Normal Fit');
hold off;

figure
subplot(1,2,1)
histogram(max_m_TPS, 25, 'FaceColor', [0.2 0.6 1], 'Normalization', 'pdf');
xlabel('Maximum Heat Shield Mass per MC Run (kg)');
ylabel('Probability Density');
grid on;
hold on;
pd1 = fitdist(max_m_TPS, 'Normal');
x1 = linspace(min(max_m_TPS), max(max_m_TPS), 100);
y1 = normpdf(x1, pd1.mu, pd1.sigma);
plot(x1, y1, 'r-', 'LineWidth', 2);
hold off;

subplot(1,2,2)
histogram(time_taken_days, 25, 'FaceColor', [0.3 0.7 0.3], 'Normalization', 'pdf');
xlabel('Time Taken (days)');
ylabel('Probability Density');
grid on;
hold on;
pd3 = fitdist(time_taken_days, 'Normal');
x3 = linspace(min(time_taken_days), max(time_taken_days), 100);
y3 = normpdf(x3, pd3.mu, pd3.sigma);
plot(x3, y3, 'r-', 'LineWidth', 2);
hold off;

fprintf('\n--- Monte Carlo Statistics ---\n');
fprintf('Max heat shield mass:      Mean = %.2f kg, 99th percentile = %.2f kg\n', mean(max_m_TPS), prctile(max_m_TPS,99));
fprintf('Max heat shield thickness: Mean = %.2f cm, 99th percentile = %.2f cm\n', mean(max_t_TPS), prctile(max_t_TPS,99));
fprintf('Min heat shield mass:      Mean = %.2f kg, 99th percentile = %.2f kg\n', mean(min_m_TPS), prctile(min_m_TPS,99));
fprintf('Min heat shield thickness: Mean = %.2f cm, 99th percentile = %.2f cm\n', mean(min_t_TPS), prctile(min_t_TPS,99));
fprintf('Time taken:                Mean = %.2f days, 99th percentile = %.2f days\n', mean(time_taken_days), prctile(time_taken_days,99));
fprintf('Number of passes:          Mean = %.2f, 99th percentile = %.2f\n', mean(num_passes), prctile(num_passes,99));

% Define threshold for failure (1 year)
duration_limit = 365.25; % days

% Compute the fraction of runs that exceed the limit
n_exceed = sum(time_taken_days > duration_limit);
fraction_exceed = n_exceed / length(time_taken_days);
percentile_exceed = 100 * (1 - fraction_exceed);

fprintf('\n%.2f%% of runs completed within 1 year.\n', 100 - 100*fraction_exceed);
fprintf('=> Aerobraking succeeds up to the %.2fth percentile.\n', percentile_exceed);


%% === Helper function: run_aerobrake_sim ===
function [m_TPS_vec, t_TPS_cm_vec, t_days, num_passes] = run_aerobrake_sim(atmos, zonal, merid, temps, beta, alt_entry, mu_N, r_N, r_N_equ, i_cap, RAAN_cap, omega_cap)
    % Initial orbit setup (as in your main script)                        
    alt_cap = 10000;                      % km above surface
    rp_cap = alt_cap + r_N;
    e_cap = 1;
    vp_cap = sqrt(2 * mu_N/rp_cap);       % km/s
    % Use sampled i_cap, RAAN_cap, omega_cap
    [rp_cap_vec, vp_cap_vec] = rv_parabolic(rp_cap, vp_cap, i_cap, RAAN_cap, omega_cap);

    dv_init = 0.22;                               % km/s
    vp_init = vp_cap - dv_init;
    rp_init_vec = rp_cap_vec;
    vp_init_vec = vp_cap_vec * (1 - dv_init/vp_cap);
    [a_init, e_init, h_init, i_init, RAAN_init, omega_init, ~] = oe_from_rv(rp_init_vec, vp_init_vec, mu_N);
    [ra_init_vec, va_init_vec] = rv_from_oe(a_init, e_init, i_init, RAAN_init, omega_init, 180, mu_N);
    ra_init = norm(ra_init_vec);
    va_init = norm(va_init_vec);

    % Use the sampled entry altitude (fixed for all passes)
    alt_brake = alt_entry;
    rp_brake = r_N + alt_brake;
    ra_brake = ra_init;
    e_brake = (ra_brake - rp_brake)./(ra_brake + rp_brake);
    a_brake = 0.5 * (ra_brake + rp_brake);
    h_brake = sqrt(mu_N * a_brake .* (1 - e_brake.^2));
    va_pre_brake = h_brake/ra_brake;
    va_post_brake = va_init;
    i_brake = i_init; RAAN_brake = RAAN_init; omega_brake = omega_init;
    [ra_brake_init, va_brake_init] = rv_from_oe(a_brake, e_brake, i_brake, RAAN_brake, omega_brake, 180, mu_N);
    Y0 = [ra_brake_init;va_brake_init];

    % Integrate until e < 0.45 at apogee (or your science orbit condition)
    opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
                  'Events', @(t, Y) stop_when_reached_science(t, Y, mu_N));
    [t_out, Y_out] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, zonal, merid, beta, ["drag"]), ...
                            [0, 2.5 * 365.25 * 86400], Y0, opts);

    % Post-process to find passes and heat shield requirements
    altitudes = vecnorm(Y_out(:,1:3), 2, 2) - r_N;
    in_atmos = altitudes <= 4000;
    entry_idx = find(diff(in_atmos) == 1);         % Just before entry
    exit_idx = find(diff(in_atmos) == -1) + 1;     % Just after exit
    N = min(length(entry_idx), length(exit_idx));
    m_TPS_vec = zeros(1, N);
    t_TPS_cm_vec = zeros(1, N);

    for i = 1:N
        idx_entry = entry_idx(i)+1; % first entry within atmosphere
        idx_exit = exit_idx(i)-1;   % last point within atmosphere

        % Extract profiles within the atmosphere for this pass
        t_pass = t_out(idx_entry:idx_exit) - t_out(idx_entry); % seconds, start from 0
        alt_pass = altitudes(idx_entry:idx_exit);              % km
        v_pass = vecnorm(Y_out(idx_entry:idx_exit,4:6),2,2);   % km/s

        % Call HeatShieldMass
        [m_TPS, t_TPS_cm] = HeatShieldMass('tiles' ,alt_pass, v_pass, t_pass, temps, atmos);

        m_TPS_vec(i) = m_TPS;
        t_TPS_cm_vec(i) = t_TPS_cm;
    end

    t_days = t_out(end)/(24*3600);
    num_passes = N;
end

%% === Event function for apogee ===
function [value, isterminal, direction] = stop_when_reached_science(~, Y, mu)
    r = Y(1:3); v = Y(4:6);
    [~, e, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu);

    % Condition 1: e < 0.45 and at apogee
    cond1 = e - 0.45;
    cond2 = abs(mod(theta,360) - 180) - 10;
    event1 = max([cond1, cond2]); % triggers when both < 0

    % Condition 2: collision with Neptune (r <= r_N)
    r_N = 24622; % km (should match your main code)
    r_mag = norm(r);
    event2 = r_mag - r_N; % triggers when <= 0

    value = [event1, event2];
    isterminal = [1, 1];
    direction = [-1, -1];
end
%% === Equations of motion with drag ===
function dY = eom_perturbed(~, Y, mu, R, R_equ, atmos, zonal, merid, beta, perturbations)
    r = Y(1:3); v = Y(4:6);
    r_norm = norm(r);
    a_total = -mu * r / r_norm^3;

    for k = 1:length(perturbations)
        switch perturbations{k}
            case "drag"
                if r_norm - R <= 4000
                    a_total = a_total + drag_acc(r, v, R, atmos, zonal, merid, beta);
                end
        end
    end

    dY = [v; a_total];
end

%% === Drag acceleration with variable beta ===
function a_drag = drag_acc(r, v, R, atmos, zonal_model, meridional_model, beta)
    omega_N = [0; 0; 2 * pi/(16.11 * 3600)];
    v_rel = v - cross(omega_N, r) - v_wind_inertial(r, R, zonal_model, meridional_model);
    rho = atmos(norm(r) - R);
    a_drag = -0.5 * rho * norm(v_rel) * v_rel / beta;
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