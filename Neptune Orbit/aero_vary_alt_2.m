clear; clc; close all

%% === SETTINGS ===
% Load atmospheric and wind statistics
load("gram_profiles.mat"); % rho_mean, ew_mean, ns_mean, alt_steps

% Nominal values
mu_N = 6.8351e6;         % km^3/s^2
r_N = 24622;             % km (mean radius)
r_N_equ = 24764;         % km (equatorial radius)
C_D_nom = 1.22;
A = pi * 3^2;          % m^2
m = 1740e6;              % kg
beta_nom = m / (C_D_nom * A);
alt_entry_nominal = 581; % km

% Build interpolant functions for nominal run
atmos_nom = @(alt_km) interp1(alt_steps, rho_mean * 1e9, alt_km, 'linear', 0); % kg/km^3
zonal_nom = @(alt_km) interp1(alt_steps, ew_mean/1000, alt_km, 'linear', 0);   % km/s
merid_nom = @(alt_km) interp1(alt_steps, ns_mean/1000, alt_km, 'linear', 0);   % km/s
temps_nom = @(alt_km) interp1(alt_steps, temps_mean, alt_km, "linear", 0);

%% --- Sweep over periapsis altitudes ---
% alt_range = 540:90:620; % Example range in km
% max_m_TPS = zeros(size(alt_range));
% time_taken_days = zeros(size(alt_range));
% num_passes_vec = zeros(size(alt_range));
% 
% for j = 1:length(alt_range)
%     alt_entry = alt_range(j);
%     fprintf('Running for periapsis altitude: %d km (%d of %d)\n', alt_entry, j, length(alt_range));
%     [m_TPS_vec, ~, t_days, num_passes, ~] = run_aerobrake_nominal(atmos_nom, zonal_nom, merid_nom, temps_nom, beta_nom, alt_entry, mu_N, r_N, r_N_equ);
%     max_m_TPS(j) = max(m_TPS_vec);
%     time_taken_days(j) = t_days;
%     num_passes_vec(j) = num_passes;
% end
% 
% % --- Plot results ---
% figure;
% subplot(2,1,1)
% plot(alt_range, max_m_TPS, '-o', 'LineWidth', 2);
% xlabel('Periapsis Altitude (km)');
% ylabel('Max Heat Shield Mass (kg)');
% title('Max Heat Shield Mass vs Periapsis Altitude');
% grid on;
% 
% subplot(2,1,2)
% plot(alt_range, time_taken_days, '-o', 'LineWidth', 2);
% xlabel('Periapsis Altitude (km)');
% ylabel('Time to Science Orbit (days)');
% title('Time to Science Orbit vs Periapsis Altitude');
% grid on;

%% Run the nominal aerobrake simulation
[m_TPS_vec, t_TPS_cm_vec, t_days, num_passes, Y_out] = run_aerobrake_nominal(atmos_nom, zonal_nom, merid_nom, temps_nom, beta_nom, alt_entry_nominal, mu_N, r_N, r_N_equ);

r_traj = Y_out(:, 1:3); % Trajectory positions
% Display results
fprintf('Nominal aerobrake results:\n');
fprintf('  Number of passes: %d\n', num_passes);
fprintf('  Time taken: %.2f days\n', t_days);
fprintf('  Max heat shield mass: %.2f kg\n', max(m_TPS_vec));
SF = 1.4;
fprintf('  Max heat shield mass with SF = %.1f: %.2f kg\n', SF, max(m_TPS_vec) * SF);
fprintf('  Max heat shield thickness: %.2f cm\n', max(t_TPS_cm_vec));

% Plot heat shield mass and thickness per pass
figure;
subplot(2,1,1)
plot(1:num_passes, m_TPS_vec, '-', 'LineWidth', 2);
xlabel('Pass'); ylabel('Heat Shield Mass (kg)');
title('Heat Shield Mass per Pass'); grid on;
subplot(2,1,2)
plot(1:num_passes, t_TPS_cm_vec, '-', 'LineWidth', 2);
xlabel('Pass'); ylabel('Heat Shield Thickness (cm)');
title('Heat Shield Thickness per Pass'); grid on;

% Plot the 3D trajectory
figure;
plot3(r_traj(:,1), r_traj(:,2), r_traj(:,3), 'b', 'LineWidth', 1.5); hold on;
[Xs, Ys, Zs] = sphere(100);
surf(r_N*Xs, r_N*Ys, r_N*Zs, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title(sprintf('Aerobraking Trajectory (Periapsis Altitude = %d km)', alt_entry_nominal));
axis equal; grid on; view(3);
legend('Trajectory','Neptune');
hold off;

%% --- Plot orbital element variation after each pass ---

% Get indices for the end of each pass
altitudes = vecnorm(Y_out(:,1:3), 2, 2) - r_N;
in_atmos = altitudes <= 4000;
entry_idx = find(diff(in_atmos) == 1) + 1;         % First index inside atmosphere
exit_idx = find(diff(in_atmos) == -1);             % Last index inside atmosphere

N = min(length(entry_idx), length(exit_idx));
a_vec = zeros(1, N);
e_vec = zeros(1, N);
rp_vec = zeros(1, N);
ra_vec = zeros(1, N);

for i = 1:N
    idx = exit_idx(i); % State at the end of each pass
    r_pass = Y_out(idx,1:3)';
    v_pass = Y_out(idx,4:6)';
    [a, e, ~, ~, ~, ~, ~] = oe_from_rv(r_pass, v_pass, mu_N);
    a_vec(i) = a;
    e_vec(i) = e;
    rp_vec(i) = a * (1 - e);
    ra_vec(i) = a * (1 + e);
end

figure;
subplot(2,2,1)
plot(1:N, a_vec, '-o', 'LineWidth', 2);
xlabel('Pass'); ylabel('Semi-major axis (km)');
title('Semi-major axis after each pass'); grid on;

subplot(2,2,2)
plot(1:N, e_vec, '-o', 'LineWidth', 2);
xlabel('Pass'); ylabel('Eccentricity');
title('Eccentricity after each pass'); grid on;

subplot(2,2,3)
plot(1:N, rp_vec - r_N, '-o', 'LineWidth', 2);
xlabel('Pass'); ylabel('Periapsis Altitude (km)');
title('Periapsis Altitude after each pass'); grid on;

subplot(2,2,4)
plot(1:N, ra_vec - r_N, '-o', 'LineWidth', 2);
xlabel('Pass'); ylabel('Apoapsis Altitude (km)');
title('Apoapsis Altitude after each pass'); grid on;

%% Comaprison of final delta v with eccentricty

% --- Final science orbit insertion at apoapsis after aerobraking ---

% 1. Get the final state after aerobraking
r_final = Y_out(end,1:3)'; % Last position vector
v_final = Y_out(end,4:6)';

% 2. Get current orbital elements
[a_curr, e_curr, h_curr, i_curr, RAAN_curr, omega_curr, theta_curr] = oe_from_rv(r_final, v_final, mu_N);

[r_final, v_final] = rv_from_oe(a_curr, e_curr, i_curr, RAAN_curr, omega_curr, 180, mu_N);
% 3. Compute current apogee radi
% us
ra_curr = a_curr * (1 + e_curr);

% 4. Set desired periapsis for science orbit
rp_science = r_N + 1800; % km

% 5. Compute required eccentricity and semi-major axis for new orbit
e_science = (ra_curr - rp_science) / (ra_curr + rp_science);
a_science = 0.5 * (ra_curr + rp_science);

% 6. Compute velocity at apoapsis before and after maneuver
v_ap_before = norm(v_final);
v_ap_after = sqrt(mu_N * (2/ra_curr - 1/a_science));

% 7. Compute required delta-v at apoapsis
delta_v_apogee = v_ap_after - v_ap_before;

fprintf('\n--- Science Orbit Insertion ---\n');
fprintf('Current apogee radius: %.2f km\n', ra_curr);
fprintf('Target periapsis: %.2f km\n', rp_science);
fprintf('Required science orbit eccentricity: %.6f\n', e_science);
fprintf('Required science orbit semi-major axis: %.2f km\n', a_science);
fprintf('Velocity at apoapsis before maneuver: %.4f km/s\n', v_ap_before);
fprintf('Velocity at apoapsis after maneuver:  %.4f km/s\n', v_ap_after);
fprintf('Required delta-v at apoapsis:         %.4f m/s\n', delta_v_apogee*1000);

% Optionally, update the velocity vector for the new orbit:
v_final_after_impulse = v_final + (delta_v_apogee / norm(v_final)) * v_final;

[a_final, e_final, h_final, i_final, RAAN_final, omega_final, ~] = oe_from_rv(r_final, v_final_after_impulse, mu_N);

t_laser = find_time_neptune(r_N, h_final, e_final, a_final, mu_N, 2000)/60;
P = 2*pi * sqrt(a_final^3/mu_N)/3600;

fprintf('Laser run time/orbit:  %.2f minutes\n', t_laser);
fprintf('Science orbit period:  %.2f hours\n', P);
fprintf('Number of science orbits = %.1f\n', (5 * 365.25 - t_days)/(P/24));

save("ScienceOrbit.mat", "a_final", "e_final", "h_final", "i_final", "RAAN_final", "omega_final", "t_days")
%% === Helper function: run_aerobrake_nominal ===
function [m_TPS_vec, t_TPS_cm_vec, t_days, num_passes, Y_out] = run_aerobrake_nominal(atmos, zonal, merid, temps, beta, alt_entry, mu_N, r_N, r_N_equ)
    % Initial orbit setup
    alt_cap = 10000;                      % km above surface
    rp_cap = alt_cap + r_N;
    vp_cap = sqrt(2 * mu_N/rp_cap);       % km/s
    i_cap = 65; RAAN_cap = 24.8; omega_cap = 99.8;
    [rp_cap_vec, vp_cap_vec] = rv_parabolic(rp_cap, vp_cap, i_cap, RAAN_cap, omega_cap);

    dv_init = 0.22;                               % km/s
    vp_init = vp_cap - dv_init;
    rp_init_vec = rp_cap_vec;
    vp_init_vec = vp_cap_vec * (1 - dv_init/vp_cap);
    [a_init, e_init, h_init, i_init, RAAN_init, omega_init, ~] = oe_from_rv(rp_init_vec, vp_init_vec, mu_N);
    [ra_init_vec, va_init_vec] = rv_from_oe(a_init, e_init, i_init, RAAN_init, omega_init, 180, mu_N);
    ra_init = norm(ra_init_vec);
    va_init = norm(va_init_vec);

    % Use the fixed entry altitude
    alt_brake = alt_entry;
    rp_brake = r_N + alt_brake;
    ra_brake = ra_init;
    e_brake = (ra_brake - rp_brake)./(ra_brake + rp_brake);
    a_brake = 0.5 * (ra_brake + rp_brake);
    h_brake = sqrt(mu_N * a_brake .* (1 - e_brake.^2));
    va_pre_brake = h_brake/ra_brake;
    va_post_brake = va_init;
    va_post_brake - va_pre_brake
    i_brake = i_init; RAAN_brake = RAAN_init; omega_brake = omega_init;
    [ra_brake_init, va_brake_init] = rv_from_oe(a_brake, e_brake, i_brake, RAAN_brake, omega_brake, 180, mu_N);
    Y0 = [ra_brake_init;va_brake_init];

    % Integrate until e < 0.45 at apogee
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

        % if i == N
        %     fileID = fopen('aerobraking_data.txt', 'w');
        % 
        %     % Write data row by row
        %     for k = 1:length(t_pass)
        %         fprintf(fileID, '%.6e\t%.6e\t%.6e\n', t_pass(k), v_pass(k), alt_pass(k));
        %     end
        % 
        %     fclose(fileID);
        % end

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
    cond2 = abs(mod(theta,360) - 180) - 5;
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

function [t] = find_time_neptune(R, h, e, a, mu, altitude)

    r_need = R + altitude;  
    % theta_need = acosd((h^2 / (mu*r_need) - 1)/e)
    % 
    % E_need = 2 * atan( sqrt((1-e)/(1+e)) * tand(theta_need/2) );

    % Clamp cos(theta) to avoid NaNs due to floating point precision
    cos_theta = (h^2 / (mu * r_need) - 1) / e;
    cos_theta = max(-1, min(1, cos_theta));
    theta_rad = acos(cos_theta);  % radians

    % Eccentric anomaly (elliptical orbit)
    E_need = 2 * atan( sqrt((1 - e)/(1 + e)) * tan(theta_rad / 2) );

    M_need = E_need - e * sin(E_need);
    if M_need < 0
        M_need = M_need + 2*pi;
    end
    P_need = 2 * pi * sqrt(a^3 / mu);
    T_need = M_need * P_need/(2*pi);
    t = T_need * 2;
    %t = t /(3600);


end 