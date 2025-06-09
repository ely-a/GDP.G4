clear
clc
close all

%% === GRAPHING SETTINGS ===
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontSize',16)
set(groot,'defaulttextfontsize',16)
set(groot,'defaultLineMarkerSize',12)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

%% === CONSTANTS FOR NEPTUNE ===
mu_N = 6.8351e6;         % km^3/s^2
m_N = 102.409e24;        % kg
r_N = 24622;             % km (mean radius)
r_N_equ = 24764;         % km (equatorial radius)

%% === LOAD ATMOSPHERIC DENSITY AND WIND MODELS ===
load("gram_profiles.mat") % provides rho_mean, ew_mean, ns_mean, alt_steps
rhos = rho_mean * 1e9;
alts = alt_steps';
atmos = @(alt_km) interp1(alts, rhos, alt_km, 'linear', 0);

winds_EW = ew_mean/1000; % km/s
winds_NS = ns_mean/1000; % km/s
zonal_model = @(alt_km) interp1(alts, winds_EW, alt_km, 'linear', 0);
meridional_model = @(alt_km) interp1(alts, winds_NS, alt_km, 'linear', 0);

%% === INITIAL PARABOLIC ORBIT ===
alt_cap = 10000;                      % km above surface
rp_cap = alt_cap + r_N;
e_cap = 1;
vp_cap = sqrt(2 * mu_N/rp_cap);       % km/s
i_cap = 23.2; RAAN_cap = 10; omega_cap = 115;
[rp_cap_vec, vp_cap_vec] = rv_parabolic(rp_cap, vp_cap, i_cap, RAAN_cap, omega_cap);

%% === APPLY INITIAL Δv AT PERIAPSIS (TO CREATE ELLIPTICAL ORBIT) ===
dv_init = 0.22;                               % km/s
vp_init = vp_cap - dv_init;
rp_init_vec = rp_cap_vec;
vp_init_vec = vp_cap_vec * (1 - dv_init/vp_cap);
[a_init, e_init, h_init, i_init, RAAN_init, omega_init, ~] = oe_from_rv(rp_init_vec, vp_init_vec, mu_N);
[ra_init_vec, va_init_vec] = rv_from_oe(a_init, e_init, i_init, RAAN_init, omega_init, 180, mu_N);
ra_init = norm(ra_init_vec);
va_init = norm(va_init_vec);

% %% === ANIMATE ORBIT FROM PERIAPSIS TO APOAPASIS OF FIRST ELLIPTICAL ORBIT ===
% % Initial state at periapsis after capture
% Y0_peri = [rp_init_vec; vp_init_vec];
% 
% % Event function to stop at apogee (theta = 180 deg)
% opts_apogee = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
%     'Events', @(t, Y) apogee_event(t, Y, mu_N));
% 
% [t1, Y1] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, zonal_model, meridional_model, ["drag"]), ...
%     [0, 1e7], Y0_peri, opts_apogee);
% 
% r_vecs1 = Y1(:,1:3); % First segment: periapsis to apogee
% 
% % Vary altitude to compare number of passes and time taken
% alt_range = 550:5:550; % km
% num_passes = zeros(size(alt_range));
% time_taken = zeros(size(alt_range));
% max_m_TPS = zeros(size(alt_range));
% max_t_TPS_cm = zeros(size(alt_range));
% 
% for j = 1:length(alt_range)
%     fprintf('Examining entry altitude: %d km\n', alt_range(j));
%     alt_brake = alt_range(j);
%     rp_brake = r_N + alt_brake;
% 
%     e_brake = (ra_init - rp_brake)./(ra_init + rp_brake);
%     a_brake = 0.5 * (ra_init + rp_brake);
%     h_brake = sqrt(mu_N * a_brake .* (1 - e_brake.^2));
%     va_pre_brake = h_brake/ra_init;
%     va_post_brake = va_init;
%     dv_brake = va_pre_brake - va_post_brake;
% 
%     [ra_brake_init, va_brake_init] = rv_from_oe(a_brake, e_brake, i_init, RAAN_init, omega_init, 180, mu_N);
%     Y0_loop = [ra_brake_init;va_brake_init];
% 
%     opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
%                   'Events', @(t, Y) stop_when_reached_science(t, Y, mu_N));
% 
%     [t_out_loop, Y_out_loop] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, zonal_model, meridional_model, ["drag"]), ...
%                                 [0, 1 * 365.25 * 86400], Y0_loop, opts);
% 
%     % Count passes
%     r_vecs_loop = Y_out_loop(:, 1:3);
%     altitudes_loop = vecnorm(r_vecs_loop, 2, 2) - r_N;
%     in_atmos_loop = altitudes_loop <= 4000;
%     entry_idx_loop = find(diff(in_atmos_loop) == 1);
%     exit_idx_loop = find(diff(in_atmos_loop) == -1) + 1;
%     num_passes(j) = min(length(entry_idx_loop), length(exit_idx_loop));
%     time_taken(j) = t_out_loop(end)/(24*3600); % convert to days
% 
%     % Initialize max values for this altitude
%     max_m = NaN;
%     max_t = NaN;
% 
%     if num_passes(j) > 0
%         max_m = 0;
%         max_t = 0;
%         for i = 1:num_passes(j)
%             idx_entry = entry_idx_loop(i)+1; % first entry within atmosphere
%             idx_exit = exit_idx_loop(i)-1;   % last point within atmosphere
% 
%             % Check for valid indices
%             if idx_exit <= idx_entry || idx_entry < 1 || idx_exit > length(t_out_loop)
%                 continue
%             end
% 
%             % Extract profiles within the atmosphere for this pass
%             t_pass = t_out_loop(idx_entry:idx_exit) - t_out_loop(idx_entry); % seconds, start from 0
%             alt_pass = altitudes_loop(idx_entry:idx_exit);                   % km
%             v_pass = vecnorm(Y_out_loop(idx_entry:idx_exit,4:6),2,2);        % km/s
% 
%             % Call HeatShieldMass
%             [m_TPS, t_TPS_cm] = HeatShieldMass(alt_pass, v_pass, t_pass);
% 
%             % Update max if this pass is greater
%             if m_TPS > max_m || isnan(max_m)
%                 max_m = m_TPS;
%             end
%             if t_TPS_cm > max_t || isnan(max_t)
%                 max_t = t_TPS_cm;
%             end
%         end
%     end
% 
%     max_m_TPS(j) = max_m;
%     max_t_TPS_cm(j) = max_t;
% end
% 
% % Altitude variation plots 
% % figure;
% % subplot(2,1,1)
% % plot(alt_range, num_passes, '-', 'LineWidth', 2)
% % xlabel('Aerobrake Aiming Periapsis Altitude (km)')
% % ylabel('Number of Passes')
% % title('Number of Passes vs Aiming Periapsis Altitude')
% % grid on
% 
% figure
% subplot(2,1,1)
% plot(alt_range, time_taken, '-', 'LineWidth', 2)
% hold on
% time_constraint = 0.6 * 365.25; % days
% h_y = yline(time_constraint, 'r--', 'LineWidth', 2, 'Label', '0.6 yr', 'FontSize', 15);
% idx_time = find(time_taken < time_constraint, 1, 'last');
% if ~isempty(idx_time)
%     h_x = xline(alt_range(idx_time), 'r--', 'LineWidth', 2, 'Label', sprintf('Alt = %d km', alt_range(idx_time)),'FontSize', 15, 'LabelVerticalAlignment', 'bottom');
% end
% hold off
% xlabel('Aerobrake Aiming Periapsis Altitude (km)')
% ylabel('Time to Reach Target Orbit (days)')
% title('Time vs Aiming Periapsis Altitude')
% grid on
% 
% mass_constraint = 325; % kg
% subplot(2,1,2)
% hold on
% h_y = yline(mass_constraint, 'r--', 'LineWidth', 2, 'Label', '325 kg', 'FontSize',15);
% idx_mass = find(max_m_TPS < mass_constraint, 1, 'first');
% if ~isempty(idx_mass)
%     h_x = xline(alt_range(idx_mass), 'r--', 'LineWidth', 2, 'Label', sprintf('Alt = %d km', alt_range(idx_mass)), 'FontSize', 15, 'LabelVerticalAlignment','top');
% end
% plot(alt_range, max_m_TPS, '-', 'LineWidth', 2)
% hold off
% xlabel('Aerobrake Aiming Periapsis Altitude (km)')
% ylabel('Max Heat Shield Mass (kg)')
% title('Max Heat Shield Mass vs Aiming Periapsis Altitude')
% grid on
% 
% % subplot(2,1,2)
% % plot(alt_range, max_t_TPS_cm, '-', 'LineWidth', 2)
% % xlabel('Aerobrake Aiming Periapsis Altitude (km)')
% % ylabel('Max Heat Shield Thickness (cm)')
% % title('Max Heat Shield Thickness vs Aiming Periapsis Altitude')
% % grid on

%% === BRAKING MANEUVER AT APOAPSIS TO DROP PERIAPSIS INTO ATMOSPHERE ===
alt_brake = 530;
rp_brake = r_N + alt_brake;
ra_brake = ra_init;
e_brake = (ra_brake - rp_brake)./(ra_brake + rp_brake);
a_brake = 0.5 * (ra_brake + rp_brake);
h_brake = sqrt(mu_N * a_brake .* (1 - e_brake.^2));
va_pre_brake = h_brake/ra_brake;
va_post_brake = va_init;
dv_brake = va_pre_brake - va_post_brake;
i_brake = i_init; RAAN_brake = RAAN_init; omega_brake = omega_init;
[ra_brake_init, va_brake_init] = rv_from_oe(a_brake, e_brake, i_brake, RAAN_brake, omega_brake, 180, mu_N);
Y0 = [ra_brake_init;va_brake_init];

%% === ANIMATE TRAJECTORY FROM APOAPSIS TO ATMOSPHERE ENTRY ===
% % Initial state at apogee after braking maneuver
% Y0_apogee = [ra_brake_init; va_brake_init];
% 
% % Event function to stop at atmosphere entry (altitude = 4000 km)
% opts_atm = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
%     'Events', @(t, Y) atmosphere_entry_event(t, Y, r_N));
% 
% [t2, Y2] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, zonal_model, meridional_model, ["drag"]), ...
%     [0, 1e7], Y0_apogee, opts_atm);
% 
% r_vecs2 = Y2(:,1:3); % Apogee to just before atmosphere
% 
% % Remove the last point of r_vecs1 and t1 to avoid overlap at apogee
% r_vecs_combined = [r_vecs1(1:end-1, :); r_vecs2];
% t_combined = [t1(1:end-1); t2 + t1(end)];
% 
% % animate_trajectory(r_vecs_combined, t_combined, r_N, 4000, "CaptureToEntry");

%% === PROPAGATE UNDER DRAG UNTIL E < 0.45 AT APOAPSIS ===
opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
              'Events', @(t, Y) stop_when_reached_science(t, Y, mu_N));
[t_out, Y_out] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, zonal_model, meridional_model, ["drag"]), ...
                        [0, 1 * 365.25 * 86400], Y0, opts);

%% === ANIMATE ALL AEROBRAKES ===
% animate_trajectory(Y_out(1:1200, 1:3), t_out(1:1200,1), r_N, 5000, "Aerobrakes")

%% === ANIMATE ONE AEROBRAKE ===
% % Find indices for the first atmospheric pass
% altitudes = vecnorm(Y_out(:,1:3), 2, 2) - r_N;
% in_atmos = altitudes <= 4000;
% entry_idx = find(diff(in_atmos) == 1, 1, 'first') + 1; % first entry within atmosphere
% exit_idx = find(diff(in_atmos) == -1, 1, 'first');     % first exit from atmosphere
% 
% % Safety check in case exit_idx is empty (e.g., s/c never leaves atmosphere)
% if isempty(exit_idx)
%     exit_idx = length(altitudes);
% end
% 
% % Extract position and time vectors for the first aerobrake
% r_vecs_ab = Y_out(entry_idx:exit_idx, 1:3);
% t_vec_ab = t_out(entry_idx:exit_idx);
% 
% % Animate the first aerobrake pass
% % animate_trajectory(r_vecs_ab, t_vec_ab, r_N, 50, 'FirstAerobrake.mp4');

%% === PROCESS OUTPUT ===
r_vecs = Y_out(:, 1:3);
altitudes = vecnorm(r_vecs, 2, 2) - r_N;
in_atmos = altitudes <= 4000;                  % Logical array: inside drag region
entry_idx = find(diff(in_atmos) == 1);         % Just before entry
exit_idx = find(diff(in_atmos) == -1) + 1;     % Just after exit
N = min(length(entry_idx), length(exit_idx));
max_atm_r = r_N + alts(end);

% Compute Δv lost and orbital elements and heat shield per pass
delta_vs = zeros(1, N); 
as=zeros(1,N); eccs = zeros(1, N); incs = zeros(1, N); RAANs = zeros(1, N); omegas = zeros(1, N);
m_TPS_vec = zeros(1, N); t_TPS_cm_vec = zeros(1, N);

for i = 1:N
    [as(i), eccs(i), hs, incs(i), RAANs(i), omegas(i), ~] = oe_from_rv(Y_out(exit_idx(i), 1:3)', Y_out(exit_idx(i), 4:6)', mu_N);
    theta_exit = acosd(((hs^2 / (mu_N*max_atm_r)) - 1)/eccs(i));
    [~, v_exit] = rv_from_oe(as(i), eccs(i), incs(i), RAANs(i), omegas(i), theta_exit, mu_N);
    v_exit = norm(v_exit);
    
    [a_entry, ecc_entry, h_entry, inc_entry, RAAN_entry, omega_entry, ~] = oe_from_rv(Y_out(entry_idx(i), 1:3)', Y_out(entry_idx(i), 4:6)', mu_N);
    theta_entry = acosd(((h_entry^2 / (mu_N*max_atm_r)) - 1)/ecc_entry);
    theta_entry = 360 - theta_entry;
    [~, v_entry] = rv_from_oe(a_entry, ecc_entry, inc_entry, RAAN_entry, omega_entry, theta_entry, mu_N);
    v_entry = norm(v_entry);    

    delta_vs(i) = v_exit - v_entry;

    idx_entry = entry_idx(i)+1; % first entry within atmosphere
    idx_exit = exit_idx(i)-1;   % last point within atmosphere

    % Extract profiles within the atmosphere for this pass
    t_pass = t_out(idx_entry:idx_exit) - t_out(idx_entry); % seconds, start from 0
    alt_pass = altitudes(idx_entry:idx_exit);              % km
    v_pass = vecnorm(Y_out(idx_entry:idx_exit,4:6),2,2);   % km/s

    % Call HeatShieldMass (expects altitude in km, velocity in km/s, time in s)
    [m_TPS, t_TPS_cm] = HeatShieldMass(alt_pass, v_pass, t_pass);

    m_TPS_vec(i) = m_TPS;
    t_TPS_cm_vec(i) = t_TPS_cm;
end

% Display results
t_days = t_out(end)/(24 * 3600);
fprintf('Number of atmosphere entries: %d\n', N);
fprintf('Time taken: %.2f days \n', t_days);

%% === SAVE OE VARIATION ===
a_brakes = [a_brake, as];
e_brakes = [e_brake, eccs];
i_brakes = [i_brake, incs];
RAAN_brakes = [RAAN_brake, RAANs];
omega_brakes = [omega_brake, omegas];

save("Aerobrakes_OE.mat", "a_brakes", "omega_brakes", "RAAN_brakes", "i_brakes", "e_brakes")

%% === PLOT TRAJECTORY IN 3D ===
r_vecs = Y_out(:, 1:3)';
figure;
hold on;
plot3(r_vecs(1,:), r_vecs(2,:), r_vecs(3,:), 'b', 'LineWidth', 1.5);
[Xs, Ys, Zs] = sphere(100);
surf(r_N*Xs, r_N*Ys, r_N*Zs, ...
    'FaceColor', [0.2 0.4 1], 'EdgeAlpha', 0.1, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('3D Trajectory Around Neptune');
axis equal; grid on; view(3);

%% === PLOT Δv PER PASS ===
pass_numbers = 1:N;
figure;
subplot(2,1,1)
plot(pass_numbers, delta_vs * 1000, '-', 'LineWidth', 2);
xlabel('Aerobrake Pass'); ylabel('Δv (m/s)');
title('Δv Lost During Each Aerobrake Pass'); grid on;

%% === PLOT PERIAPSIS ALTITUDE PER PASS ===
is_periapsis = islocalmin(altitudes);
periapsis_altitudes = altitudes(is_periapsis);
pass_numbers = 1:length(periapsis_altitudes);
subplot(2,1,2)
plot(pass_numbers, periapsis_altitudes, '-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Aerobrake Pass Number'); ylabel('Periapsis Altitude (km)');
title('Periapsis Altitude vs Aerobrake Pass'); grid on;

%% === Plot evolution of orbital elements after each pass ===
figure;
subplot(2,2,1);
plot(pass_numbers, eccs, '-', 'LineWidth', 2);
xlabel('Pass'); ylabel('Eccentricity'); title('Eccentricity vs Pass'); grid on;

subplot(2,2,2);
plot(pass_numbers, incs, '-', 'LineWidth', 2);
xlabel('Pass'); ylabel('Inclination (deg)'); title('Inclination vs Pass'); grid on;

subplot(2,2,3);
plot(pass_numbers, RAANs, '-', 'LineWidth', 2);
xlabel('Pass'); ylabel('RAAN (deg)'); title('RAAN vs Pass'); grid on;

subplot(2,2,4);
plot(pass_numbers, omegas, '-', 'LineWidth', 2);
xlabel('Pass'); ylabel('\omega (deg)'); title('Argument of Periapsis vs Pass'); grid on;

%% === Plot heat shield info vs pass number ===
figure;
subplot(2,1,1)
plot(pass_numbers, m_TPS_vec, '-', 'LineWidth', 2);
xlabel('Aerobrake Pass');
ylabel('Heat Shield Mass (kg)');
title('Heat Shield Mass per Pass');
grid on;

subplot(2,1,2)
plot(pass_numbers, t_TPS_cm_vec, '-', 'LineWidth', 2);
xlabel('Aerobrake Pass');
ylabel('Heat Shield Thickness (cm)');
title('Heat Shield Thickness per Pass');
grid on;

%% === SAVE FINAL AEROBRAKE DATA ===
% --- Save time, velocity, and altitude for the final aerobrake pass ---
% Find indices for all atmospheric entries/exits
altitudes = vecnorm(Y_out(:,1:3), 2, 2) - r_N;
in_atmos = altitudes <= 4000;
entry_idx = find(diff(in_atmos) == 1);         % Just before entry
exit_idx = find(diff(in_atmos) == -1) + 1;     % Just after exit

% Use the last pass
if ~isempty(entry_idx) && ~isempty(exit_idx)
    idx_entry = entry_idx(end) + 1; % first point inside atmosphere
    idx_exit = exit_idx(end) - 1;   % last point inside atmosphere

    t_pass = t_out(idx_entry:idx_exit); % time stamps (s)
    v_pass = vecnorm(Y_out(idx_entry:idx_exit,4:6),2,2); % velocity (km/s)
    alt_pass = altitudes(idx_entry:idx_exit); % altitude (km)

    % Combine into one matrix for saving
    data = [t_pass, v_pass, alt_pass];

    % Write to text file
    writematrix(data, 'aerobraking_data.txt', 'Delimiter', 'tab');
    fprintf('Final aerobrake data saved to aerobraking_data.txt\n');
else
    warning('No atmospheric pass found to save.');
end

%% === REACH SCIENCE ORBIT ===
% Assume rp_final is already defined (in km, measured from Neptune's center)
% Get final state after last aerobrake
r_final = Y_out(end,1:3)';
v_final = Y_out(end,4:6)';

% Get current orbital elements
[a_curr, e_curr, h_curr, i_curr, RAAN_curr, omega_curr, theta_curr] = oe_from_rv(r_final, v_final, mu_N);

% Current apogee radius
ra_curr = a_curr * (1 + e_curr);
rp_final = 2.620767634816547e+04;
% Required eccentricity for desired periapsis
e_req = (ra_curr - rp_final) / (ra_curr + rp_final);

% Required semi-major axis for new orbit
a_req = 0.5 * (ra_curr + rp_final);

% Velocity at apogee before maneuver
v_ap_before = norm(v_final);

% Velocity at apogee after maneuver (for new orbit)
v_ap_after = sqrt(mu_N * (2/ra_curr - 1/a_req));

% Required delta-v at apogee
delta_v_apogee = v_ap_after - v_ap_before;

% Apply the final Δv at apogee to the velocity vector
v_final_after_impulse = v_final + (delta_v_apogee / norm(v_final)) * v_final;

% Get the new orbital elements after the impulse
[a_final, e_final, h_final, i_final, RAAN_final, omega_final, theta_final] = oe_from_rv(r_final, v_final_after_impulse, mu_N);

% Calculate the period of the final orbit
T_final = 2 * pi * sqrt(a_final^3 / mu_N); % seconds

fprintf('\nFinal Science Orbit Elements after Δv:\n');
fprintf('  Semi-major axis (a):      %.2f km\n', a_final);
fprintf('  Eccentricity (e):         %.5f\n', e_final);
fprintf('  Inclination (i):          %.2f deg\n', i_final);
fprintf('  RAAN:                     %.2f deg\n', RAAN_final);
fprintf('  Argument of Periapsis:    %.2f deg\n', omega_final);
fprintf('  True Anomaly:             %.2f deg\n', theta_final);
fprintf('  Period:                   %.2f hours\n', T_final / 3600); % convert to hours

% Initial state after impulse
Y0_final = [r_final; v_final_after_impulse];

% Propagate for one period (no drag, no J2)
[t_orb, Y_orb] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, zonal_model, meridional_model, {}), ...
    [0, T_final], Y0_final);

% Extract position vectors
r_orb = Y_orb(:,1:3);

% Plot the orbit
figure;
hold on;
[Xs, Ys, Zs] = sphere(100);
surf(r_N*Xs, r_N*Ys, r_N*Zs, ...
    'FaceColor', [0.2 0.4 1], 'EdgeAlpha', 0.1, ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot3(r_orb(:,1), r_orb(:,2), r_orb(:,3), 'r', 'LineWidth', 2);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('One Orbit of Final Science Trajectory Around Neptune');
axis equal; grid on; view(3);
legend('Neptune','Final Science Orbit');
hold off;

save("ScienceOrbit.mat", "a_final", "e_final", "h_final", "i_final", "RAAN_final", "omega_final", "t_days")

%% === EVENT FUNCTIONS: STOP AT APOGEE WITH E < 0.45 ===
function [value, isterminal, direction] = stop_when_reached_science(~, Y, mu)
    r = Y(1:3); v = Y(4:6);
    [~, e, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu);

    % Condition 1: e < 0.45 and at apogee
    cond1 = e - 0.45;
    cond2 = abs(mod(theta,360) - 180) - 0.1;
    event1 = max([cond1, cond2]); % triggers when both < 0

    % Condition 2: collision with Neptune (r <= r_N)
    r_N = 24622; % km (should match your main code)
    r_mag = norm(r);
    event2 = r_mag - r_N; % triggers when <= 0

    value = [event1, event2];
    isterminal = [1, 1];
    direction = [-1, -1];
end

function [value, isterminal, direction] = apogee_event(~, Y, mu)
    r = Y(1:3); v = Y(4:6);
    [~, ~, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu);
    value = abs(mod(theta,360) - 180) - 0.01; % within 2 deg of apogee
    isterminal = 1;
    direction = 0;
end

function [value, isterminal, direction] = atmosphere_entry_event(~, Y, r_N)
    r = Y(1:3);
    alt = norm(r) - r_N;
    value = alt - 4000; % triggers when altitude = 4000 km
    isterminal = 1;
    direction = -1;
end

%% === WIND INERTIAL VELOCITY FUNCTION ===
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

%% === DRAG ACCELERATION FUNCTION ===
function a_drag = drag_acc(r, v, R, atmos, zonal_model, meridional_model)
    beta = 1800e6/(1.15 * pi * 2.5^2);
    omega_N = [0; 0; 2 * pi/(16.11 * 3600)];
    v_rel = v - cross(omega_N, r) - v_wind_inertial(r, R, zonal_model, meridional_model);
    rho = atmos(norm(r) - R);
    a_drag = -0.5 * rho * norm(v_rel) * v_rel / beta;
end

%% === J2 PERTURBATION ACCELERATIONS ===
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

%% === EQUATIONS OF MOTION INCLUDING DRAG (AND OPTIONAL J2) ===
function dY = eom_perturbed(~, Y, mu, R, R_equ, atmos, zonal_model, meridional_model, perturbations)
    r = Y(1:3); v = Y(4:6);
    r_norm = norm(r);
    a_total = -mu * r / r_norm^3;

    for k = 1:length(perturbations)
        switch perturbations{k}
            case "J2"
                a_total = a_total + J2_acc(r, mu, R_equ);
            case "drag"
                if r_norm - R <= 4000
                    a_total = a_total + drag_acc(r, v, R, atmos, zonal_model, meridional_model);
                end
        end
    end

    dY = [v; a_total];
end

%% === ANIMATE TRAJECTORY ===
function animate_trajectory(r_vecs, t_vec, r_N, resolution, filename)
% animate_trajectory Animate a 3D trajectory around Neptune and save as MP4
%   r_vecs: Nx3 array of position vectors (km)
%   t_vec:  Nx1 array of time stamps (s)
%   r_N:    Neptune radius (km)
%   filename: name of the mp4 file (string), e.g. 'my_animation.mp4'

    % Interpolate at fixed time steps
    t_uniform = (t_vec(1):resolution:t_vec(end))';
    r_vecs_uniform = interp1(t_vec, r_vecs, t_uniform, 'spline');

    % Ensure output folder exists
    out_folder = 'Aerobraking Figures';
    if ~exist(out_folder, 'dir')
        mkdir(out_folder);
    end
    out_path = fullfile(out_folder, filename);

    % Set up video writer
    v = VideoWriter(out_path, 'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 60;
    open(v);

    figure;
    set(gcf, 'Position', [100, 100, 1600, 900]); % Set figure size to 1600x900
    hold on;
    [Xs, Ys, Zs] = sphere(100);
    surf(r_N*Xs, r_N*Ys, r_N*Zs, ...
        'FaceColor', [0.2 0.4 1], 'EdgeAlpha', 0.1, ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
    title('Trajectory Animation Around Neptune');

    % Fixed frame
    max_radius = max([r_N; vecnorm(r_vecs_uniform,2,2)]);
    lim = max_radius * 1.1;
    xlim([-lim, lim]);
    ylim([-lim, lim]);
    zlim([-lim, lim]);
    axis equal; grid on; view(3);

    h_traj = plot3(r_vecs_uniform(1,1), r_vecs_uniform(1,2), r_vecs_uniform(1,3), 'b-', 'LineWidth', 2);
    % Use a star marker for the spacecraft, smaller size
    h_sc = plot3(r_vecs_uniform(1,1), r_vecs_uniform(1,2), r_vecs_uniform(1,3), '*', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

    % Draw legend once and keep it for all frames
    legend({'Neptune','Trajectory','Spacecraft'},'Location','best','AutoUpdate','off');

    % Animation loop
    for k = 2:length(t_uniform)
        set(h_traj, 'XData', r_vecs_uniform(1:k,1), 'YData', r_vecs_uniform(1:k,2), 'ZData', r_vecs_uniform(1:k,3));
        set(h_sc, 'XData', r_vecs_uniform(k,1), 'YData', r_vecs_uniform(k,2), 'ZData', r_vecs_uniform(k,3));
        drawnow;
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    hold off;
    close(v);
    disp(['Animation saved to ', out_path]);
end