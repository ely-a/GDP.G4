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

%% === LOAD ATMOSPHERIC DENSITY MODEL ===
density_data = readmatrix("densitymodel.txt");
rhos = density_data(:, 1) * 1e9;      % Convert to kg/km^3
alts = density_data(:, 2) / 1e3;      % Convert to km
atmos = @(alt_km) interp1(alts, rhos, alt_km, 'linear', 0);  % Interpolation function

%% === INITIAL PARABOLIC ORBIT ===
alt_cap = 10000;                      % km above surface
rp_cap = alt_cap + r_N;
e_cap = 1;
vp_cap = sqrt(2 * mu_N/rp_cap);       % km/s
i_cap = 10; RAAN_cap = 40; omega_cap = 50;
[rp_cap_vec, vp_cap_vec] = rv_parabolic(rp_cap, vp_cap, i_cap, RAAN_cap, omega_cap);

%% === LOAD SCIENCE ORBIT ===
load("ScienceOrbit.mat")

%% === APPLY INITIAL Δv AT PERIAPSIS (TO CREATE ELLIPTICAL ORBIT) ===
dv_init = 0.22;                               % km/s
vp_init = vp_cap - dv_init;
rp_init_vec = rp_cap_vec;
vp_init_vec = vp_cap_vec * (1 - dv_init/vp_cap);
[a_init, e_init, h_init, i_init, RAAN_init, omega_init, ~] = oe_from_rv(rp_init_vec, vp_init_vec, mu_N);
[ra_init_vec, va_init_vec] = rv_from_oe(a_init, e_init, i_init, RAAN_init, omega_init, 180, mu_N);
ra_init = norm(ra_init_vec);
va_init = norm(va_init_vec);

% %% Vary altitude to compare number of passes and time taken
% alt_range = 400:1:550; % km
% num_passes = zeros(size(alt_range));
% time_taken = zeros(size(alt_range));
% max_m_TPS = zeros(size(alt_range));
% max_t_TPS_cm = zeros(size(alt_range));

% for j = 1:length(alt_range)
%     fprintf('Examining entry altitude: %d km\n', alt_range(j));
%     alt_brake = alt_range(j);
%     rp_brake = r_N + alt_brake;

%     e_brake = (ra_init - rp_brake)./(ra_init + rp_brake);
%     a_brake = 0.5 * (ra_init + rp_brake);
%     h_brake = sqrt(mu_N * a_brake .* (1 - e_brake.^2));
%     va_pre_brake = h_brake/ra_init;
%     va_post_brake = va_init;
%     dv_brake = va_pre_brake - va_post_brake;

%     [ra_brake_init, va_brake_init] = rv_from_oe(a_brake, e_brake, i_init, RAAN_init, omega_init, 180, mu_N);
%     Y0_loop = [ra_brake_init;va_brake_init];

%     opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
%                   'Events', @(t, Y) stop_when_reached_science(t, Y, mu_N));

%     [t_out_loop, Y_out_loop] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, ["drag"]), ...
%                                 [0, 0.5 * 365.25 * 86400], Y0_loop, opts);

%     % Count passes
%     r_vecs_loop = Y_out_loop(:, 1:3);
%     altitudes_loop = vecnorm(r_vecs_loop, 2, 2) - r_N;
%     in_atmos_loop = altitudes_loop <= 4000;
%     entry_idx_loop = find(diff(in_atmos_loop) == 1);
%     exit_idx_loop = find(diff(in_atmos_loop) == -1) + 1;
%     num_passes(j) = min(length(entry_idx_loop), length(exit_idx_loop));
%     time_taken(j) = t_out_loop(end)/(24*3600); % convert to days

%     % Initialize max values for this altitude
%     max_m = 0;
%     max_t = 0;

%     for i = 1:num_passes(j)
%         idx_entry = entry_idx_loop(i)+1; % first entry within atmosphere
%         idx_exit = exit_idx_loop(i)-1;   % last point within atmosphere

%         % Extract profiles within the atmosphere for this pass
%         t_pass = t_out_loop(idx_entry:idx_exit) - t_out_loop(idx_entry); % seconds, start from 0
%         alt_pass = altitudes_loop(idx_entry:idx_exit);                   % km
%         v_pass = vecnorm(Y_out_loop(idx_entry:idx_exit,4:6),2,2);        % km/s

%         % Call HeatShieldMass
%         [m_TPS, t_TPS_cm] = HeatShieldMass(alt_pass, v_pass, t_pass);

%         % Update max if this pass is greater
%         if m_TPS > max_m
%             max_m = m_TPS;
%         end
%         if t_TPS_cm > max_t
%             max_t = t_TPS_cm;
%         end
%     end

%     max_m_TPS(j) = max_m;
%     max_t_TPS_cm(j) = max_t;
% end

% figure;
% subplot(2,1,1)
% plot(alt_range, num_passes, '-o', 'LineWidth', 2)
% xlabel('Aerobrake Entry Altitude (km)')
% ylabel('Number of Passes')
% title('Number of Passes vs Entry Altitude')
% grid on

% subplot(2,1,2)
% plot(alt_range, time_taken, '-s', 'LineWidth', 2)
% xlabel('Aerobrake Entry Altitude (km)')
% ylabel('Time to Reach Target Orbit (days)')
% title('Time vs Entry Altitude')
% grid on

% figure;
% subplot(2,1,1)
% plot(alt_range, max_m_TPS, '-o', 'LineWidth', 2)
% xlabel('Aerobrake Entry Altitude (km)')
% ylabel('Max Heat Shield Mass (kg)')
% title('Max Heat Shield Mass vs Entry Altitude')
% grid on

% subplot(2,1,2)
% plot(alt_range, max_t_TPS_cm, '-o', 'LineWidth', 2)
% xlabel('Aerobrake Entry Altitude (km)')
% ylabel('Max Heat Shield Thickness (cm)')
% title('Max Heat Shield Thickness vs Entry Altitude')
% grid on

%% === BRAKING MANEUVER AT APOAPSIS TO DROP PERIAPSIS INTO ATMOSPHERE ===
% --- Set desired periapsis altitude ---
alt_brake = 460;
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

% --- Loop for repeated aerobrake passes with periapsis correction ---
max_passes = 3;
Y_hist = [];
t_hist = [];
m_TPS_vec = [];
t_TPS_cm_vec = [];
delta_vs = [];
eccs = [];
incs = [];
RAANs = [];
omegas = [];
pass_count = 0;

while true
    % Propagate until apogee (theta ~ 180)
    opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
        'Events', @(t, Y) apogee_event(t, Y, mu_N));
    [t_out, Y_out] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, ["drag"]), ...
        [0, 10 * 86400], Y0, opts);

    % Store history
    Y_hist = [Y_hist; Y_out];
    if isempty(t_hist)
        t_hist = t_out;
    else
        t_hist = [t_hist; t_out + t_hist(end)];
    end

    % Get orbital elements at apogee
    r_ap = Y_out(end,1:3)';
    v_ap = Y_out(end,4:6)';
    [a, e, ~, inc, RAAN, omega, theta] = oe_from_rv(r_ap, v_ap, mu_N);

    % Store elements
    eccs(end+1) = e;
    incs(end+1) = inc;
    RAANs(end+1) = RAAN;
    omegas(end+1) = omega;

    % Calculate required velocity at apogee for constant periapsis
    [~, v_ap_req] = rv_from_oe(a_brake, e_brake, inc, RAAN, omega, 180, mu_N);
    dv_vec = v_ap_req - v_ap;
    delta_vs(end+1) = norm(dv_vec);

    % Apply impulse at apogee
    v_ap_new = v_ap + dv_vec;
    Y0 = [r_ap; v_ap_new];

    % Check stop condition (e < 0.42 at apogee)
    if e < 0.42 || pass_count >= max_passes
        break;
    end

    % --- Heat shield calculation for this pass ---
    % Find indices for this pass's in-atmosphere segment
    r_vecs = Y_out(:, 1:3);
    altitudes = vecnorm(r_vecs, 2, 2) - r_N;
    in_atmos = altitudes <= 4000;
    entry_idx = find(diff(in_atmos) == 1) + 1;
    exit_idx = find(diff(in_atmos) == -1);

    if ~isempty(entry_idx) && ~isempty(exit_idx)
        idx_entry = entry_idx(1);
        idx_exit = exit_idx(1);
        t_pass = t_out(idx_entry:idx_exit) - t_out(idx_entry);
        alt_pass = altitudes(idx_entry:idx_exit);
        v_pass = vecnorm(Y_out(idx_entry:idx_exit,4:6),2,2);
        [m_TPS, t_TPS_cm] = HeatShieldMass(alt_pass, v_pass, t_pass);
        m_TPS_vec(end+1) = m_TPS;
        t_TPS_cm_vec(end+1) = t_TPS_cm;
    else
        m_TPS_vec(end+1) = NaN;
        t_TPS_cm_vec(end+1) = NaN;
    end

    pass_count = pass_count + 1;
end

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
eccs = zeros(1, N); incs = zeros(1, N); RAANs = zeros(1, N); omegas = zeros(1, N);
m_TPS_vec = zeros(1, N); t_TPS_cm_vec = zeros(1, N);

for i = 1:N
    [as, eccs(i), hs, incs(i), RAANs(i), omegas(i), ~] = oe_from_rv(Y_out(exit_idx(i), 1:3)', Y_out(exit_idx(i), 4:6)', mu_N);
    theta_exit = acosd(((hs^2 / (mu_N*max_atm_r)) - 1)/eccs(i));
    [~, v_exit] = rv_from_oe(as, eccs(i), incs(i), RAANs(i), omegas(i), theta_exit, mu_N);
    v_exit = norm(v_exit);
    
    [a_entry, ecc_entry, h_entry, inc_entry, RAAN_entry, omega_entry, ~] = oe_from_rv(Y_out(entry_idx(i), 1:3)', Y_out(entry_idx(i), 4:6)', mu_N);
    theta_entry = acosd(((h_entry^2 / (mu_N*max_atm_r)) - 1)/ecc_entry);
    theta_entry = 360 - theta_entry;
    [~, v_entry] = rv_from_oe(a_entry, ecc_entry, inc_entry, RAAN_entry, omega_entry, theta_entry, mu_N);
    v_entry = norm(v_entry);    

    delta_vs(i) = v_exit - v_entry;

    idx_entry = entry_idx(i);
    idx_exit = exit_idx(i);

    idx_entry = entry_idx(i)+1; % first entry within atmosphere
    idx_exit = exit_idx(i)-1; % last point within atmosphere

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
plot(pass_numbers, delta_vs * 1000, '-o', 'LineWidth', 2);
xlabel('Aerobrake Pass'); ylabel('Δv (m/s)');
title('Δv Lost During Each Aerobrake Pass'); grid on;

%% === Plot evolution of orbital elements after each pass ===
figure;
subplot(2,2,1);
plot(pass_numbers, eccs, '-o', 'LineWidth', 2);
xlabel('Pass'); ylabel('Eccentricity'); title('Eccentricity vs Pass'); grid on;

subplot(2,2,2);
plot(pass_numbers, incs, '-o', 'LineWidth', 2);
xlabel('Pass'); ylabel('Inclination (deg)'); title('Inclination vs Pass'); grid on;

subplot(2,2,3);
plot(pass_numbers, RAANs, '-o', 'LineWidth', 2);
xlabel('Pass'); ylabel('RAAN (deg)'); title('RAAN vs Pass'); grid on;

subplot(2,2,4);
plot(pass_numbers, omegas, '-o', 'LineWidth', 2);
xlabel('Pass'); ylabel('\omega (deg)'); title('Argument of Periapsis vs Pass'); grid on;

%% === PLOT PERIAPSIS ALTITUDE PER PASS ===
is_periapsis = islocalmin(altitudes);
periapsis_altitudes = altitudes(is_periapsis);
pass_numbers = 1:length(periapsis_altitudes);
figure;
plot(pass_numbers, periapsis_altitudes, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Aerobrake Pass Number'); ylabel('Periapsis Altitude (km)');
title('Periapsis Altitude vs Aerobrake Pass'); grid on;

%% === Plot heat shield info vs pass number ===
figure;
subplot(2,1,1)
plot(pass_numbers, m_TPS_vec, '-o', 'LineWidth', 2);
xlabel('Aerobrake Pass');
ylabel('Heat Shield Mass (kg)');
title('Heat Shield Mass per Pass');
grid on;

subplot(2,1,2)
plot(pass_numbers, t_TPS_cm_vec, '-o', 'LineWidth', 2);
xlabel('Aerobrake Pass');
ylabel('Heat Shield Thickness (cm)');
title('Heat Shield Thickness per Pass');
grid on;

%% === EVENT FUNCTION: STOP AT APOGEE WITH E < 0.42 ===
function [value, isterminal, direction] = stop_when_reached_science(~, Y, mu)
    r = Y(1:3); v = Y(4:6);
    [~, e, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu);

    % Condition 1: e < 0.42 and at apogee
    cond1 = e - 0.42;
    cond2 = abs(mod(theta,360) - 180) - 2;
    event1 = max([cond1, cond2]); % triggers when both < 0

    % Condition 2: collision with Neptune (r <= r_N)
    r_N = 24622; % km (should match your main code)
    r_mag = norm(r);
    event2 = r_mag - r_N; % triggers when <= 0

    value = [event1, event2];
    isterminal = [1, 1];
    direction = [-1, -1];
end

%% --- Event function for apogee ---
function [value, isterminal, direction] = apogee_event(~, Y, mu)
    r = Y(1:3); v = Y(4:6);
    [~, ~, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu);
    value = abs(mod(theta,360) - 180) - 2; % within 2 deg of apogee
    isterminal = 1;
    direction = 0;
end

%% === DRAG ACCELERATION FUNCTION ===
function a_drag = drag_acc(r, v, R, atmos)
    beta = 1800e6/(1.15 * pi * 2.5^2);
    omega_N = [0; 0; 2 * pi/(16.11 * 3600)];
    v_rel = v - cross(omega_N, r);
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
function dY = eom_perturbed(~, Y, mu, R, R_equ, atmos, perturbations)
    r = Y(1:3); v = Y(4:6);
    r_norm = norm(r);
    a_total = -mu * r / r_norm^3;

    for k = 1:length(perturbations)
        switch perturbations{k}
            case "J2"
                a_total = a_total + J2_acc(r, mu, R_equ);
            case "drag"
                if r_norm - R <= 4000
                    a_total = a_total + drag_acc(r, v, R, atmos);
            end
        end
    end

    dY = [v; a_total];
end
