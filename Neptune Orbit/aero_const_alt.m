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
load("gram_profiles.mat")
rhos = rho_mean * 1e9;
alts = alt_steps';
atmos = @(alt_km) interp1(alts, rhos, alt_km, 'linear', 0);  % Interpolation function

%% === LOAD WIND MODELS ===

winds_EW = ew_mean/1000;
winds_NS = ns_mean/1000;

zonal_model =  @(alt_km) interp1(alts, winds_EW, alt_km, 'linear', 0);
meridonial_model =  @(alt_km) interp1(alts, winds_NS,alt_km, 'linear', 0);

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

%% === Loop to iterate passes === 

states = [];
pass = 0;
delta_v_brake_total = 0;

while true
    pass = pass + 1;

    opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
              'Events', @(t, Y) apogee_event(t, Y, mu_N));
    [t_out, Y_out] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, zonal_model, meridonial_model, ["drag"]), ...
                        [0, 25 * 86400], Y0, opts);

    % --- Apply impulse at apogee to restore periapsis radius (rp_brake) ---
    r_ap = Y_out(end, 1:3)';
    v_ap = Y_out(end, 4:6)';

    % Store orbital elements before impulse at apogee
    [a_pass(pass), e_pass(pass), ~, i_pass(pass), RAAN_pass(pass), omega_pass(pass), ~] = oe_from_rv(r_ap, v_ap, mu_N);

    % Get current orbital elements
    [a, e, ~, i, RAAN, omega, theta] = oe_from_rv(r_ap, v_ap, mu_N);
    ra = a * (1 + e);

    % --- Extract atmospheric segment for this pass ---
    altitudes_pass = vecnorm(Y_out(:,1:3),2,2) - r_N;
    in_atmos = altitudes_pass <= 4000;
    idx_entry = find(diff(in_atmos)==1,1,'first')+1;
    idx_exit = find(diff(in_atmos)==-1,1,'first');
    if isempty(idx_entry) || isempty(idx_exit)
        idx_entry = 1;
        idx_exit = length(altitudes_pass);
    end
    alt_profile = altitudes_pass(idx_entry:idx_exit);
    vel_profile = vecnorm(Y_out(idx_entry:idx_exit,4:6),2,2);
    t_profile = t_out(idx_entry:idx_exit);

    % --- Calculate heat shield mass and thickness for this pass ---
    [m_TPS_vec(pass), t_TPS_cm_vec(pass)] = HeatShieldMass(alt_profile, vel_profile, t_profile);

    if pass == 1
        times = t_out;
    else
        times = [times; t_out + times(end)];
    end
    states = [states; Y_out];

    % Optional: break condition to avoid infinite loop
    if times(end) > 1 * 365.25 * 86400 || e < 0.45
        break
    end

    [r_ap, v_ap] = rv_from_oe(a, e, i, RAAN, omega, 180, mu_N);

    % Compute required semi-major axis and eccentricity for target periapsis
    a_new = 0.5 * (ra + rp_brake);
    e_new = (ra - rp_brake) / (ra + rp_brake);

    % Compute required velocity at apogee for new orbit
    v_ap_required = sqrt(mu_N * (2/ra - 1/a_new));
    v_ap_current = norm(v_ap);

    delta_v_brake = abs(v_ap_required - v_ap_current);

    % Apply the impulse in the direction of current velocity
    v_ap_new = v_ap + (delta_v_brake / v_ap_current) * v_ap;

    % Use [r_ap; v_ap_new] as the initial state for the next propagation
    Y0 = [r_ap; v_ap_new];

    if norm(r_ap) <= r_N
        warning('Spacecraft has impacted Neptune. Exiting loop.');
        break
    end

    delta_v_brake_total = delta_v_brake_total + delta_v_brake;
end

% Display results
t_days = times(end)/(24 * 3600);
fprintf('Number of atmosphere entries: %d\n', pass);
fprintf('Time taken: %.2f days \n', t_days);

%% === PLOT TRAJECTORY ===
r_traj = states(:,1:3);

figure;
hold on;
[Xs, Ys, Zs] = sphere(100);
surf(r_N*Xs, r_N*Ys, r_N*Zs, ...
    'FaceColor', [0.2 0.4 1], 'EdgeAlpha', 0.1, ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot3(r_traj(:,1), r_traj(:,2), r_traj(:,3), 'r', 'LineWidth', 2);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Aerobraking Trajectory with Constant Periapsis Altitude');
axis equal; grid on; view(3);
legend('Neptune','Trajectory');
hold off;

%% === PLOT ALTITUDE PER PASS ===
% Calculate altitude at each state
altitudes = vecnorm(r_traj, 2, 2) - r_N;

% Find local minima (periapsis altitudes for each pass)
is_periapsis = islocalmin(altitudes);
periapsis_altitudes = altitudes(is_periapsis);
pass_numbers = 1:length(periapsis_altitudes);

figure;
plot(pass_numbers, periapsis_altitudes, '-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Pass Number');
ylabel('Periapsis Altitude (km)');
title('Minimum Altitude (Periapsis) of Each Pass');
grid on;

%% === Plot orbital elements at apogee before each impulse ===
passes = 1:pass;

figure;
subplot(3,2,1);
plot(passes, a_pass, '-');
xlabel('Pass'); ylabel('a (km)'); title('Semi-major Axis');

subplot(3,2,2);
plot(passes, e_pass, '-');
xlabel('Pass'); ylabel('Eccentricity'); title('Eccentricity');

subplot(3,2,3);
plot(passes, i_pass, '-');
xlabel('Pass'); ylabel('Inclination (deg)'); title('Inclination');

subplot(3,2,4);
plot(passes, RAAN_pass, '-');
xlabel('Pass'); ylabel('RAAN (deg)'); title('RAAN');

subplot(3,2,5);
plot(passes, omega_pass, '-');
xlabel('Pass'); ylabel('\omega (deg)'); title('Argument of Periapsis');

sgtitle('Orbital Elements Before Each Apogee Impulse');

%% === PLOT HEAT SHIELD MASS AND THICKNESS PER PASS ===
figure;
subplot(2,1,1)
plot(1:pass, m_TPS_vec, '-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Pass Number');
ylabel('Heat Shield Mass (kg)');
title('Heat Shield Mass per Pass');
grid on;

subplot(2,1,2)
plot(1:pass, t_TPS_cm_vec, '-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Pass Number');
ylabel('Heat Shield Thickness (cm)');
title('Heat Shield Thickness per Pass');
grid on;

%% === REACH SCEIENCE ORBIT ===

r_final= states(end, 1:3);
v_final = states(end, 4:6);

% Get current orbital elements
[a_curr, e_curr, h_curr, i_curr, RAAN_curr, omega_curr, theta_curr] = oe_from_rv(r_final, v_final, mu_N);

% Current apogee radius
ra_curr = a_curr * (1 + e_curr);
rp_final = 2.620767634816547e+04;
% rp_final = r_N + 4000;
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

%% === GET TOTAL DELTA V ===

dv_total = abs(dv_init) + abs(dv_brake) + abs(delta_v_brake_total) + abs(delta_v_apogee);

fprintf('Total delta v: %.4f km/s \n', dv_total);


%% === APOGEE EVENT ACHEIVED === 

function [value, isterminal, direction] = apogee_event(~, Y, mu)
    r = Y(1:3); v = Y(4:6);
    [~, ~, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu);

    % Event triggers at apogee (theta = 180 deg)
    value = mod(theta,360) - 180 + 0.01; % within 0.5 deg of apogee
    isterminal = 1;
    direction = 1;
end

%% === WIND FUNCTION ===

function v_wind_inertial = v_wind_inertial(r_vec, R, zonal_model, meridional_model)
% Computes the wind vector in inertial frame from local zonal/meridional wind models.
%
% Inputs:
%   r_vec          - 3x1 position vector of spacecraft [km] in inertial frame
%   v_vec          - 3x1 velocity vector of spacecraft [km/s] (not used, included for interface compatibility)
%   r_planet       - Mean radius of Neptune [km]
%   zonal_model    - function handle for zonal wind (altitude [km] → km/s)
%   meridional_model - function handle for meridional wind (altitude [km] → km/s)
%
% Output:
%   v_wind_inertial - 3x1 atmospheric wind vector in inertial frame [km/s]

    % Position components
    x = r_vec(1); y = r_vec(2); z = r_vec(3);
    r_mag = norm(r_vec);

    % Altitude above mean planetary surface
    alt = r_mag - R;

    % Convert position to latitude and longitude (assumes Z is spin axis)
    lat = asind(z / r_mag);                   % degrees
    lon = atan2d(y, x);                       % degrees

    % Convert to radians for unit vectors
    lat_rad = deg2rad(lat);
    lon_rad = deg2rad(lon);

    % East (zonal) unit vector in inertial frame
    e_hat = [-sind(lon_rad); cosd(lon_rad); 0];

    % North (meridional) unit vector in inertial frame
    n_hat = [-sind(lat_rad)*cosd(lon_rad);
             -sind(lat_rad)*sind(lon_rad);
              cosd(lat_rad)];

    % Wind speeds from model [km/s]
    v_zonal = zonal_model(alt);       % positive eastward
    v_meridional = meridional_model(alt); % positive northward

    % Wind vector in inertial frame
    v_wind_inertial = v_zonal * e_hat + v_meridional * n_hat;
end


%% === DRAG ACCELERATION FUNCTION ===
function a_drag = drag_acc(r, v, R, atmos, zonal_model, meridonial_model)
    beta = 1800e6/(1.15 * pi * 2.5^2);
    omega_N = [0; 0; 2 * pi/(16.11 * 3600)];
    v_rel = v - cross(omega_N, r) - v_wind_inertial(r, R, zonal_model, meridonial_model);
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

%% === EQUATIONS OF MOTION ===
function dY = eom_perturbed(~, Y, mu, R, R_equ, atmos, zonal, meridonial, perturbations)
    r = Y(1:3); v = Y(4:6);
    r_norm = norm(r);
    a_total = -mu * r / r_norm^3;

    for k = 1:length(perturbations)
        switch perturbations{k}
            case "J2"
                a_total = a_total + J2_acc(r, mu, R_equ);
            case "drag"
                if r_norm - R <= 4000
                    a_total = a_total + drag_acc(r, v, R, atmos, zonal, meridonial);
                end
        end
    end

    dY = [v; a_total];
end