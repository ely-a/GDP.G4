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

%% Original position

if exist('ScienceOrbit.mat', 'file') == 2
    load("ScienceOrbit.mat")
else
    error('no intitial conditions')
end

% [r0, v0] = rv_from_oe(3e5, 0.9, 90, 90, 250, 0, mu_neptune);
[r0, v0] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, 180, mu_N);
ra_final = a_final * (1 + e_final);

aiming_periapsis = 0:1:0;
rp_deorbit = aiming_periapsis + r_N;
a_deorbits = 0.5 * (ra_final + rp_deorbit);
va_final = norm(v0);
va_deorbits = sqrt(mu_N * (2/ra_final - 1./a_deorbits));
dv_deorbit = va_deorbits - va_final;
alt_entry = 4000 + r_N;

impact_times = NaN(size(aiming_periapsis));  % Use NaN to mark cases with no impact

for idx = 1:length(aiming_periapsis)
    fprintf('Aiming periapsis:  %.1f km \n', aiming_periapsis(idx));
    va_deorbit = v0 + dv_deorbit(idx) * (v0/va_final);
    [a_deorbit, e_deorbit, h_deorbit, i_deorbit, RAAN_deorbit, omega_deorbit, ~] = oe_from_rv(r0, va_deorbit, mu_N);
    theta_entry = 360 - acosd((h_deorbit^2/(mu_N * alt_entry) - 1)/e_deorbit);
    [r_entry, v_entry] = rv_from_oe(a_deorbit, e_deorbit, i_deorbit, RAAN_deorbit, omega_deorbit, theta_entry, mu_N);

    Y0 = [r_entry; v_entry];

    tspan = [0, 6.6 * 3600];         % time steps
    
    % Integrate perturbed motion
    opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
              'Events', @(t,Y) impact_event(t, Y, r_N));
    perturbs = ["drag"];
    
    
    [t_out, Y_out, te, Ye, ie] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, perturbs, ...
                                                 atmos_nom, zonal_nom, merid_nom, beta_nom), ...
                                                 tspan, Y0, opts);

    r_perturbed = Y_out(:,1:3)';
    v_perturbed = Y_out(:,4:6)';
    
    % === Check if impact occurred ===
    if ~isempty(te)
        fprintf('   → Impact detected at t = %.2f s, altitude = %.2f km\n', ...
                te(end), norm(Ye(end,1:3)) - r_N);
        impact_times(idx) = te(end);  % Store impact time
    else
        fprintf('   → No impact: spacecraft remained above surface for full duration.\n');
    end

end

figure;
plot(abs(dv_deorbit), impact_times / 60, 'b-');  % Convert seconds to hours
xlabel('\DeltaV_{deorbit} [km/s]');
ylabel('Impact Time [hours]');
title('Impact Time vs Deorbit \DeltaV');
% axis([])
grid on;

%% Trajectories

dv_last = 0.17;
va_last_vec = v0 - dv_last * (v0/va_final);
ra_last = r0;

Y0 = [ra_last; va_last_vec];

tspan = [0, 24 * 3600];         % time steps

% Integrate perturbed motion
opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
    'Events', @(t,Y) impact_event(t, Y, r_N));
perturbs = ["drag"];

[t_out, Y_out, te, Ye, ie] = ode113(@(t,Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, perturbs, ...
    atmos_nom, zonal_nom, merid_nom, beta_nom), ...
    tspan, Y0, opts);

r_traj = Y_out(:,1:3);
v_traj = Y_out(:,4:6);

figure
plot3(r_traj(:, 1), r_traj(:, 2), r_traj(:, 3))
hold on
[Xs, Ys, Zs] = sphere(100);
surf(r_N*Xs, r_N*Ys, r_N*Zs, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');

%% Event function - collision
function [value, isterminal, direction] = impact_event(t, Y, R)
    r = Y(1:3);
    value = norm(r) - R;     % Distance above Neptune's surface
    isterminal = 1;          % Stop integration
    direction = -1;          % Only trigger when decreasing through zero
end

%% Drag
function a_drag = drag_acc(r, v, R, atmos, zonal_model, meridional_model, beta)
    omega_N = [0; 0; 2 * pi/(16.11 * 3600)];
    v_rel = v - cross(omega_N, r) - v_wind_inertial(r, R, zonal_model, meridional_model);
    rho = atmos(norm(r) - R);
    a_drag = -0.5 * rho * norm(v_rel) * v_rel / beta;
end

%% EOM

function dY = eom_perturbed(t, Y, mu_N, R, R_equ, perturbations, atmos, zonal, merid, beta)
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