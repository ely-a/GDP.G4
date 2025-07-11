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

dv_init = 0.22;
vp_init = vp_cap - dv_init;
rp_init_vec = rp_cap_vec;
vp_init_vec = vp_cap_vec * (1 - dv_init/vp_cap);

[a_init, e_init, h_init, i_init, RAAN_init, omega_init, ~] = oe_from_rv(rp_init_vec, vp_init_vec, mu_N);
[ra_init_vec, va_init_vec] = rv_from_oe(a_init, e_init, i_init, RAAN_init, omega_init, 180, mu_N);

ra_init = norm(ra_init_vec);
va_init = norm(va_init_vec);

%% Aerobrake propogation for densty uncertatinities

max_atm_alt = alts(end);
max_atm_r = r_N + max_atm_alt;

ra_brake = ra_init;
va_post_brake = va_init;
i_brake = i_init;
RAAN_brake = RAAN_init;
omega_brake = omega_init;

aerobrakes = 1;

for no_brake = 1:aerobrakes

    alt_brake = 450;
    rp_brake = r_N + alt_brake;

    e_brake = (ra_brake - rp_brake)./(ra_brake + rp_brake);
    a_brake = 0.5 * (ra_brake + rp_brake);
    h_brake = sqrt(mu_N * a_brake .* (1 - e_brake.^2));
    va_pre_brake = h_brake/ra_brake;

    dv_brake = va_pre_brake - va_post_brake;
    
    theta_entry = acosd(((h_brake^2 / (mu_N*max_atm_r)) - 1)/e_brake);
    theta_entry = 360 - theta_entry;
    [r_entry, v_entry] = rv_from_oe(a_brake, e_brake, i_brake, RAAN_brake, omega_brake, theta_entry, mu_N);
    Y0 = [r_entry; v_entry];
    opts = odeset('RelTol',1e-9, 'AbsTol',1e-9, ...
        'Events', @(t, Y) exit_atmosphere_event(t, Y, max_atm_r));
    [t_out, Y_out, te, ye, ie] = ode113(@(t, Y) eom_perturbed(t, Y, mu_N, r_N, r_N_equ, atmos, ["drag"]), ...
        [0, 2e4], Y0, opts);
    
    [a_brake, e_brake, h_brake, i_brake, Omega_brake, omega_brake, ~] = oe_from_rv(Y_out(end, 1:3)', Y_out(end, 4:6)', mu_N);
    ra_brake = a_brake * (1 + e_brake);
    % return_alt = a_brake * (1 - e_brake) - r_N
    va_post_brake = h_brake/ra_brake;
    e_brake
end

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
    beta = 325e6;

    sidereal_period_N = 16.11 * 3600; % 16.11 hours
    omega_N = [0;0;2 * pi/sidereal_period_N]; % rad/s
    v_atm = cross(omega_N, r);
    v_rel = v - v_atm;

    alt = norm(r) - R;
    rho = atmos(alt);

    a_drag = -0.5 * rho * norm(v_rel) * v_rel / beta;
end


%% Pertubed propogation


function dY = eom_perturbed(~, Y, mu, R, R_equ, atmos, perturbations)
    r = Y(1:3);
    v = Y(4:6);
    r_norm = norm(r);

    a_total = -mu * r / r_norm^3;

    for k = 1:length(perturbations)
        switch perturbations{k}
            case "J2"
                a_total = a_total + J2_acc(r, mu, R_equ);
            case "drag"
                a_total = a_total + drag_acc(r, v, R, atmos);
            otherwise
                warning("Unknown perturbation: %s", perturbations{k});
        end
    end

    dY = [v; a_total];
end