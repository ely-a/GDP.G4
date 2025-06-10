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

%% Earth Constants

mu_Earth = 398600; % km^3/s^2
mass_Earth = 5.972e24; % kg
radius_Earth = 6378; % km

%% Unperturbed Trajectory

if exist('vectors_Earth.mat', 'file') == 2
    load("vectors_Earth.mat")
    r0 = r0_sc;
    v0 = v0_sc;
else
    error('no intitial conditions')
end


rp = 200 + radius_Earth;
e = 0.96;
a = rp / (1 - e);

N_AR = 20;
e_AR = linspace(0,e,N_AR);
a_AR = rp ./ (1 - e_AR);


h = sqrt(mu_Earth * a * (1 - e^2));
i = 28.5;
Omega = 18.3;
omega = 359.2;
theta = 0;

r_perturbed_big = [];
r_unperturbed_big = [];
i_big = [];
Omega_big = [];
omega_big = [];
t_days_big = [0];

for idx = 1:N_AR
    [r0, v0] = rv_from_oe(a_AR(idx), e_AR(idx), i, Omega, omega, theta, mu_Earth);
    %[a, e, h, i, Omega, omega, theta] = find_OE(r0, v0, mu_Earth);
    P = 2 * pi * sqrt(a_AR(idx)^3 / mu_Earth);  % seconds
    N_orbits = 1;
    P_AR = 2 * pi * sqrt(a_AR(idx)^3 / mu_Earth);  % seconds

    delta_theta = 1;
    thetas = mod(theta:delta_theta:theta+360, 360);
    
    r_unperturbed = zeros(3, length(thetas));
    v_unperturbed = zeros(3, length(thetas));
    
    r_unperturbed(:,1) = r0;
    v_unperturbed(:,1) = v0;
    
    for j = 2:length(thetas)
        [r_j, v_j] = rv_from_oe(a_AR(idx),  e_AR(idx), i, Omega, omega, thetas(j), mu_Earth);
        r_unperturbed(:,j) = r_j;
        v_unperturbed(:,j) = v_j;
    end
    
    %% Perturbed trajectory
    
    % Initial state vector (km and km/s)
    Y0 = [r_unperturbed(:,1); v_unperturbed(:,1)];
    
    % Orbital period (estimate from semi-major axis)
    tspan = [0, P];         % time steps
    
    % Integrate perturbed motion
    opts = odeset('RelTol',1e-9, 'AbsTol',1e-9);
    perturbs = ["J2"];
    [t_out, Y_out] = ode113(@(t,Y) eom_perturbed(t, Y, mu_Earth, radius_Earth, perturbs) ...
                            , tspan, Y0, opts);
    
    r_perturbed = Y_out(:,1:3)';
    v_perturbed = Y_out(:,4:6)';
    
    v_perturbed_peri = max(vecnorm(v_perturbed(:,30:end)));
    dV(idx) = abs(v_perturbed_peri - norm(v_perturbed(:,1)))*1e3;

    r_perturbed_big = [r_perturbed_big r_perturbed];
    r_unperturbed_big = [r_unperturbed_big r_unperturbed];

    
    % %% Plot trajectories
    % 
    % figure;
    % plot3(r_unperturbed(1,:), r_unperturbed(2,:), r_unperturbed(3,:), 'b-');  % unperturbed
    % hold on;
    % plot3(r_perturbed(1,:), r_perturbed(2,:), r_perturbed(3,:), 'r--', 'LineWidth', 1.5);  % perturbed
    % 
    % % Mark start and end points
    % plot3(r_unperturbed(1,1), r_unperturbed(2,1), r_unperturbed(3,1), 'bo', 'MarkerFaceColor', 'b');
    % text(r_unperturbed(1,1), r_unperturbed(2,1), r_unperturbed(3,1), '  Start', 'Color', 'b');
    % 
    % plot3(r_unperturbed(1,end), r_unperturbed(2,end), r_unperturbed(3,end), 'bs');
    % text(r_unperturbed(1,end), r_unperturbed(2,end), r_unperturbed(3,end), '  End (Unperturbed)', 'Color', 'b');
    % 
    % plot3(r_perturbed(1,end), r_perturbed(2,end), r_perturbed(3,end), 'rs');
    % text(r_perturbed(1,end), r_perturbed(2,end), r_perturbed(3,end), '  End (Perturbed)', 'Color', 'r');
    % 
    % legend('Unperturbed', 'Perturbed (J2)', 'Location', 'best');
    % grid on;
    % axis equal;
    % xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    % title('3D Trajectory with J2 Perturbation');
    
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
        %[a_k, e_k, h_k, i_k, Omega_k, omega_k, theta_k] = find_OE(r_k, v_k, mu_Earth);
        [coe] = sv_2_oe(r_k, v_k, mu_Earth);
        a_vals(k) = coe(7);
        e_vals(k) = coe(2);
        i_vals(k) = coe(4);
        Omega_vals(k) = coe(3);
        omega_vals(k) = coe(5);
        theta_vals(k) = coe(6);
        % a_vals(k) = a_k;
        % e_vals(k) = e_k;
        % i_vals(k) = i_k;
        % Omega_vals(k) = Omega_k;
        % omega_vals(k) = omega_k;
        % theta_vals(k) = theta_k;
    end
    
    % Convert time to days
    t_days = t_out / 86400;
    t_days_big = [t_days_big; t_days+t_days_big(end)];
    omega_vals(omega_vals > 180) = omega_vals(omega_vals > 180) - 360;
    
    i_big = [i_big i_vals];
    Omega_big = [Omega_big Omega_vals];
    omega_big = [omega_big omega_vals];
    
    
    % Plot
    figure;
    subplot(3,2,1); plot(t_days, a_vals); hold on;       ylabel('a [km]'); title('Semi-major Axis');xlabel('Time [days]');
    subplot(3,2,2); plot(t_days, e_vals);       ylabel('e'); title('Eccentricity');xlabel('Time [days]');
    subplot(3,2,3); plot(t_days, (i_vals));      ylabel('i [deg]');xlabel('Time [days]');
    subplot(3,2,4); plot(t_days, (Omega_vals));  ylabel('\Omega [deg]');xlabel('Time [days]');
    subplot(3,2,5); plot(t_days, (omega_vals));  ylabel('\omega [deg]');xlabel('Time [days]');
    subplot(3,2,6); plot(t_days, (theta_vals));  ylabel('\theta [deg]'); xlabel('Time [days]');
    sgtitle('Orbital Elements Over Time (Perturbed)');
end

dV_mag = sum(dV)

% Plot trajectories

figure;
plot3(r_unperturbed_big(1,:), r_unperturbed_big(2,:), r_unperturbed_big(3,:), 'b-','LineWidth', 0.75);  % unperturbed
hold on;
plot3(r_perturbed_big(1,:), r_perturbed_big(2,:), r_perturbed_big(3,:), 'r--', 'LineWidth', 1.5);  % perturbed

%Mark start and end points
% plot3(r_unperturbed_big(1,1), r_unperturbed_big(2,1), r_unperturbed_big(3,1), 'bo', 'MarkerFaceColor', 'b');
% text(r_unperturbed_big(1,1), r_unperturbed_big(2,1), r_unperturbed_big(3,1), '  Start', 'Color', 'b');

% plot3(r_unperturbed_big(1,end), r_unperturbed_big(2,end), r_unperturbed_big(3,end), 'bs');
% text(r_unperturbed_big(1,end), r_unperturbed_big(2,end), r_unperturbed_big(3,end), '  End (Unperturbed)', 'Color', 'b');
% 
% plot3(r_perturbed_big(1,end), r_perturbed_big(2,end), r_perturbed_big(3,end), 'rs');
% text(r_perturbed_big(1,end), r_perturbed_big(2,end), r_perturbed_big(3,end), '  End (Perturbed)', 'Color', 'r');

legend('Unperturbed', 'Perturbed (J2)', 'Location', 'best');
grid on;
axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
%title('3D Trajectory with J2 Perturbation');


% Plot
figure;
subplot(1,3,1); plot(t_days_big(2:end), (i_big));      ylabel('i [deg]');xlabel("Days since launch");
subplot(1,3,2); plot(t_days_big(2:end), (Omega_big));  ylabel('\Omega [deg]');xlabel("Days since launch")
subplot(1,3,3); plot(t_days_big(2:end), (omega_big));  ylabel('\omega [deg]');xlabel("Days since launch")
sgtitle('Orbital Elements Over Time (Perturbed)');

%% J2 Omega and omega rate of change comparison
% === Analytical J2 Drift Rates ===

J2 = 1.082636e-3;           % Earth's J2 coefficient
mu = mu_Earth;           % km^3/s^2
R = radius_Earth;        % km

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

%% J2 Perturbations

function a_J2 = J2_acc(r_vec, mu, R)
    x = r_vec(1);
    y = r_vec(2);
    z = r_vec(3);
    r = norm(r_vec);

    J2 = 1.082636e-3; 
    J2_const = (1.5 * J2 * mu * R^2)/(r^5);

    a_J2 = J2_const * [
        x * (5 * (z/r)^2 - 1);
        y * (5 * (z/r)^2 - 1);
        z * (5 * (z/r)^2 - 3)
    ];
end

%% Perturbed equations of motion

function dY = eom_perturbed(~, Y, mu, R, perturbations)
    r = Y(1:3);
    v = Y(4:6);
    r_norm = norm(r);

    a_total = -mu * r / r_norm^3;

    for pertubation = perturbations
        if pertubation == "J2"
            a_total = a_total + J2_acc(r, mu, R);
        end
    end

    dY = [v; a_total];
end
