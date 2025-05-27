clear
clc
close all

%% Neptune Constants

mu_neptune = 6.8351e6; % km^3/s^2
mass_neptune = 102.409e24; % kg
radius_neptune = 24764; % km

%% Unperturbed Trajectory

if exist('vectors.mat', 'file') == 2
    load("vectors.mat")
    r0 = r/1000;
    v0 = v/1000;
else
    error('no intitial conditions')
end

% [r0, v0] = rv_from_oe(3e5, 0.9, 90, 90, 250, 0, mu_neptune);
[a, e, h, i, Omega, omega, theta] = find_OE(r0, v0, mu_neptune);


delta_theta = 1;
thetas = mod(theta:delta_theta:theta+360, 360);

r_unperturbed = zeros(3, length(thetas));
v_unperturbed = zeros(3, length(thetas));

r_unperturbed(:,1) = r0;
v_unperturbed(:,1) = v0;

for j = 2:length(thetas)
    [r_j, v_j] = rv_from_oe(a, e, i, Omega, omega, thetas(j), mu_neptune);
    r_unperturbed(:,j) = r_j;
    v_unperturbed(:,j) = v_j;
end

%% Perturbed trajectory

% Initial state vector (km and km/s)
Y0 = [r_unperturbed(:,1); v_unperturbed(:,1)];

% Orbital period (estimate from semi-major axis)
P = 2 * pi * sqrt(a^3 / mu_neptune);  % seconds
tspan = [0, P];         % time steps

% Integrate perturbed motion
opts = odeset('RelTol',1e-9, 'AbsTol',1e-9);
[t_out, Y_out] = ode45(@(t,Y) eom_perturbed(t, Y, mu_neptune, radius_neptune), tspan, Y0, opts);

r_perturbed = Y_out(:,1:3)';
v_perturbed = Y_out(:,4:6)';

%% Plot trajectories

figure;
plot3(r_unperturbed(1,:), r_unperturbed(2,:), r_unperturbed(3,:), 'b-');  % unperturbed
hold on;
plot3(r_perturbed(1,:), r_perturbed(2,:), r_perturbed(3,:), 'r--', 'LineWidth', 1.5);  % perturbed

% Mark start and end points
plot3(r_unperturbed(1,1), r_unperturbed(2,1), r_unperturbed(3,1), 'bo', 'MarkerFaceColor', 'b');
text(r_unperturbed(1,1), r_unperturbed(2,1), r_unperturbed(3,1), '  Start', 'Color', 'b');

plot3(r_unperturbed(1,end), r_unperturbed(2,end), r_unperturbed(3,end), 'bs');
text(r_unperturbed(1,end), r_unperturbed(2,end), r_unperturbed(3,end), '  End (Unperturbed)', 'Color', 'b');

plot3(r_perturbed(1,end), r_perturbed(2,end), r_perturbed(3,end), 'rs');
text(r_perturbed(1,end), r_perturbed(2,end), r_perturbed(3,end), '  End (Perturbed)', 'Color', 'r');

legend('Unperturbed', 'Perturbed (J2)', 'Location', 'best');
grid on;
axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('3D Trajectory with J2 Perturbation');

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
    [a_k, e_k, h_k, i_k, Omega_k, omega_k, theta_k] = find_OE(r_k, v_k, mu_neptune);
    a_vals(k) = a_k;
    e_vals(k) = e_k;
    i_vals(k) = i_k;
    Omega_vals(k) = Omega_k;
    omega_vals(k) = omega_k;
    theta_vals(k) = theta_k;
end

% Convert time to days
t_days = t_out / 86400;

% Plot
figure;
subplot(3,2,1); plot(t_days, a_vals);       ylabel('a [km]'); title('Semi-major Axis');
subplot(3,2,2); plot(t_days, e_vals);       ylabel('e'); title('Eccentricity');
subplot(3,2,3); plot(t_days, (i_vals));      ylabel('i [deg]');
subplot(3,2,4); plot(t_days, (Omega_vals));  ylabel('\Omega [deg]');
subplot(3,2,5); plot(t_days, (omega_vals));  ylabel('\omega [deg]');
subplot(3,2,6); plot(t_days, (theta_vals));  ylabel('\theta [deg]'); xlabel('Time [days]');
sgtitle('Orbital Elements Over Time (Perturbed)');

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

%% Perturbed equations of motion

function dY = eom_perturbed(~, Y, mu, R)
    r = Y(1:3);
    v = Y(4:6);
    r_norm = norm(r);

    a_gravity = -mu * r / r_norm^3;
    a_J2 = J2_acc(r, mu, R);

    a_total = a_gravity + a_J2;

    dY = [v; a_total];
end
