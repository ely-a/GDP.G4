clear
clc
close all

%% Neptune Constants

mu_neptune = 6.8351e6; % km^3/s^2
mass_neptune = 102.409e24; % kg
radius_neptune = 24764; % km

%% Unperturbed Trajectory

r0 = [-6.0955;-11.8452;11.0082];   % import from file later                        
v0 = [-0.1653;-0.0018;-0.1474]; 

mu_neptune = 1;

[a, e, h, i, Omega, omega, theta] = find_OE(r0, v0, mu_neptune);

P = 2*pi*sqrt(a^3/mu_neptune);
timesteps = linspace(0, P, 1000);
delta_t = timesteps(2) - timesteps(1);

r_unperturbed = zeros(3, length(timesteps));
v_unperturbed = zeros(3, length(timesteps));

r_unperturbed(:,1) = r0;
v_unperturbed(:,1) = v0;

for j = 2:length(timesteps)
    [r0, v0] = propogate_lagrangian(r0, v0, mu_neptune, delta_t);
    r_unperturbed(:,j) = r0;
    v_unperturbed(:,j) = v0;
end

%% Plot trajectories

plot3(r_unperturbed(:,1), r_unperturbed(:,2), r_unperturbed(:,3), 'b-');  % blue line
grid on;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Trajectory');

%% Enckes Method
% r0 and v0 get from elya import % km

% Propogte Lagrangian
% prop_lagrangian(t)

% Calculate perturbing acceleration
% J2_acc(r_0)

% Integrate delta
% ode45

% calcaulte perturbed trajectory
% r = r0 + delta
% r0 = r and repeat

%% J2 Perturbations

function a_J2 = J2_acc(r_vec, mu, R)
    a_J2 = zeros(3,1);
    x = r(1);
    y = r(2);
    z = r(3);
    r = vecnorm(r_vec);

    J2 = 3.411e-3;
    J2_const = (1.5 * J2 * mu * R^2)/(r^5);

    a_J2(1,1) = J2_const * x * (5 * (z/r)^2 - 1);
    a_J2(2,1) = J2_const * y * (5 * (z/r)^2 - 1);
    a_J2(3,1) = J2_const * z * (5 * (z/r)^2 - 3);
end