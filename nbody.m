clear all
clc

% extract data
data = load("trajectory.mat");
r_sc = data.r_out;
r_planets = data.r;
n = size(r_sc, 1);
r_planets = r_planets(1:n, :);

% constants
mu = 1.32712e11; % for the sun
mu_planets = [3.24859e5, 3.98600e5, 4.28284e4, 1.26687e8, 3.7912e7, 5.79394e6, 6.83653e6];

% =========================================================================

% actual code
% Load trajectory data
data = load("trajectory.mat");  % assumes r_planets and r_sc are in here
r_planets = data.r;
r_sc = data.r_out;

% Constants
mu_sun = 1.32712e11; % km^3/s^2, gravitational parameter of the Sun
mu_planets = [3.24859e5, 3.98600e5, 4.28284e4, 1.26687e8, ...
              3.7912e7, 5.79394e6, 6.83653e6];  % km^3/s^2

% Initial spacecraft state from existing trajectory
r1 = r_sc(2, :)';          % initial position (km). avoid starting inside earth!
v1 = (r_sc(2, :) - r_sc(1, :))' / 86400;  % finite-diff estimate of velocity (km/s)

% Propagation duration in days
tf_days = 1000;

% Call n-body propagator
r_propagated = propagate_orbit_nbody(r1, v1, tf_days, mu_sun, mu_planets, r_planets);

% Plot comparison
figure;
plot3(r_sc(1:tf_days+1, 1), r_sc(1:tf_days+1, 2), r_sc(1:tf_days+1, 3), 'b', 'LineWidth', 1.5);
hold on;
plot3(r_propagated(:,1), r_propagated(:,2), r_propagated(:,3), 'r--', 'LineWidth', 1.5);
legend('Original r_{sc}', 'Propagated (n-body)');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
grid on;
axis equal;
title('Spacecraft Orbit Propagation (N-body vs Original)');

% =========================================================================

% propagation function
function rv_out = propagate_orbit_nbody(r1, v1, tf_days, mu_sun, mu_planets, r_planets)
    % tf is in days
    % Time vector in seconds
    tspan = (0:1:tf_days) * 86400;
    y0 = [r1(:); v1(:)];
    num_planets = length(mu_planets);
    
    % dynamics model
    function dydt = nbody_dynamics(t, y)
        r_curr = y(1:3);  % this is the current (propagated) position
        v_curr = y(4:6);
    
        epsilon = 1e-6;
    
        % Sun gravity
        r_norm = norm(r_curr);
        if r_norm < epsilon
            a_total = [0; 0; 0];
        else
            a_total = -mu_sun * r_curr / r_norm^3;
        end
    
        % Planet positions at this time step
        day_idx = floor(t / 86400) + 1;
        day_idx = min(day_idx, size(r_planets, 1));
        r_p = r_planets(day_idx, :);
    
        for i = 1:num_planets
            r_planet_i = r_p(1, 3*i-2:3*i);
            diff_r = r_planet_i' - r_curr;
            d_norm = norm(diff_r);
            if d_norm < epsilon
                warning('Spacecraft is nearly at planet %d position — skipping acceleration.', i);
                continue;
            end
            a_total = a_total + mu_planets(i) * diff_r / d_norm^3;
        end
    
        dydt = [v_curr; a_total];
    end

    options = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep', 86400);
    [~, Y] = ode113(@nbody_dynamics, tspan, y0, options);

    rv_out = Y;
end