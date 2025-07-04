clear all
clc
close all

% aiming radius
a = 638534.3331; % km
e = 2.0125;
aiming_radius = a * sqrt(e^2 - 1);
turn_angle = 2 * asin(1 / e);
r_p = a * (e - 1);
theta_inf = acos(-1 / e);
beta = pi - theta_inf;
x = aiming_radius /  sin(beta);

% extract data
data = load("trajectory.mat");
r_sc = data.r_out;
r_planets = data.r;
vector = r_sc(467, :) - r_sc(466, :);
r_sc(467:end, :) = r_sc(467:end, :) + repmat(vector,length(r_sc)-466,1);
r_sc(466, :) = [];
n = size(r_sc, 1);
r_planets = r_planets(1:n, :);
vec1 = r_sc(465, :) - r_sc(464, :);
vec2 = r_sc(467, :) - r_sc(466, :);
x_vector = (vec1 + vec2) / norm(vec1 + vec2) * x;
r_planets(:, 10:12) = r_planets(:, 10:12) - x_vector;

% constants
mu = 1.32712e11; % for the sun
mu_planets = [3.24859e5, 3.98600e5, 4.28284e4, 1.26687e8, 3.7912e7, 5.79394e6, 6.83653e6];
radius = [6051.8, 6378.0, 3389.5, 69911, 58232, 25362, 24622]; %in km
c = 2.998e8;                % speed of light m/s^2

% Cannonball model values
S0  = 63.15e6;           % surface radiated power intensity (W/m^2)
R_S  = 696340;           % radius of the Sun (km)
R    = 3;                % radius of spacecraft (m)
C_R  = 1.25;             % radiation pressure coefficient (1–2)
m = 1740 + 140 + 21115;           % mass of spacecraft (orbiterx2 + cruise stage)
m = 42500;
v    = zeros(n, 7);      % shadow function, size n×7
A_s  = pi * R^2;         % cross-sectional area of spacecraft (m^2)


% =========================================================================

% Constants
mu_sun = 1.32712e11; % km^3/s^2, gravitational parameter of the Sun
mu_planets = [3.24859e5, 3.98600e5, 4.28284e4, 1.26687e8, ...
              3.7912e7, 5.79394e6, 6.83653e6];  % km^3/s^2

% Initial spacecraft state from existing trajectory
r1 = r_sc(2, :)';          % initial position (km). avoid starting inside earth!
v1 = (r_sc(2, :) - r_sc(1, :))' / 86400;  % finite-diff estimate of velocity (km/s)

% Propagation duration in days
tcm_interval = 180; % TCM once every 30 days
time_elapsed = 0;
dv_total = 0;
big_rv_propagated = [r1' v1'];
dv_list = [];
tcm_list = [];
targets = [];
tcm_locs = [];

% Call n-body propagator
while time_elapsed < n
    tf_days = min(tcm_interval, n - time_elapsed - 1); % so the propagation length is not too long
    if tf_days == 0
        break
    end
    if 465 - time_elapsed < tcm_interval && time_elapsed < 465
        tf_days = 465 - time_elapsed;
    end
    rv_propagated = propagate_orbit_nbody(r1, v1, tf_days, mu_sun, mu_planets, r_planets,radius,S0,R_S,c,C_R,A_s,m);
    r1 = rv_propagated(end, 1:3);
    v1 = rv_propagated(end, 4:6);
    big_rv_propagated = [big_rv_propagated; rv_propagated];
    tcm_list = [tcm_list time_elapsed + 1];

    % test
    % correction_vector = r_sc(time_elapsed, :) - r1; % vector in correction direction
    % dv = correction_vector / sqrt(norm(correction_vector)) * 1e-6;
    % dv_total = dv_total + norm(dv);
    % v1 = v1 + dv;

    % lambert solver test
    time_elapsed = time_elapsed + tf_days;
    if 466 - time_elapsed < tcm_interval && time_elapsed < 465
        tf_days2 = 465 - time_elapsed;
    else
        tf_days2 = tcm_interval;
    end
    endpos = r_sc(min(n, time_elapsed + tf_days2), :);
    targets = [targets; endpos];
    tcm_locs = [tcm_locs; r1];
    [V1, ~] = lambert2(r1, endpos, tf_days2, 0, mu_sun);
    dv = (V1 - v1);
    if norm(dv) < 0.1
        dv_total = dv_total + norm(dv);
        dv_list = [dv_list; dv];
    end
    v1 = v1 + dv;
end

% Plot comparison
figure;
plot3(r_sc(:, 1), r_sc(:, 2), r_sc(:, 3), "magenta", 'LineWidth', 1.5);
hold on;
plot3(big_rv_propagated(:,1), big_rv_propagated(:,2), big_rv_propagated(:,3), "cyan", 'LineWidth', 1.5);
% scatter3(r_sc(tcm_list,1), ...
%     r_sc(tcm_list,2), ...
%     r_sc(tcm_list,3), 'kx', 'LineWidth', 1.5)
scatter3(tcm_locs(:, 1), tcm_locs(:, 2), tcm_locs(:, 3), 'ro', 'filled')
legend('Unperturbed trajectory', 'Perturbed trajectory', "TCMs");
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
grid on;
axis equal;
%title('Spacecraft Orbit Propagation (N-body vs Original)');
set(gca, 'Clipping', 'off');
% =========================================================================

% plot dv scatter plot
dv_list_mag = vecnorm(dv_list, 2, 2) * 1000; % in m/s
figure
semilogy(1:14, dv_list_mag, "kx", "LineWidth", 5)
grid on
xlabel("TCM number")
ylabel("TCM Δv (m/s)")
set(gca, "FontSize", 16)
xlim([0, 15])

% propagation function
function rv_out = propagate_orbit_nbody(r1, v1, tf_days, mu_sun, mu_planets, r_planets, radius, S0,R_S,c,C_R,A_s,m)
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
            a_nbody = [0; 0; 0];
        else
            a_nbody = -mu_sun * r_curr / r_norm^3;
        end
    
        % Planet positions at this time step
        day_idx = floor(t / 86400) + 1;
        day_idx = min(day_idx, size(r_planets, 1));
        r_p = r_planets(day_idx, :);
    
        for i = 1:num_planets
            r_planet_i = r_p(1, 3*i-2:3*i);
            diff_r = r_planet_i' - r_curr;
            d_norm = norm(diff_r);
            R_curr = norm(r_curr);
            R_planet_i = norm(r_planet_i);
            if d_norm < epsilon
                warning('Spacecraft is nearly at planet %d position — skipping acceleration.', i);
                continue;
            end
            a_nbody = a_nbody + mu_planets(i) * diff_r / d_norm^3;
            
            R_b = vecnorm(r_planet_i, 2 , 2 );       %in km

            theta = acosd((dot(r_curr',r_planet_i , 2 ))./(R_curr*R_planet_i));
            theta1 = acosd(radius(i)./R_curr);
            theta2 = acosd(radius(i)./R_planet_i);

            sight = (theta1 + theta2) >= theta;  % n×1 logical
            v(sight, i) = 1;

            S = S0*(R_S./R_curr).^2;      % solar constant in W/m^2
            P = S/c;                    % solar radiation pressure N/m^2

        end
        v = double(any(v, 2));

        F = -v.*P.*C_R.*A_s;                    %perturbing force in N
        a_srp = -(F/m).*((r_curr./R_curr)*0.001);   %perturbing acceleration vector in km/s^2

        a_total = a_nbody + 0*a_srp;
        dydt = [v_curr; a_total];
    end

    options = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep', 86400);
    [~, Y] = ode113(@nbody_dynamics, tspan, y0, options);

    rv_out = Y;
end