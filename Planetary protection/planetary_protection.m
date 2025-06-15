clc

t_start = juliandate(2032,3,23);
t_end = t_start + 100 * 365;

% Define the list of variables to keep
dataVars = {'r', 'v'};
missingData = ~all(cellfun(@(v) evalin('base', sprintf('exist(''%s'', ''var'')', v)), dataVars));
if missingData
    clearvars -except t_start t_end
    t_query = (t_start:t_end)';
    r = zeros(length(t_query), 21);
    v = zeros(length(t_query), 21);
    disp("Getting Venus ephimeris")
    [r(:, 1:3), v(:, 1:3)] = planetEphemeris(t_query, "SolarSystem", "Venus", "430");
    disp("Getting Earth ephimeris")
    [r(:, 4:6), v(:, 4:6)] = planetEphemeris(t_query, "SolarSystem", "Earth", "430");
    disp("Getting Mars ephimeris")
    [r(:, 7:9), v(:, 7:9)] = planetEphemeris(t_query, "SolarSystem", "Mars", "430");
    disp("Getting Jupiter ephimeris")
    [r(:, 10:12), v(:, 10:12)] = planetEphemeris(t_query, "SolarSystem", "Jupiter", "430");
    disp("Getting Saturn ephimeris")
    [r(:, 13:15), v(:, 13:15)] = planetEphemeris(t_query, "SolarSystem", "Saturn", "430");
    disp("Getting Uranus ephimeris")
    [r(:, 16:18), v(:, 16:18)] = planetEphemeris(t_query, "SolarSystem", "Uranus", "430");
    disp("Getting Neptune ephimeris")
    [r(:, 19:21), v(:, 19:21)] = planetEphemeris(t_query, "SolarSystem", "Neptune", "430");
else
    % Keep data variables, clear everything else
    clearvars -except r v t_start t_end
end

r_planets = r;
trajectory = load("trajectory.mat");
r_sc = trajectory.r_out;

% constants
mu = 1.32712e11; % for the sun
mu_planets = [3.24859e5, 3.98600e5, 4.28284e4, 1.26687e8, 3.7912e7, 5.79394e6, 6.83653e6];
radius = [6051.8, 6378.0, 3389.5, 69911, 58232, 25362, 24622]; %in km
c = 2.998e8;                % speed of light m/s^2
n = 100 * 365;

% Cannonball model values
S0  = 63.15e6;           % surface radiated power intensity (W/m^2)
R_S  = 696340;           % radius of the Sun (km)
R    = 3;                % radius of spacecraft (m)
C_R  = 1.25;             % radiation pressure coefficient (1–2)
m    = 18000;            % mass of spacecraft (kg)
v    = zeros(n, 7);      % shadow function, size n×7
A_s  = pi * R^2;         % cross-sectional area of spacecraft (m^2)


% =========================================================================

% Constants
mu_sun = 1.32712e11; % km^3/s^2, gravitational parameter of the Sun
mu_planets = [3.24859e5, 3.98600e5, 4.28284e4, 1.26687e8, ...
              3.7912e7, 5.79394e6, 6.83653e6];  % km^3/s^2

% Initial spacecraft state from existing trajectory
v1 = (r_sc(2, :) - r_sc(1, :))' / 86400;  % finite-diff estimate of velocity (km/s)
nCases = 27;
r_plot = cell(nCases, 1);  % Column cell array

% Range of offsets for each component
delta = [-1, 0, 1] * 5e-4;
[dx, dy, dz] = ndgrid(delta, delta, delta);
offsets = [dx(:), dy(:), dz(:)];  % 27×3

v1 = v1(:)';  % Convert 3×1 to 1×3 row vector

% Add central velocity to each offset
velocities = bsxfun(@plus, offsets, v1);  % Result is 27×3

for i = 1:size(velocities, 1)
    disp(i)
    v1 = velocities(i, :);
    r1 = r_sc(2, :);          % initial position (km). avoid starting inside earth!

    % Propagation duration in days
    time_elapsed = 0;
    dv_total = 0;
    big_rv_propagated = [r1 v1];
    dv_list = [];
    tcm_list = [];
    targets = [];
    
    planetList = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"];
    
    xrange = [-5e8 5e8]; % this will be for neptune's orbit
    yrange = [-5e8 5e8];
    zrange = [-2e9 2e9];
    
    % Define planet colors
    colors = [
        0.95, 0.75, 0.10;  % Venus
        0.00, 0.45, 0.74;  % Earth
        0.85, 0.33, 0.10;  % Mars
        0.80, 0.60, 0.40;  % Jupiter
        0.93, 0.75, 0.50;  % Saturn
        0.40, 0.80, 0.90;  % Uranus
        0.30, 0.40, 0.85;  % Neptune
        0, 0, 0 % s/c
    ];
    
    rv_propagated = propagate_orbit_nbody(r1, v1, 100 * 365, mu_sun, mu_planets, r_planets,radius,S0,R_S,c,C_R,A_s,m);
    r1 = rv_propagated(end, 1:3);
    v1 = rv_propagated(end, 4:6);
    big_rv_propagated = [big_rv_propagated; rv_propagated];
    
    % add s/c to plotting array
    r_plot_array = rv_propagated(:, 1:3);
    r_plot{i} = r_plot_array;

end

figure;
hold on;
grid on;
axis equal;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Planetary Orbits and Multiple Spacecraft Trajectories');
view(3);

planetNames = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"];
colors = [
    0.95, 0.75, 0.10;  % Venus
    0.00, 0.45, 0.74;  % Earth
    0.85, 0.33, 0.10;  % Mars
    0.80, 0.60, 0.40;  % Jupiter
    0.93, 0.75, 0.50;  % Saturn
    0.40, 0.80, 0.90;  % Uranus
    0.30, 0.40, 0.85;  % Neptune
];

nPlanets = length(planetNames);
planet_handles = gobjects(nPlanets,1);

% Plot planet orbits and keep handles for legend
for i = 1:nPlanets
    idx = (i-1)*3 + 1;
    x = r_planets(:, idx);
    y = r_planets(:, idx+1);
    z = r_planets(:, idx+2);
    
    % Plot orbit as line
    planet_handles(i) = plot3(x, y, z, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Plot final position as scatter (no legend entry)
    scatter3(x(end), y(end), z(end), 50, colors(i,:), 'filled', 'HandleVisibility','off');
    
    text(x(end), y(end), z(end), ['      ' planetNames(i)], 'FontSize', 10, 'Color', colors(i,:));
end

% Plot spacecraft trajectories
nTraj = length(r_plot);
sc_handles = gobjects(1);

for j = 1:nTraj
    traj = r_plot{j};
    nCols = size(traj, 2);
    sc_x = traj(:, nCols-2);
    sc_y = traj(:, nCols-1);
    sc_z = traj(:, nCols);
    
    % Plot spacecraft trajectory (only create handle once for legend)
    if j == 1
        sc_handles = plot3(sc_x, sc_y, sc_z, 'Color', [0 0 0 0.2], 'LineWidth', 1);
    else
        plot3(sc_x, sc_y, sc_z, 'Color', [0 0 0 0.2], 'LineWidth', 1, 'HandleVisibility','off');
    end
end

% Create a clean legend: planets + one entry for all spacecraft trajectories
legend([planet_handles; sc_handles], [planetNames'; "Spacecraft Trajectories"], 'Location', 'bestoutside');

hold off;
% nSteps = size(r_plot, 1);
% nPlanets = 8;
% trail = gobjects(nPlanets,1);
% dot = gobjects(nPlanets,1);
% 
% figure;
% hold on;
% axis equal;
% 
% for i = 1:8
%     trail(i) = plot3(NaN, NaN, NaN, 'Color', colors(i,:), 'LineWidth', 1.5);
%     dot(i) = plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 6, ...
%                    'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
% end
% 
% % Animation loop
% % % Create a VideoWriter object
% % v = VideoWriter('trajectory.mp4', 'MPEG-4'); % You can also use 'Motion JPEG AVI'
% % v.FrameRate = 60; % Set frame rate (optional)
% % open(v); % Open the video file for writing
% 
% for k = 1:100:nSteps  % Skip some steps to speed up animation
%     for i = 1:8
%         idx = (i-1)*3 + 1;
%         % Update trails
%         set(trail(i), 'XData', r_plot(1:k, idx), ...
%                       'YData', r_plot(1:k, idx+1), ...
%                       'ZData', r_plot(1:k, idx+2));
%         % Update dots
%         set(dot(i), 'XData', r_plot(k, idx), ...
%                     'YData', r_plot(k, idx+1), ...
%                     'ZData', r_plot(k, idx+2));
%     end
%     drawnow;
%     axis equal
%     grid on
%     view(2)
%     xlim(xrange);
%     ylim(yrange);
%     zlim(zrange);
%     pause(1e-6)
% 
%     % Capture the plot as a frame
%     % frame = getframe(gcf);
%     % writeVideo(v, frame);
% 
% end
% % close(v);

% =========================================================================

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

    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [~, Y] = ode45(@nbody_dynamics, tspan, y0, options);
    rv_out = Y;
end