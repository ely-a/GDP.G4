% REQUIRED ADD-ONS
% type "aeroDataPackage" into the command window and get the planetary data

close all
clc

t_start = juliandate(2031,2,9);
t_end = t_start + 1e4;

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

planetList = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"];

xrange = [-5e9 5e9]; % this will be for neptune's orbit
yrange = [-5e9 5e9];
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



% iterate lambert's problem to get graphs
mu = 1.32712e11; % for the sun
mu_planets = [3.24859e5, 3.98600e5, 4.28284e4, 1.26687e8, 3.7912e7, 5.79394e6, 6.83653e6];
r_planets = [6051.8, 6378.1, 3396.2, 71492, 60268, 25559, 24764];
initial_v_inf = 8.59; % km/s
sequence = [2, 3, 4, 7]; % 1 = Venus --> 7 = Neptune
max_lengths = [1000, 2000, 5000]; % transfer times (days)
initial_LEO_alt = 300; %km
final_neptune_alt = 1000; %km
dv_launch = sqrt(initial_v_inf^2 + 2 * mu_planets(sequence(1)) /...
    (initial_LEO_alt + r_planets(sequence(1)))) - sqrt(mu_planets(sequence(1)) /...
    (initial_LEO_alt + r_planets(sequence(1))));
for time_elapsed = 0:1:0
    r_out = zeros(1, 3);
    disp(strcat("Time past initial launch date: ", num2str(time_elapsed), " days"))
    failed = false;
    initial_t = time_elapsed;
    flyby_vlist = zeros(1, 3);
    
    for num = 1:length(sequence)-1
        planet = sequence(num);
        tf_list = (10:max_lengths(num))'; % times of flight to check;
        visviva = getVisViva(r, v, sequence(num), sequence(num+1), tf_list, mu, time_elapsed);
        % vis-viva matching
        y1 = visviva(:, 1);
        y2 = visviva(:, 2);
        if num == 1
            y_target = initial_v_inf ^ 2;
        else
            y_target = y2_target; % use the old v_inf
        end
        % find the tof which gives the closest v_inf
        % Initialize result arrays
        x_targets = [];
        y2_targets = [];
        for i = 1:length(tf_list)-1
            if (y1(i) - y_target) * (y1(i+1) - y_target) < 0
                % Linear interpolation factor
                t = (y_target - y1(i)) / (y1(i+1) - y1(i));
                
                % Interpolate x and y2 at crossing point
                x_target = tf_list(i) + t * (tf_list(i+1) - tf_list(i));
                y2_target = y2(i) + t * (y2(i+1) - y2(i));
                
                % Store results
                x_targets(end+1) = x_target;
                y2_targets(end+1) = y2_target;
            end
            if length(x_targets) > 1
                x_target = x_targets(y2_targets == max(y2_targets));
                y2_target = max(y2_targets);
            end
            if i == length(tf_list)-1 && isempty(x_targets)
                disp(strcat("Spacecraft could not reach ", planetList(sequence(num+1))))
                failed = true;
                break
            end
        end

        % ========= comment out if iterating ==========
        % figure
        % semilogy(tf_list, visviva(:, 1))
        % hold on
        % semilogy(tf_list, visviva(:, 2))
        % grid on
        % xlabel("Time of Flight (days)")
        % ylabel("Vis-viva at intercept (km^2/s^2)")
        % legend("Departure", "Arrival")
        % xline(x_target)
        % 
        %data output
        disp("--------------------------------")
        disp(strcat(planetList(sequence(num)), " to ", planetList(sequence(num+1))))
        disp(strcat("Time of flight: ", num2str(x_target), " days"))
        disp(strcat("Arrival V_inf: ", num2str(y2_target ^ 0.5), " km/s"))
        % ===============================================
    
        % plotting
        p1 = sequence(num);
        p2 = sequence(num+1);
        idx1 = (3*p1-2):3*p1; % to find the right ephimeris, e.g. 4:6 for planet 2
        idx2 = (3*p2-2):3*p2;
        tf = floor(x_target);
        % position and velocity of the flyby planet
        t_floor = floor(time_elapsed);
        t_ceil = ceil(time_elapsed);
        alpha = time_elapsed - t_floor;  % Fractional part for interpolation
        % Get lower and upper vectors
        r1_lower = r(1 + t_floor, idx1);
        r1_upper = r(1 + t_ceil, idx1);
        r2_lower = r(1 + tf + t_floor, idx2);
        r2_upper = r(1 + tf + t_ceil, idx2);
        v1_lower = v(1 + t_floor, idx1);
        v1_upper = v(1 + t_ceil, idx1);
        v2_lower = v(1 + tf + t_floor, idx2);
        v2_upper = v(1 + tf + t_ceil, idx2);
        % Linearly interpolate each
        r1 = (1 - alpha) * r1_lower + alpha * r1_upper;
        r2 = (1 - alpha) * r2_lower + alpha * r2_upper;
        v1 = (1 - alpha) * v1_lower + alpha * v1_upper;
        v2 = (1 - alpha) * v2_lower + alpha * v2_upper;
        % propagation
        [V1, V2] = lambert2(r1, r2, tf, 0, mu);
        r_leg = propagate_orbit(r1, V1, tf, mu);
        r_out = [r_out; r_leg];
        % relative velocity to calculate flyby parameters
        flyby_vlist = [flyby_vlist; V1-v1; V2-v2];

        time_elapsed = time_elapsed + x_target;

        % if failed, move on
        if failed
            break
        end
    end
    if ~failed
        disp("====================== Successful trajectory! ======================")
        total_time = time_elapsed - initial_t;
        years = floor(total_time / 365);
        days = round(total_time - years*365);
        disp(strcat("Total time: ", num2str(years), "y ", num2str(days), "d"))
        launch_JD = t_start + initial_t;
        launch_date = datetime(launch_JD, 'ConvertFrom', 'juliandate');
        disp(strcat("Launch Date: ", string(launch_date)));
        dv_capture = sqrt(y2_target + 2 * mu_planets(sequence(end)) /...
        (final_neptune_alt + r_planets(sequence(end)))) - sqrt(mu_planets(sequence(end)) /...
        (final_neptune_alt + r_planets(sequence(end))));
        total_dv = dv_launch + dv_capture;
        disp(strcat("Total ΔV: ", num2str(total_dv), " km/s"))
        disp("====================== Flyby information ======================")

        % flyby calculations
        flyby_vlist(1, :) = []; % remove initial zeroes
        for flyby = 1:length(sequence)-2 % as no flyby occurs at Earth or Neptune
            initialV = flyby_vlist(2*flyby, :);
            finalV = flyby_vlist(2*flyby+1, :);
            turn_angle = acos(dot(initialV, finalV) / norm(initialV) / norm(finalV));
            ecc = 1 / sin(turn_angle / 2);
            sma = mu_planets(sequence(flyby+1)) / norm(finalV)^2;
            r_p = sma * ecc-1;
            radius = r_planets(sequence(flyby+1));
            alt = r_p - radius;
            radius_number = r_p / radius;
            planet_name = planetList(sequence(flyby+1));
            disp(strcat(planet_name, " flyby: Altitude = ", num2str(alt), " km (", num2str(radius_number), " radii)"))
        end

    end
end

figure
r_out(1, :) = [];
plot3(r_out(:, 1), r_out(:, 2), r_out(:, 3), color="black", LineWidth=1, DisplayName="Spacecraft")
axis equal
grid on

% Planet names for labeling

hold on
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');

% Plot each planet's trajectory
for i = 1:length(planetList)
    idx = (i - 1) * 3 + 1; % Column index for X
    plot3(r(:, idx), r(:, idx+1), r(:, idx+2), 'DisplayName', planetList(i), ...
          'Color', colors(i,:), 'LineWidth', 1.5);
end
legend show;

% Initialize the figure
figure;
axis equal;

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
view(3);
hold on;

% add s/c to plotting array
r_truncated = r(1:size(r_out,1), :);
r_plot = [r_truncated r_out];

nSteps = size(r_plot, 1);
nPlanets = 8;
trail = gobjects(nPlanets,1);
dot = gobjects(nPlanets,1);

for i = 1:8
    trail(i) = plot3(NaN, NaN, NaN, 'Color', colors(i,:), 'LineWidth', 1.5);
    dot(i) = plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 6, ...
                   'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
end

% Animation loop
% % Create a VideoWriter object
% v = VideoWriter('trajectory.mp4', 'MPEG-4'); % You can also use 'Motion JPEG AVI'
% v.FrameRate = 60; % Set frame rate (optional)
% open(v); % Open the video file for writing

for k = 1:20:nSteps  % Skip some steps to speed up animation
    for i = 1:8
        idx = (i-1)*3 + 1;
        % Update trails
        set(trail(i), 'XData', r_plot(1:k, idx), ...
                      'YData', r_plot(1:k, idx+1), ...
                      'ZData', r_plot(1:k, idx+2));
        % Update dots
        set(dot(i), 'XData', r_plot(k, idx), ...
                    'YData', r_plot(k, idx+1), ...
                    'ZData', r_plot(k, idx+2));
    end
    drawnow;
    axis equal
    grid on
    view(2)
    xlim(xrange);
    ylim(yrange);
    zlim(zrange);
    pause(1e-6)

    % Capture the plot as a frame
    % frame = getframe(gcf);
    % writeVideo(v, frame);

end
% close(v);


% for plotting only
function r_out = propagate_orbit(r1, v1, tf, mu)
    % Inputs:
    % r1 - 3x1 position vector (km)
    % v1 - 3x1 velocity vector (km/s)
    % tf - total propagation time (s)
    % m  - gravitational parameter μ (km^3/s^2)

    % Time vector from 0 to tf (1-JD intervals)
    tspan = (0:1:tf) * 86400;

    % Initial state vector [r; v]
    y0 = [r1(:); v1(:)];

    % Define two-body dynamics
    function dydt = two_body(~, y)
        r = y(1:3);
        v = y(4:6);
        r_norm = norm(r);
        a = -mu * r / r_norm^3;
        dydt = [v; a];
    end

    % Integrate using ODE113 (accurate for orbital dynamics)
    options = odeset('RelTol',1e-10,'AbsTol',1e-12);
    [~, Y] = ode113(@two_body, tspan, y0, options);

    % Extract position vectors
    r_out = Y(:, 1:3);
end



% vis viva function
function visviva = getVisViva(r, v, p1, p2, tf_list, mu, t_elapsed)
% r, v are the planet ephimeris
% p1 is the planet number for the first planet (Venus - Neptune)
% p2 is the planet number for the second planet (Venus - Neptune)
% tf_list is an array of the flight times to check
    % iterate lambert's problem to get graphs
    visviva = zeros(length(tf_list), 2);
    i = 1;
    for idx = 1:length(tf_list)
        tf = tf_list(idx);
        idx1 = (3*p1-2):3*p1; % to find the right ephimeris, e.g. 4:6 for planet 2
        idx2 = (3*p2-2):3*p2;
        r1 = r(1+floor(t_elapsed), idx1);
        r2 = r(1+tf+floor(t_elapsed), idx2);
        [V1, V2] = lambert2(r1, r2, tf, 0, mu);
        % change to relative velocities
        VV1 = norm(V1 - v(1+floor(t_elapsed), idx1))^2;
        VV2 = norm(V2 - v(1+tf+floor(t_elapsed), idx2))^2;
        visviva(i, :) = [VV1 VV2];
        i = i + 1;
    end
end