% REQUIRED ADD-ONS
% type "aeroDataPackage" into the command window and get the planetary data

close all
clc

t_start = juliandate(2033,6,1);
t_end = t_start + 2e4;

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
mu_planets = [3.24859e5, 3.98600e5, 4.28284e4, 1.26687e8, 3.7912e7, 5.79394e6, 6.83653e6]; % km^3/s^2
r_planets = [6051.8, 6378.1, 3396.2, 71492, 60268, 25559, 24764]; % km
sequence = [2, 7]; % 1 = Venus --> 7 = Neptune
max_lengths = 5000; % transfer times (days)
initial_LEO_alt = 300; %km
final_neptune_alt = 10000; %km
max_delay = 1; % days past initial launch window
big_list = zeros(max_lengths - 9, max_delay);
big_list_launch = zeros(max_lengths - 9, max_delay);
big_list_capture = zeros(max_lengths - 9, max_delay);
for time_elapsed = 0:max_delay-1 % adjust to test multiple launch days
    r_out = zeros(1, 3);
    disp(strcat("Time past initial launch date: ", num2str(time_elapsed), " days"))
    failed = false;
    initial_t = time_elapsed;
    flyby_vlist = zeros(1, 3);
    
    for num = 1:length(sequence)-1
        planet = sequence(num);
        tf_list = (10:max_lengths(num))'; % times of flight to check;
        [v1_list, visviva] = getVisViva(r, v, sequence(num), sequence(num+1), tf_list, mu, time_elapsed);
        % vis-viva matching
        y1 = visviva(:, 1);
        y2 = visviva(:, 2);
        dv_launch = sqrt(y1 + 2 * mu_planets(sequence(1)) /...
        (initial_LEO_alt + r_planets(sequence(1)))) - sqrt(mu_planets(sequence(1)) /...
        (initial_LEO_alt + r_planets(sequence(1))));
        r_p = final_neptune_alt + r_planets(sequence(end));     % Periapsis radius
        v_arrival = sqrt(y2 + 2 * mu_planets(sequence(end)) / r_p);  % Hyperbolic arrival speed at periapsis
        v_parabolic = sqrt(2 * mu_planets(sequence(end)) / r_p);     % Parabolic speed at same radius
        dv_capture = v_arrival - v_parabolic;   % Δv to just capture (into parabolic orbit)
        dv_total = dv_launch + dv_capture;

        %vis-viva matching graphs
        % figure
        % semilogy(tf_list, visviva(:, 1))
        % hold on
        % semilogy(tf_list, visviva(:, 2))
        % grid on
        % xlabel("Time of Flight (days)")
        % ylabel("Vis-viva at intercept (km^2/s^2)")
        % legend("Departure", "Arrival")
    end
    big_list(:, time_elapsed+1) = dv_total;
    big_list_launch(:, time_elapsed+1) = dv_launch;
    big_list_capture(:, time_elapsed+1) = dv_capture;
end
% figure
% semilogy(tf_list, dv_launch)
% hold on
% semilogy(tf_list, dv_capture)
% semilogy(tf_list, dv_capture+dv_launch)
% grid on
% legend("Launch dV", "Capture dV", "Total dV")
% xlabel("Days to reach Neptune")
% ylabel("dV (km/s)")

%[~, best_idx] = min(sum(visviva, 2));
best_idx = 365*11+2;
v1 = v1_list(best_idx, :);
tf_best = tf_list(best_idx);

% Get planet departure position
idx1 = (3*sequence(num)-2):3*sequence(num);
r1 = r(1+floor(time_elapsed), idx1);

% Propagate Lambert arc using two-body motion
r_sc_traj = propagate_orbit(r1, v1, tf_best, mu);  % Output: Nx3

% Grab planetary positions during transfer
r_planet_interp = r(1+floor(time_elapsed):1+floor(time_elapsed)+tf_best, :);  % Nx21

% Combine into r_plot
r_plot = [r_planet_interp, r_sc_traj];  % Nx(21+3)

% Plot
[nRows, nCols] = size(big_list);
x = 1:nCols;   % Time past initial launch date (days)
y = 1:nRows;   % Transfer time (days)

function coe = orbital_elements(r1, v1, mu)
% INPUTS:
% r1 - position vector [km]
% v1 - velocity vector [km/s]
% mu - gravitational parameter [km^3/s^2]

% OUTPUT:
% coe - structure containing orbital elements:
%       a     - semi-major axis [km]
%       e     - eccentricity
%       i     - inclination [deg]
%       RAAN  - right ascension of ascending node [deg]
%       omega - argument of periapsis [deg]
%       theta - true anomaly [deg]

% Magnitudes
r = norm(r1);
v = norm(v1);

% Specific angular momentum
h = cross(r1, v1);
h_mag = norm(h);

% Eccentricity vector
e_vec = (1/mu) * ((v^2 - mu/r)*r1 - dot(r1,v1)*v1);
e = norm(e_vec);

% Semi-major axis
energy = v^2/2 - mu/r;
if abs(e - 1) < 1e-6
    a = Inf;  % Parabolic (technically semi-latus rectum should be used)
else
    a = -mu / (2*energy);
end

% Inclination
i = acosd(h(3)/h_mag);

% Node vector
K = [0 0 1];
n = cross(K, h);
n_mag = norm(n);

% Right ascension of ascending node (RAAN)
if n_mag ~= 0
    RAAN = acosd(n(1)/n_mag);
    if n(2) < 0
        RAAN = 360 - RAAN;
    end
else
    RAAN = 0;
end

% Argument of periapsis (omega)
if n_mag ~= 0 && e > 1e-8
    omega = acosd(dot(n, e_vec)/(n_mag*e));
    if e_vec(3) < 0
        omega = 360 - omega;
    end
else
    omega = 0;
end

% True anomaly (theta)
if e > 1e-8
    theta = acosd(dot(e_vec, r1)/(e*r));
    if dot(r1, v1) < 0
        theta = 360 - theta;
    end
else
    % Circular orbit
    cp = cross(n, r1);
    if cp(3) >= 0
        theta = acosd(dot(n, r1)/(n_mag*r));
    else
        theta = 360 - acosd(dot(n, r1)/(n_mag*r));
    end
end

% Return as structure
coe = struct('a', a, 'e', e, 'i', i, ...
             'RAAN', RAAN, 'omega', omega, 'theta', theta);
end

figure;
imagesc(x, y, big_list);   % Specify x and y axes
colormap hot;
colorbar;
clim([0 50]);              % Focus on low dV
set(gca, 'YDir', 'normal');
xlabel("Time past initial launch date (days)");
ylabel("Transfer time (days)");
title("Total dV");

% Add contour lines
hold on;
levels = 0:2:50;  % Contour levels every 5 units between 0 and 50
[C, h] = contour(x, y, big_list, levels, 'LineColor', 'k');
clabel(C, h);  % Add labels to contours

% Plot
figure;
imagesc(big_list_launch);        % Display the array as an image
colorbar;          % Add colorbar to show value scale
colormap hot;
xlabel("Time past initial launch date (days)")
ylabel("Transfer time (days)")
clim([0 30]);  % Set your desired color range
set(gca, 'YDir', 'normal');
title("Earth departure dV")

% Add contour lines
hold on;
levels = 0:2:30;  % Contour levels every 5 units between 0 and 50
[C, h] = contour(x, y, big_list_launch, levels, 'LineColor', 'k');
clabel(C, h);  % Add labels to contours

% Plot
figure;
imagesc(big_list_capture);        % Display the array as an image
colorbar;          % Add colorbar to show value scale
colormap hot;
xlabel("Time past initial launch date (days)")
ylabel("Transfer time (days)")
clim([0 30]);  % Set your desired color range
set(gca, 'YDir', 'normal');
title("Neptune capture dV")

% Add contour lines
hold on;
levels = [0:0.5:5 6:2:30];  % Contour levels every 5 units between 0 and 50
[C, h] = contour(x, y, big_list_capture, levels, 'LineColor', 'k');
clabel(C, h);  % Add labels to contours
hold off;


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
function [v1, visviva] = getVisViva(r, v, p1, p2, tf_list, mu, t_elapsed)
% r, v are the planet ephimeris
% p1 is the planet number for the first planet (Venus - Neptune)
% p2 is the planet number for the second planet (Venus - Neptune)
% tf_list is an array of the flight times to check
    % iterate lambert's problem to get graphs
    visviva = zeros(length(tf_list), 2);
    v1 = zeros(length(tf_list), 3);
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
        v1(i, :) = V1;
        i = i + 1;
    end
end