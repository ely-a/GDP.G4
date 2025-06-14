%% Apogee Raising Maneuver Optimization
% This script designs a multi-impulse orbit-raising strategy to move from
% a low Earth orbit (LEO) to a highly elliptical orbit (HEO).
%
% Problem Statement:
% - Initial Orbit: Circular, 300 km altitude.
% - Final Orbit: Eccentricity e = 0.96.
% - Constraint 1: Perigee altitude must remain constant at 300 km.
% - Constraint 2: Total maneuver time must be less than 3 weeks (21 days).
% - Objective: Optimize for minimal delta-v, which in this impulsive
%   model means finding the maximum number of burns possible within the
%   time limit.

clear; clc; close all;

%% 1. Define Constants and Parameters
% -------------------------------------------------------------------------
% Gravitational and Time Constants
mu = 3.986004418e5; % Earth's standard gravitational parameter [km^3/s^2]
R_earth = 6378;   % Earth's equatorial radius [km]
T_max_days = 15;      % Maximum allowed maneuver duration [days]
T_max_sec = T_max_days * 24 * 3600; % Convert max time to seconds

% Initial and Final Orbit Parameters
h_perigee = 300;    % Constant perigee altitude [km]
e_final = 0.96;       % Desired final eccentricity

%% 2. Initial and Final Orbit Calculations
% -------------------------------------------------------------------------
% The perigee radius (closest point) is constant for all orbits.
r_p = R_earth + h_perigee; % Perigee radius [km]

% --- Initial Circular Orbit ---
% For a circle, radius is the semi-major axis, and eccentricity is 0.
r_initial = r_p;
a_initial = r_initial;
e_initial = 0;

% Calculate the velocity in the initial circular orbit.
v_initial_circ = sqrt(mu / r_initial);

% --- Final Highly Elliptical Orbit ---
% Calculate the semi-major axis of the final orbit using the formula for
% perigee radius: r_p = a * (1 - e)
a_final = r_p / (1 - e_final);

% Calculate the velocity required at perigee for the final orbit using the
% Vis-viva equation: v^2 = mu * (2/r - 1/a)
v_p_final = sqrt(mu * (2/r_p - 1/a_final));

% --- Total Delta-V Calculation ---
% In an impulsive model, the total delta-v is the final perigee velocity
% minus the initial perigee velocity, regardless of the number of burns.
total_delta_v = v_p_final - v_initial_circ;

fprintf('--- ORBIT PARAMETERS ---\n');
fprintf('Initial Circular Velocity: %.2f km/s\n', v_initial_circ );
fprintf('Final Perigee Velocity:    %.2f km/s\n', v_p_final );
fprintf('Total Impulsive Delta-V:   %.2f km/s (%.2f m/s)\n\n', total_delta_v , total_delta_v);


%% 3. Find Optimal Number of Burns (N) for Time Constraint
% -------------------------------------------------------------------------
% We will iterate, increasing the number of burns (N), and calculate the
% total time for each N. The optimal N is the largest value for which the
% total maneuver time is still under the 21-day limit.

N = 1; % Start with one burn
time_history = [];
n_history = [];
optimal_N = 0; % Initialize in case the first orbit is already too long
final_time_days = 0;

fprintf('--- FINDING OPTIMAL NUMBER OF BURNS ---\n');
fprintf('N \t Total Time (days)\n');
fprintf('--------------------------\n');

while true
    % We assume the change in specific energy is constant for each burn.
    % The specific orbital energy is -mu/(2a), so we make the change in 1/a
    % constant for each of the N steps.
    inv_a_initial = 1 / a_initial;
    inv_a_final = 1 / a_final;
    delta_inv_a = (inv_a_final - inv_a_initial) / N;

    % Calculate the semi-major axis for each intermediate transfer orbit
    a_transfer = zeros(1, N);
    for i = 1:N
        inv_a_i = inv_a_initial + i * delta_inv_a;
        a_transfer(i) = 1 / inv_a_i;
    end

    % The total time is the sum of the periods of the N-1 transfer orbits.
    % After the N-th burn, we are in the final orbit and the maneuver is over.
    total_time_sec = 0;
    if N > 1
        for i = 1:(N-1)
            % Period of an orbit T = 2*pi*sqrt(a^3/mu)
            T_i = 2 * pi * sqrt(a_transfer(i)^3 / mu);
            total_time_sec = total_time_sec + T_i;
        end
    end
    
    total_time_days = total_time_sec / (24 * 3600);
    fprintf('%d \t %.4f\n', N, total_time_days);

    % Store history for plotting
    time_history(end+1) = total_time_days;
    n_history(end+1) = N;

    % Check if the calculated time exceeds our maximum allowed time
    if total_time_days > T_max_days
        % If it exceeds, the previous N was the last valid one.
        optimal_N = N - 1;
        if optimal_N > 0
          final_time_days = time_history(end-1);
        end
        break;
    end

    % If time is somehow zero or negative, break to prevent infinite loop
    if total_time_days <= 0 && N > 1
        warning('Time calculation resulted in non-positive value. Breaking loop.');
        optimal_N = N-1;
        final_time_days = time_history(end-1);
        break;
    end
    
    N = N + 1;
end

fprintf('\n--- OPTIMIZATION RESULTS ---\n');
if optimal_N > 0
    fprintf('Optimal Number of Burns (N): %d\n', optimal_N);
    fprintf('Total Maneuver Time:         %.2f days (Limit: %d days)\n', final_time_days, T_max_days);
else
    fprintf('No solution found. Even a single transfer orbit exceeds the time limit.\n');
end
fprintf('Total Delta-V Required:      %.2f km/s\n', total_delta_v);


%% 4. Plotting the Results
% -------------------------------------------------------------------------

% Plot 1: Maneuver Time vs. Number of Burns
figure('Name', 'Maneuver Time vs. Number of Burns');
plot(n_history, time_history, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
yline(T_max_days, 'r--', 'LineWidth', 2, 'Label', '21-Day Time Limit');
if optimal_N > 0
    plot(optimal_N, final_time_days, 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'LineWidth', 2);
    text(optimal_N, final_time_days, sprintf('  Optimal Point\n  N = %d', optimal_N), 'VerticalAlignment', 'bottom');
end
grid on;
title('Total Maneuver Time vs. Number of Burns');
xlabel('Number of Burns (N)');
ylabel('Total Time (days)');
legend('Maneuver Time', 'Time Limit', 'Optimal N', 'Location', 'northwest');
set(gca, 'FontSize', 12);

% Plot 2: Visualization of the Orbits
figure('Name', 'Orbit Trajectory Visualization');
hold on;
axis equal;

% Plot Earth
[x_e, y_e, z_e] = sphere;
surf(x_e*R_earth, y_e*R_earth, z_e*R_earth, 'FaceColor', 'blue', 'EdgeColor', 'none');

% Define angles for plotting once to be efficient
theta = linspace(0, 2*pi, 500);

% SYNTAX FIX: Anonymous functions can only contain a single expression.
% This function now calculates r and converts to x,y in one line.
plot_orbit = @(a, e, style) plot( ...
    ((a * (1 - e^2)) ./ (1 + e * cos(theta))) .* cos(theta) , ...
    ((a * (1 - e^2)) ./ (1 + e * cos(theta))) .* sin(theta) , ...
    style, 'LineWidth', 0.5);

% Plot initial and final orbits
h_initial_plot = plot_orbit(a_initial, e_initial, 'g--');
h_final_plot = plot_orbit(a_final, e_final, 'r-');

% MODIFIED: Plot ALL intermediate orbits from the optimal solution
h_intermediate_plots = [];
if optimal_N > 1
    % Calculate the parameters needed for the optimal N transfer orbits
    inv_a_initial_opt = 1 / a_initial;
    inv_a_final_opt = 1 / a_final;
    delta_inv_a_opt = (inv_a_final_opt - inv_a_initial_opt) / optimal_N;

    % Loop through and plot ALL intermediate transfer orbits
    fprintf('\nPlotting all %d intermediate transfer orbits...\n', optimal_N - 1);
    for i = 1:(optimal_N - 1) % We plot N-1 transfer orbits
        inv_a_i_opt = inv_a_initial_opt + i * delta_inv_a_opt;
        a_i = 1 / inv_a_i_opt;
        e_i = 1 - r_p / a_i; % eccentricity of transfer orbit, r_p is constant
        
        % Plot the orbit. Store the handle only for the first one for the legend.
        plot_handle = plot_orbit(a_i, e_i, 'k-');
        if i == 1
            h_intermediate_plots = plot_handle;
        end
    end
end


title('Orbit Trajectory Visualization');
xlabel('X [km]');
ylabel('Y [km]');

% LOGIC FIX: Build legend dynamically to prevent errors if no intermediate
% orbits are plotted.
legend_handles = [h_initial_plot, h_final_plot];
legend_labels = {'Initial LEO', 'Final HEO'};
if ~isempty(h_intermediate_plots)
    legend_handles = [legend_handles(1), h_intermediate_plots(1), legend_handles(2)];
    legend_labels = {'Initial LEO', 'Intermediate Orbits', 'Final HEO'};
end
legend(legend_handles, legend_labels, 'Location', 'northwest');

grid on;
set(gca, 'FontSize', 12);
view(2); % 2D view