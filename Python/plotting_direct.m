close all; clear; clc;

% Settings
files = {'spacecraft_trajectory_direct_2031.txt', 'spacecraft_trajectory_direct_2032.txt', ... 
'spacecraft_trajectory_direct_2033.txt'};
labels = {'01/06/2031', '05/06/2032', '11/06/2033'};
colors = [
    0.0, 0.45, 0.74;   % Launch 2031 - blue
    0.85, 0.33, 0.10;  % Launch 2032 - orange/red
    0.47, 0.67, 0.19;  % Launch 2033 - olive green
];
planet_names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"];
planet_colors = turbo(8);

% Create figure
figure('Position', [100 100 1400 800]);
hold on; grid on; axis equal;
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
set(gca, 'FontSize', 14);
axis([-2e9 4.5e9 -1e9 3e9 -1e9 1e9]);

% Data storage
planet_trails = cell(1, 8);
planet_final_positions = zeros(8, 3);

traj_handles = gobjects(1, 3);  % to store spacecraft trajectory line handles


% Loop through each file
for i = 1:3
    data = readmatrix(files{i});
    time = data(:, 1);
    sc_xyz = data(:, end-2:end);
    
    % Plot spacecraft trajectory
    traj_handles(i) = plot3(sc_xyz(:,1), sc_xyz(:,2), sc_xyz(:,3), '-', ...
    'Color', colors(i,:), 'LineWidth', 2);
    
    % On first file, store planet trails and mark arrival epoch positions
    if i == 1
        for p = 1:8
            col_idx = 2 + (p-1)*3;
            xyz = data(:, col_idx:col_idx+2);
            planet_trails{p} = xyz;
        end
        arrival_index = size(data, 1);
        for p = 1:8
            planet_final_positions(p, :) = planet_trails{p}(arrival_index, :);
        end
    end
end

% Plot planet trails and arrival positions
for p = 1:8
    plot3(planet_trails{p}(:,1), planet_trails{p}(:,2), planet_trails{p}(:,3), ...
        'Color', planet_colors(p,:), 'LineStyle', '--', 'LineWidth', 1.2, ...
        'DisplayName', planet_names(p) + " Trail");
    
    % Mark final position
    plot3(planet_final_positions(p,1), planet_final_positions(p,2), planet_final_positions(p,3), ...
        'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', planet_colors(p,:), ...
        'DisplayName', planet_names(p) + " @ Arrival");
end

legend(traj_handles, labels, 'Location', 'southeast', 'FontSize', 18);
