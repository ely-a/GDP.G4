% Spacecraft and Planet Animation using Exported Trajectory and Planet Positions

close all
clc

% --- Load data from Python ---
sc_data = readmatrix('spacecraft_trajectory.txt'); % skip header if present

sc_time = sc_data(:,1);

planet_names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"];
colors = [
    0.70, 0.70, 0.70;  % Mercury (gray)
    0.95, 0.75, 0.10;  % Venus
    0.00, 0.45, 0.74;  % Earth
    0.85, 0.33, 0.10;  % Mars
    0.80, 0.60, 0.40;  % Jupiter
    0.93, 0.75, 0.50;  % Saturn
    0.40, 0.80, 0.90;  % Uranus
    0.30, 0.40, 0.85;  % Neptune
];

% Extract planet positions (each planet: 3 columns)
planet_xyz = cell(1,8);
for i = 1:8
    idx = 2 + (i-1)*3;
    planet_xyz{i} = sc_data(:, idx:idx+2); % Nx3 matrix for planet i
end

% Extract spacecraft position (last 3 columns)
sc_x = sc_data(:, end-2);
sc_y = sc_data(:, end-1);
sc_z = sc_data(:, end);

% --- Plot and animate ---
figure('Position', [100, 100, 1600, 900]);
hold on;
axis equal;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Spacecraft Trajectory and Planet Orbits');

% Plot full planet orbits (for the whole time span)
for i = 1:length(planet_names)
    plot3(planet_xyz{i}(:,1), planet_xyz{i}(:,2), planet_xyz{i}(:,3), ...
        'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', planet_names(i));
end

% Plot initial spacecraft position
sc_traj = plot3(sc_x(1), sc_y(1), sc_z(1), 'k-', 'LineWidth', 2, 'DisplayName', 'Spacecraft');
sc_dot = plot3(sc_x(1), sc_y(1), sc_z(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% legend show;

% --- Video Writer setup ---
v = VideoWriter('spacecraft_animation.mp4', 'MPEG-4');
v.FrameRate = 30;
open(v);

% Animation loop
for k = 1:10:length(sc_time)
    % Update spacecraft trail and dot
    set(sc_traj, 'XData', sc_x(1:k), 'YData', sc_y(1:k), 'ZData', sc_z(1:k));
    set(sc_dot, 'XData', sc_x(k), 'YData', sc_y(k), 'ZData', sc_z(k));
    % Optionally, update planet positions with moving dots
    for i = 1:length(planet_names)
        if k == 1
            planet_dot(i) = plot3(planet_xyz{i}(k,1), planet_xyz{i}(k,2), planet_xyz{i}(k,3), 'o', ...
                'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
        else
            set(planet_dot(i), 'XData', planet_xyz{i}(k,1), 'YData', planet_xyz{i}(k,2), 'ZData', planet_xyz{i}(k,3));
        end
    end
    drawnow;
    axis equal
    grid on
    xlim([-5e9 5e9]);
    ylim([-5e9 5e9]);
    zlim([-1e9 1e9]);
    frame = getframe(gcf);
    writeVideo(v, frame);
    pause(0.01)
end

close(v);
disp('Video saved as spacecraft_animation.mp4');

% --- Calculate minimum distance from each planet ---
min_dist = zeros(1,8);
min_epoch = zeros(1,8);

for i = 1:8
    % Compute distance at each time step
    d = sqrt(sum((planet_xyz{i} - [sc_x sc_y sc_z]).^2, 2));
    [min_dist(i), idx] = min(d);
    min_epoch(i) = sc_time(idx);
end

% --- Output results ---
fprintf('\nMinimum distance from spacecraft to each planet:\n');
for i = 1:8
    fprintf('%s: %.3e km at epoch (MJD2000) %.2f\n', planet_names(i), min_dist(i), min_epoch(i));
end

% --- Plot positions at closest approach to Jupiter ---
[~, idx_jup] = min(sqrt(sum((planet_xyz{5} - [sc_x sc_y sc_z]).^2, 2))); % Jupiter is planet 5
epoch_jup = sc_time(idx_jup);

figure;
hold on;
axis equal;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title(sprintf('Positions at Closest Approach to Jupiter (MJD2000 = %.2f)', epoch_jup));

% Plot all planets at this epoch
for i = 1:length(planet_names)
    plot3(planet_xyz{i}(idx_jup,1), planet_xyz{i}(idx_jup,2), planet_xyz{i}(idx_jup,3), 'o', ...
        'MarkerSize', 10, 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
        'DisplayName', planet_names(i));
end

% Plot spacecraft at this epoch
plot3(sc_x(idx_jup), sc_y(idx_jup), sc_z(idx_jup), 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'DisplayName', 'Spacecraft');

legend show;