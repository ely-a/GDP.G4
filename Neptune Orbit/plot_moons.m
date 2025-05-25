clear
clc 
close all

if exist('moon_orbits.mat', 'file') ~= 2
    load_moon_data()
end

load('moon_orbits.mat', 'orbits_mat', 'moon_names', 'max_len');

    num_moons = length(moon_names);
    trail_length = 100;
    colors = lines(num_moons);

    % Set up figure
    figure;
    hold on;
    axis equal;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    xlim([-0.5e6, 0.5e6]); ylim([-0.5e6, 0.5e6]); zlim([-0.3e6, 0.3e6]);
    grid on; view(3);
    title('Neptune and its Moons');

    % Neptune
    plot3(0, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    text(0, 0, 0, ' Neptune', 'Color', 'b', 'FontSize', 10);

    % Initialize plot handles
    moon_markers = gobjects(1, num_moons);
    moon_trails = gobjects(1, num_moons);
    for i = 1:num_moons
        moon_markers(i) = plot3(NaN, NaN, NaN, 'o', 'MarkerFaceColor', colors(i,:), 'DisplayName', moon_names{i});
        moon_trails(i) = plot3(NaN, NaN, NaN, '-', 'Color', colors(i,:));
    end
    legend;

    % Animate
    for frame = 1:max_len
        for i = 1:num_moons
            pos = orbits_mat(frame, :, i);
            set(moon_markers(i), 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));

            start = max(1, frame - trail_length + 1);
            trail = orbits_mat(start:frame, :, i);
            set(moon_trails(i), 'XData', trail(:,1), 'YData', trail(:,2), 'ZData', trail(:,3));
        end
        drawnow;
        pause(0.01);
    end