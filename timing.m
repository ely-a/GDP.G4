t_start = juliandate(2025,5,12);
t_end = t_start + 8e4;
t_query = (t_start:10:t_end)';
i = 1;
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

xrange = [min(r(:, 19)), max(r(:, 19))]; % this will be for neptune's orbit
yrange = [min(r(:, 20)), max(r(:, 20))];
zrange = [min(r(:, 21)), max(r(:, 21))];

% % Initialize the figure
% figure;
% axis equal;
% 
% xlabel('X (km)');
% ylabel('Y (km)');
% zlabel('Z (km)');
% view(3);
% hold on;
% 
% % Define planet colors
% planet_colors = [
%     0.95, 0.75, 0.10;  % Venus
%     0.00, 0.45, 0.74;  % Earth
%     0.85, 0.33, 0.10;  % Mars
%     0.80, 0.60, 0.40;  % Jupiter
%     0.93, 0.75, 0.50;  % Saturn
%     0.40, 0.80, 0.90;  % Uranus
%     0.30, 0.40, 0.85   % Neptune
% ];
% 
% nSteps = length(t_query);
% nPlanets = 7;
% trail = gobjects(nPlanets,1);
% dot = gobjects(nPlanets,1);
% 
% for i = 1:nPlanets
%     trail(i) = plot3(NaN, NaN, NaN, 'Color', planet_colors(i,:), 'LineWidth', 1.5);
%     dot(i) = plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 6, ...
%                    'MarkerFaceColor', planet_colors(i,:), 'MarkerEdgeColor', 'k');
% end
% 
% % Animation loop
% for k = 1:10:nSteps  % Skip some steps to speed up animation
%     for i = 1:nPlanets
%         idx = (i-1)*3 + 1;
%         % Update trails
%         set(trail(i), 'XData', r(1:k, idx), ...
%                       'YData', r(1:k, idx+1), ...
%                       'ZData', r(1:k, idx+2));
%         % Update dots
%         set(dot(i), 'XData', r(k, idx), ...
%                     'YData', r(k, idx+1), ...
%                     'ZData', r(k, idx+2));
%     end
%     drawnow;
%     axis equal
%     grid on
%     xlim(xrange);
%     ylim(yrange);
%     zlim(zrange);
%     pause(1e-6)
% end