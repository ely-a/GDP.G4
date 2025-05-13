% REQUIRED ADD-ONS
% type "aeroDataPackage" into the command window and get the planetary data

t_start = juliandate(2030,5,12);
t_end = t_start + 1e3;
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

xrange = [min(r(:, 7)), max(r(:, 7))]; % this will be for mars's orbit
yrange = [min(r(:, 8)), max(r(:, 8))];
zrange = [min(r(:, 9)), max(r(:, 9))];
% xrange = [min(r(:, 19)), max(r(:, 19))]; % this will be for neptune's orbit
% yrange = [min(r(:, 20)), max(r(:, 20))];
% zrange = [min(r(:, 21)), max(r(:, 21))];



% iterate lambert's problem to get graphs
tf_list = (10:400)'; % times of flight to check;
mu = 1.327e11;
visviva = getVisViva(r, v, 2, 1, tf_list, mu);

figure
semilogy(tf_list, visviva(:, 1))
hold on
semilogy(tf_list, visviva(:, 2))
grid on
xlabel("Time of Flight (days)")
ylabel("Vis-viva at intercept (km^2/s^2)")
legend("Departing Earth", "Arriving at Venus")

% % Lambert's problem from Earth at t=0 to Venus at t=100 days
% tf = 100; % days
% r1 = r(1, 4:6);
% r2 = r(1+tf, 1:3);
% mu = 1.327e11;
% [V1, V2] = lambert2(r1, r2, "retro", tf * 86400, mu);
% r_out = propagate_orbit(r1, V1, tf, mu);
% % pad to make the same length as the ephimeris data
% r(1:tf+1, 22:24) = r_out;
% r(tf+2:end, 22:24) = repmat(r_out(end, :), t_end-t_start-tf, 1);
%
%
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
%     0.30, 0.40, 0.85;  % Neptune
%     0, 0, 0 % s/c
% ];
% 
% nSteps = length(t_query);
% nPlanets = 8;
% trail = gobjects(nPlanets,1);
% dot = gobjects(nPlanets,1);
% 
% for i = [1 2 3 8]
%     trail(i) = plot3(NaN, NaN, NaN, 'Color', planet_colors(i,:), 'LineWidth', 1.5);
%     dot(i) = plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 6, ...
%                    'MarkerFaceColor', planet_colors(i,:), 'MarkerEdgeColor', 'k');
% end
% 
% % Animation loop
% for k = 1:10:nSteps  % Skip some steps to speed up animation
%     for i = [1 2 3 8]
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
%     pause(1e-1)
% end

% sequence
% [V E J N]
sequence = [1 2 4 7];
initialTOF = 140; %days
for flyby = 1:length(sequence)
    
end


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
function visviva = getVisViva(r, v, p1, p2, tf_list, mu)
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
        r1 = r(1, idx1);
        r2 = r(1+tf, idx2);
        [V1, V2] = lambert2(r1, r2, "pro", tf * 86400, mu);
        [V3, V4] = lambert2(r1, r2, "retro", tf * 86400, mu);
        % change to relative velocities
        VV1 = norm(V1 - v(1, idx1))^2;
        VV2 = norm(V2 - v(1+tf, idx2))^2;
        VV3 = norm(V3 - v(1, idx1))^2; % also try going the other way round the sun to see if that's better
        VV4 = norm(V4 - v(1+tf, idx2))^2;
        visviva(i, :) = [min(VV1, VV3) min(VV2, VV4)];
        i = i + 1;
    end
end