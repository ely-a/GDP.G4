clear all
close all
clc

% atmospheric model
R_Neptune = 24622000; % m
heights = linspace(0, 2000, 100) * 1e3;
% density vs alt.
rho_vals = [3.31533855e-01 1.52964747e-01 2.54401400e-01 1.31939086e-02
 2.48900219e-03 1.22881657e-03 4.82506482e-04 3.00252442e-04
 1.98527623e-04 1.06203031e-04 4.96398560e-05 3.04623928e-05
 2.31915279e-05 1.75300299e-05 1.34847500e-05 1.02300151e-05
 7.85979220e-06 4.76112300e-06 1.60950582e-06 1.50608165e-06
 1.53042080e-06 1.51820245e-06 1.16067723e-06 8.89418675e-07
 2.15200431e-07 1.11560507e-07 7.49421287e-08 4.89378619e-08
 3.20027160e-08 2.32021889e-08 1.76842078e-08 1.33148516e-08
 1.00341358e-08 7.67965586e-09 5.81665459e-09 4.41008506e-09
 3.35839284e-09 2.56241290e-09 1.93778878e-09 1.46061480e-09
 1.10186120e-09 8.40054036e-10 6.49919044e-10 5.12139214e-10
 4.10170536e-10 3.33308595e-10 2.82138470e-10 2.37226314e-10
 1.98042650e-10 1.64119929e-10 1.35590988e-10 1.11458298e-10
 9.42071848e-11 7.96684703e-11 6.76206403e-11 5.79172149e-11
 4.98117293e-11 4.28649025e-11 3.75397747e-11 3.33176310e-11
 2.96837623e-11 2.63495741e-11 2.33715817e-11 2.07352438e-11
 1.85502752e-11 1.69721238e-11 1.55838521e-11 1.43115557e-11
 1.30970067e-11 1.15762097e-11 1.02977605e-11 9.30877656e-12
 8.42905143e-12 7.63646306e-12 6.91907207e-12 6.26967004e-12
 5.67248686e-12 5.17788096e-12 4.72671932e-12 4.31515381e-12
 3.93968049e-12 3.58638787e-12 3.25119240e-12 2.94750818e-12
 2.67235335e-12 2.42303033e-12 2.23364266e-12 2.06082379e-12
 1.90142840e-12 1.75440929e-12 1.61880139e-12 1.50434796e-12
 1.40566585e-12 1.31349140e-12 1.22739290e-12 1.14696747e-12
 1.06492368e-12 9.75065638e-13 8.95724624e-13 8.22854111e-13];
rho_vals = reshape(rho_vals.', [100, 1]);
rho_vals_min = 10 .^ (log10(rho_vals) - 0.1);
rho_vals_max = 10 .^ (log10(rho_vals) + 0.1);

% inputs
mu = 6.836529e15; % Neptune gravitational parameter, km^3/s^2
h_target = 2e6; % target altitude in m
v_pe = 29000;  % m/s

% orbiter properties
beta = 895;

% simulation setup
dt = 5; % seconds

% exit condition lists
aimList = []; % over the aiming periapses which exit the atmosphere
apoList = [];
maxaccList = [];
aimList_min = [];
apoList_min = [];
maxaccList_min = [];
aimList_max = [];
apoList_max = [];
maxaccList_max = [];

% plotting setup
figure(1); hold on; title("Acceleration vs Time"); xlabel("Time (s)"); ylabel("Acceleration (m/s^2)");
figure(2); hold on; title("Velocity vs Time"); xlabel("Time (s)"); ylabel("Velocity (m/s)");
figure(3); hold on; title("Altitude vs Time"); xlabel("Time (s)"); ylabel("Altitude (m)");
figure(4); hold on; title("Trajectory"); xlabel("Distance (km)"); ylabel("Distance (km)");
plot_h = 177500; % periapsis aiming height in m to plot the trajectory for

for h_pe = plot_h:plot_h %1.3e5:1e2:2e5

    disp(strcat("Aiming Periapsis: ", num2str(h_pe/1e3), " km"))
    % plotting lists
    aList = [];
    vList = [];
    altList = [];
    trajList = [];
    tList = [];

    % calculations
    r_pe = R_Neptune + h_pe;
    r_target = R_Neptune + h_target;
    energy = v_pe^2 / 2 - mu / r_pe;
    h = r_pe * v_pe;
    v_target = sqrt(2*(energy + mu / r_target));
    
    % entry conditions
    v_tan = h / r_target;
    v_rad = sqrt(v_target^2 - v_tan^2);

    % simulation setup
    t = 0; % seconds
    r = [0; h_target + R_Neptune]; % m, [x, y]
    v = [v_tan; -v_rad];
    
    while norm(r) <= h_target + R_Neptune && norm(r) > R_Neptune...
            && t < 1e3 % check if safe
    
        a_grav = -mu * r ./ norm(r)^3;
    
        % Interpolate atmospheric density
        h = norm(r) - R_Neptune;
        rho = interp1(heights, rho_vals, h);
    
        % drag
        a_drag = -0.5 * rho * v * norm(v) / beta;
    
        % Total acceleration
        a = a_grav + a_drag;
    
        % Update state
        v = v + a * dt;
        r = r + v * dt;
    
        % Save for plotting
        aList = [aList norm(a_drag)];
        vList = [vList norm(v)];
        altList = [altList norm(r) - R_Neptune];
        trajList = [trajList r/1000];
        tList = [tList t];
        t = t + dt;
    end

    labelStr = sprintf('%.0f km', h_pe / 1e3);
    
    if h_pe == plot_h
        % Plot acceleration
        figure(1)
        plot(tList, aList, LineWidth=2, color="black");
        text(tList(end), aList(end), labelStr, 'FontSize', 8);
        grid on
    
        % Plot velocity
        figure(2)
        plot(tList, vList, LineWidth=2, color="black");
        text(tList(end), vList(end), labelStr, 'FontSize', 8);
        grid on
    
        % Plot altitude
        figure(3)
        plot(tList, altList, LineWidth=2, color="black");
        text(tList(end), altList(end), labelStr, 'FontSize', 8);
        grid on
    end

    % exit conditions
    if norm(r) > R_Neptune && dot(v, r) > 0 % in this case we're leaving again
        % orbital stuff
        energy = norm(v)^2 / 2 - mu / norm(r);
        sma = -mu / (2 * energy);
        h_vec = cross([r; 0], [v; 0]);
        e = norm(cross([v; 0], h_vec)) / mu;
        r_a = sma * (1 + e);
        if r_a > 0
            aimList = [aimList h_pe / 1e3]; % aiming periapsis
            apoList = [apoList r_a / 1e3]; % exit apoapsis
            % other
            maxacc = max(aList);
            maxaccList = [maxaccList maxacc];
        end
    end

    % ====================== MIN DENSITY CASE ===========================

    % calculations
    r_pe = R_Neptune + h_pe;
    r_target = R_Neptune + h_target;
    energy = v_pe^2 / 2 - mu / r_pe;
    h = r_pe * v_pe;
    v_target = sqrt(2*(energy + mu / r_target));
    
    % entry conditions
    v_tan = h / r_target;
    v_rad = sqrt(v_target^2 - v_tan^2);

    % simulation setup
    t = 0; % seconds
    r = [0; h_target + R_Neptune]; % m, [x, y]
    v = [v_tan; -v_rad];
    
    while norm(r) <= h_target + R_Neptune && norm(r) > R_Neptune...
            && t < 1e3 % check if safe
    
        a_grav = -mu * r ./ norm(r)^3;
    
        % Interpolate atmospheric density
        h = norm(r) - R_Neptune;
        rho = interp1(heights, rho_vals_min, h);
    
        % drag
        a_drag = -0.5 * rho * v * norm(v) / beta;
    
        % Total acceleration
        a = a_grav + a_drag;
    
        % Update state
        v = v + a * dt;
        r = r + v * dt;
        t = t + dt;
    end

    labelStr = sprintf('%.0f km', h_pe / 1e3);

    % exit conditions
    if norm(r) > R_Neptune && dot(v, r) > 0 % in this case we're leaving again
        % orbital stuff
        energy = norm(v)^2 / 2 - mu / norm(r);
        sma = -mu / (2 * energy);
        h_vec = cross([r; 0], [v; 0]);
        e = norm(cross([v; 0], h_vec)) / mu;
        r_a = sma * (1 + e);
        if r_a > 0
            aimList_min = [aimList_min h_pe / 1e3]; % aiming periapsis
            apoList_min = [apoList_min r_a / 1e3]; % exit apoapsis
            % other
            maxacc = max(aList);
            maxaccList_min = [maxaccList_min maxacc];
        end
    end

    % ===================== MAX DENSITY CASE ========================

    % calculations
    r_pe = R_Neptune + h_pe;
    r_target = R_Neptune + h_target;
    energy = v_pe^2 / 2 - mu / r_pe;
    h = r_pe * v_pe;
    v_target = sqrt(2*(energy + mu / r_target));
    
    % entry conditions
    v_tan = h / r_target;
    v_rad = sqrt(v_target^2 - v_tan^2);

    % simulation setup
    t = 0; % seconds
    r = [0; h_target + R_Neptune]; % m, [x, y]
    v = [v_tan; -v_rad];
    
    while norm(r) <= h_target + R_Neptune && norm(r) > R_Neptune...
            && t < 1e3 % check if safe
    
        a_grav = -mu * r ./ norm(r)^3;
    
        % Interpolate atmospheric density
        h = norm(r) - R_Neptune;
        rho = interp1(heights, rho_vals_max, h);
    
        % drag
        a_drag = -0.5 * rho * v * norm(v) / beta;
    
        % Total acceleration
        a = a_grav + a_drag;
    
        % Update state
        v = v + a * dt;
        r = r + v * dt;
        t = t + dt;
    end

    labelStr = sprintf('%.0f km', h_pe / 1e3);

    % exit conditions
    if norm(r) > R_Neptune && dot(v, r) > 0 % in this case we're leaving again
        % orbital stuff
        energy = norm(v)^2 / 2 - mu / norm(r);
        sma = -mu / (2 * energy);
        h_vec = cross([r; 0], [v; 0]);
        e = norm(cross([v; 0], h_vec)) / mu;
        r_a = sma * (1 + e);
        if r_a > 0
            aimList_max = [aimList_max h_pe / 1e3]; % aiming periapsis
            apoList_max = [apoList_max r_a / 1e3]; % exit apoapsis
            % other
            maxacc = max(aList);
            maxaccList_max = [maxaccList_max maxacc];
        end
        finalV = v;
        finalR = r;
    end
end

% plot Neptune
figure(4)
theta = linspace(0, 2*pi, 500);
circleX = R_Neptune * cos(theta) / 1000;
circleY = R_Neptune * sin(theta) / 1000;
plot(circleX, circleY)
fill(circleX, circleY, "blue", "FaceAlpha", 0.5);
hold on
% plot pre-aerocapture section
initialV = [-v_tan; v_rad] / 1000;
initialR = [0; h_target + R_Neptune] / 1000;
tf = 1e4;
r_out = propagate_orbit(initialR, initialV, tf, mu/1e9);
plot(r_out(:, 1), r_out(:, 2), color="black", LineWidth=2);
% plot aerocapture section
plot(trajList(1, :), trajList(2, :), color="black", LineWidth=2);
% plot post-aerocapture section
initialV = finalV / 1000;
initialR = trajList(:, end);
tf = 1e4;
r_out = propagate_orbit(initialR, initialV, tf, mu/1e9);
plot(r_out(:, 1), r_out(:, 2), color="black", LineWidth=2);
% formatting
xlabel("X distance (km)")
ylabel("Y distance (km)")
grid on
axis equal

figure(5)
semilogy(aimList, apoList, color="black", LineWidth=2)
hold on
semilogy(aimList_min, apoList_min, color="black", Linestyle="--")
semilogy(aimList_max, apoList_max, color="black", Linestyle="--")
xlabel("Aiming periapsis (km)")
ylabel("Exit apoapsis (km)")
grid on
xline(max(aimList_min), LineStyle=":", label="Max. aiming periapsis", LineWidth=2)
xline(min(aimList_max), LineStyle=":", label="Min. aiming periapsis", LineWidth=2)

figure(6)
plot(aimList, maxaccList, color="black", LineWidth=2)
hold on
plot(aimList_min, maxaccList_min, color="black", LineStyle="--")
plot(aimList_max, maxaccList_max, color="black", LineStyle="--")
xlabel("Aiming periapsis (km)")
ylabel("Maximum acceleration (m/s^2)")
grid on
xline(max(aimList_min), LineStyle=":", label="Max. aiming periapsis", LineWidth=2)
xline(min(aimList_max), LineStyle=":", label="Min. aiming periapsis", LineWidth=2)


% extra functions
% for plotting only
function r_out = propagate_orbit(r1, v1, tf, mu)
    % Inputs:
    % r1 - 2x1 position vector (km)
    % v1 - 2x1 velocity vector (km/s)
    % tf - total propagation time (s)
    % m  - gravitational parameter μ (km^3/s^2)

    % Time vector from 0 to tf
    tspan = (0:1:tf);

    % Initial state vector [r; v]
    y0 = [r1(:); v1(:)];

    % Define two-body dynamics
    function dydt = two_body(~, y)
        r = y(1:2);
        v = y(3:4);
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