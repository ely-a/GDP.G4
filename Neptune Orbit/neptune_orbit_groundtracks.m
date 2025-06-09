clear
clc
close all

%% SET GRAPHING FONTS AND SIZES
set(groot,'defaultLineLineWidth',4) % Set all line widths to 2 unless stated otherwise
set(groot,'defaultAxesFontSize',25) % Set all axes font size to 16 unless stated otherwise
set(groot,'defaulttextfontsize',36) % Set all text font sizes to 16 unless stated otherwise
set(groot,'defaultLineMarkerSize',12) % Set all marker sizes to 16 unless stated otherwise
set(groot,'defaultAxesXGrid','on') % Turn x grid lines on
set(groot,'defaultAxesYGrid','on') % Turn y grid lines on

%% Ground Track Simulation Over Neptune

% % === ORBIT SETUP ===
% mu_N = 6835100;                % Neptune gravitational parameter [km^3/s^2]
% R_orbit = 30000;               % km radius (example circular orbit)
% inclination_deg = 45;         % orbital inclination in degrees
% i_rad = deg2rad(inclination_deg);
% 
% T_orbit = 2 * pi * sqrt(R_orbit^3 / mu_N);  % orbital period in seconds
% n_orbits = 3;                                % number of orbits to simulate
% t = linspace(0, n_orbits * T_orbit, 1000);   % simulate N orbits
% jd0 = 2451545.0;
% jd_array = jd0 + t / 86400;                  % convert time to Julian date
% 
% % Generate inclined ECI orbit
% r_ECI_all = zeros(3, length(t));
% for k = 1:length(t)
%     theta = 2*pi * t(k) / T_orbit;  % angle swept since epoch  ahhh
%     slight dff method for calc theta and reci below
%     x = R_orbit * cos(theta);
%     y = R_orbit * sin(theta) * cos(i_rad);
%     z = R_orbit * sin(theta) * sin(i_rad);
%     r_ECI_all(:,k) = [x; y; z];
% end


% === LOAD FINAL ORBITAL ELEMENTS ===
load('ScienceOrbit.mat') 

mu_N = 6835100;                          % Neptune gravitational parameter [km^3/s^2]
R_N = 24622;

rp_newnew = R_N + 1600; 
ra_newnew = 61098.5950141566;
a_final = (ra_newnew + rp_newnew)/2;
e_final = (ra_newnew - rp_newnew)/(ra_newnew + rp_newnew);
h_final = sqrt(rp_newnew * mu_N * (1 + e_final));
i_final = 65;

% === USER INPUT: INSERTION EPOCH AND TIME OFFSET TO FINAL ORBIT ===
epoch_insertion = 2466307.86031975;             % Julian date of insertion
t_since_insertion = 86400 * t_days;          % seconds from insertion to final orbit apoapsis [example: 10 days]


jd0 = epoch_insertion + t_since_insertion / 86400;  % Start date of science orbit
T_orbit = 2 * pi * sqrt(a_final^3 / mu_N);           % Orbital period
n_orbits = 3;                                         

% === SIMULATION TIME SETUP ===
t = linspace(0, n_orbits * T_orbit, 1500);  % seconds
jd_array = jd0 + t / 86400;                 % Julian dates

% === CALCULATE r_ECI_all USING ORBITAL ELEMENTS ===
r_ECI_all = zeros(3, length(t));
v_ECI_all = zeros(3,length(t));
for k = 1:length(t)
    M = 2 * pi * t(k) / T_orbit;            % Mean anomaly
    E = kepler_E(M, e_final);              % Solve Kepler's Equation
    theta = 2 * atan2(sqrt(1 + e_final) * sin(E/2), sqrt(1 - e_final) * cos(E/2));  % True anomaly
    [r_vec, v_vec] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, rad2deg(theta), mu_N);
    r_ECI_all(:,k) = r_vec;
    v_ECI_all(:,k) = v_vec;
end

% === CALL FUNCTION TO PLOT GROUND TRACKS ===
%plot_ground_track_neptune_modes(r_ECI_all, jd_array, 'neptune_map.png', n_orbits);

% == PLANE CHANGE ==
i_target_deg = 80;

[r_ECI_inclined, v_ECI_inclined, jd_array_inclined] = rvj_from_pc( ...
    r_ECI_all(:,end), v_ECI_all(:,end), i_target_deg, jd_array(:,end), n_orbits); % needs to take in final positions from inital science

% === PLOT NEW GROUND TRACK ===
%plot_ground_track_neptune_modes(r_ECI_inclined, jd_array_inclined, 'neptune_map.png', n_orbits);

%% == Number of Science Orbits == 
% -- This excludes final de-orbit --
% -- Also its excluding aerobrake --

t_in_sc = 5*365 - t_days;
p_sc = 2 * pi * a_final^(3/2) / sqrt(mu_N);
p_sc_days = p_sc /(60*60*24);
n_orbits_in_sc = t_in_sc/p_sc_days;

fprintf('Number of Orbits in Science (excluding aerobrake) = %.3f \n', ...
        n_orbits_in_sc);

% -- Plotting Science Orbit (Post Aerobrake) -- 
t = linspace(0, p_sc, 1000);              % 1 full orbit

% === Generate Orbits ===
r_science = get_orbit(a_final, e_final, i_final, RAAN_final, omega_final, mu_N, t);
r_inclined = get_orbit(a_final, e_final, i_target_deg, RAAN_final, omega_final, mu_N, t);

% === Plotting ===
figure
hold on
axis equal
grid on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Science Orbit and Inclined Orbit around Neptune')

% Plot Neptune
[Xs, Ys, Zs] = sphere(100);
surf(R_N*Xs, R_N*Ys, R_N*Zs, 'FaceColor', [0 0.3 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Neptune');

% Plot Orbits
%plot3(r_science(1,:), r_science(2,:), r_science(3,:), 'k', 'DisplayName', 'Science Orbit', 'LineWidth', 2)
%plot3(r_inclined(1,:), r_inclined(2,:), r_inclined(3,:), 'm--', 'DisplayName', sprintf('Inclined Orbit (%d°)', i_target_deg), 'LineWidth', 2)

% === Magnetic L-shell Calculation Along Science Orbit ===
tilt_deg = 47;
tilt_rad = deg2rad(tilt_deg);
R_tilt = [1 0 0;
          0 cos(tilt_rad), -sin(tilt_rad);
          0 sin(tilt_rad),  cos(tilt_rad)];

% Rotate orbit into magnetic frame
r_magnetic = R_tilt' * r_science;
r_mag = vecnorm(r_magnetic);
z_m = r_magnetic(3, :);
mag_lat = asin(z_m ./ r_mag);
cos2_lat = cos(mag_lat).^2;
L_shell = r_mag ./ (R_N * cos2_lat);

% === Plot Neptune
[Xs, Ys, Zs] = sphere(100);
surf(R_N*Xs, R_N*Ys, R_N*Zs, ...
     'FaceColor', [0 0.3 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
     'DisplayName', 'Neptune');

% === Plot Science Orbit color-coded by L-shell ===
scatter3(r_science(1,:), r_science(2,:), r_science(3,:), ...
         15, L_shell, 'filled');  % color by L-shell
colormap('turbo'); 
cb = colorbar; 
cb.Label.String = 'L-shell Value';
caxis([1 2]);  % Adjust if your orbit goes higher



% === Plot Radiation Belt (L = 1.2 Shell) with Magnetic Tilt ===
L_thresh = 1.2;
R_rad = R_N * L_thresh;
tilt_deg = 47;
tilt_rad = deg2rad(tilt_deg);

% Generate unrotated sphere
[xs, ys, zs] = sphere(100);
xs = xs * R_rad;
ys = ys * R_rad;
zs = zs * R_rad;

% Flatten matrices to vectors
pts = [xs(:), ys(:), zs(:)]';

% Rotation about x-axis
R_tilt = [1, 0, 0;
          0, cos(tilt_rad), -sin(tilt_rad);
          0, sin(tilt_rad),  cos(tilt_rad)];

pts_rotated = R_tilt * pts;

% Reshape to 2D surface format
xs_rot = reshape(pts_rotated(1,:), size(xs));
ys_rot = reshape(pts_rotated(2,:), size(ys));
zs_rot = reshape(pts_rotated(3,:), size(zs));

% Plot the tilted radiation belt shell
surf(xs_rot, ys_rot, zs_rot, ...
     'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
     'FaceColor', [1 0 0], 'DisplayName', 'L = 1.2 Radiation Belt (Tilted)');


legend show
view(35, 25)


%% IMPORTANT NOTE ==

% We are at low altitude, but also at a higher magnetic latitude due to Neptune's 47° tilt
% 
% so That cos²(λ) term in the denominator reduces, pushing L slightly above
% 1.2 (from radiation equation)
% So:
% 
% Even though you're physically < 1.2 × Rₙ (as shown in figure) radially, you’re outside L = 1.2 due to magnetic geometry
%https://en.wikipedia.org/wiki/L-shell?utm_source=chatgpt.com

%% == End of Life == 

[r_before, v_before] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, 180, mu_N);
ra_end = a_final * (1 + e_final);
rp_end = R_N + 250;
e_end = (ra_end - rp_end)/(ra_end + rp_end);
h_end = sqrt(rp_end * mu_N * (1 + e_end));
v_end = h_end / ra_end;

deltaV_end = abs(norm(v_before) - v_end);
fprintf('End of Life ΔV = %.3f km/s\n', ...
        deltaV_end);

%% == Probe Release ==

rp_probe = R_N + 4000;
cos_theta = (h_final^2 / (mu_N * rp_probe) - 1) / e_final;
theta_probe = acosd(cos_theta);  % in degrees

[r_probe, v_probe] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, theta_probe, mu_N);

fprintf('Velocity at Probe Eject = %.3f km/s\n', ...
        norm(v_probe));

%% == Radiation belt == 
rad_tilt_deg = 47;      % Magnetic tilt [deg]

% Derived period
T_final = 2 * pi * sqrt(a_final^3 / mu_N);  % [s]

% === Discretise True Anomaly ===
theta_array = linspace(0, 360, 1000);  % degrees
dt_array = zeros(1, length(theta_array));  % placeholder for time

r_vecs = zeros(3, length(theta_array));  % position vectors

for k = 1:length(theta_array)
    theta_k = theta_array(k);
    [r_k, ~] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, theta_k, mu_N);
    r_vecs(:, k) = r_k;
end

r_mag = vecnorm(r_vecs);  % magnitude of radius vectors

% === Time Calculation from Kepler's Equation ===
theta_rad = deg2rad(theta_array);
E = 2 * atan(tan(theta_rad/2) .* sqrt((1 - e_final) / (1 + e_final)));
E = mod(E, 2*pi);
M = E - e_final .* sin(E);
t_array = (M / (2 * pi)) * T_final / 60;  % time in minutes

% === Magnetic Field Tilt and L-shell Calculation ===
tilt_rad = deg2rad(rad_tilt_deg);
R_tilt = [1 0 0;
          0 cos(tilt_rad) -sin(tilt_rad);
          0 sin(tilt_rad)  cos(tilt_rad)];

r_magnetic = R_tilt' * r_vecs;
z_m = r_magnetic(3, :);
mag_lat = asin(z_m ./ r_mag);
cos2_lat = cos(mag_lat).^2;

L_shell = r_mag ./ (R_N * cos2_lat);

%
% % === Check L-shell value specifically at periapsis ===
% [r_peri, ~] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, 0, mu_N);
% r_mag_peri = norm(r_peri);
% 
% % Rotate to magnetic frame
% tilt_rad = deg2rad(rad_tilt_deg);
% R_tilt = [1 0 0;
%           0 cos(tilt_rad) -sin(tilt_rad);
%           0 sin(tilt_rad)  cos(tilt_rad)];
% r_m_peri = R_tilt' * r_peri;
% lat_peri = asin(r_m_peri(3) / r_mag_peri);
% L_peri = r_mag_peri / (R_N * cos(lat_peri)^2);
% 
% fprintf('\n🔎 L-shell at Periapsis = %.4f\n', L_peri);
% %%
% % === Diagnostics: how many points are below L_thresh? ===
% num_below = sum(L_shell < L_thresh);
% fprintf('Samples with L < %.1f: %d out of %d\n', L_thresh, num_below, length(L_shell));

% === Radiation Analysis ===
L_thresh = 1.2;
in_radiation = L_shell < L_thresh;
time_in_radiation = trapz(t_array, double(in_radiation));  % minutes


% Output diagnostics
fprintf('Minimum L-shell value in orbit: %.4f\n', min(L_shell));

% UNCOMMMENTT
% === Plotting ===
% figure;
% plot(t_array, L_shell, 'b', 'LineWidth', 2); hold on;
% yline(L_thresh, 'r--', 'LineWidth', 1.5, 'Label', 'L = 1.2');
% area(t_array, L_shell .* double(in_radiation), ...
%     'FaceColor', [1 0.6 0.6], 'EdgeColor', 'none');
% 
% xlabel('Time since periapsis [min]');
% ylabel('L-shell value');
% title('L-shell Crossings During 1 Orbit');
% legend('L-shell value', 'L = 1.2 threshold', 'In radiation zone');
% grid on;
% set(gca, 'FontSize', 14);

% === Output ===
fprintf('Time in high-radiation zone (L < %.1f): %.2f minutes per orbit.\n', ...
        L_thresh, time_in_radiation);

% Call function for TID calcs
% Time spent in L < 1.2 per orbit [minutes]
dose_rate = 100;      % rad/hour for Si unshielded

compute_TID(n_orbits_in_sc, time_in_radiation, dose_rate);

compute_TID_material('Laser Altimeter', n_orbits_in_sc, time_in_radiation, dose_rate, 'Si');

%dose_rate_GaAs = dose_rate_Si / GaAs_factor;

compute_TID_material('Microwave Radiometer', n_orbits_in_sc, time_in_radiation, dose_rate, 'GaAs');

%% == Rings of Neptune == 

% === Load Orbit and Ring Data ===
load('neptune_ring_data_for_analysis.mat'); 

tolerance_km = 100;  % Distance tolerance for intersection

% === Generate Orbit Positions ===
T_orbit = 2*pi*sqrt(a_final^3 / mu_N);
theta_array = linspace(0, 360, 1500); % degrees
r_orbit = zeros(3, length(theta_array));
r_mag_array = zeros(1, length(theta_array));

for k = 1:length(theta_array)
    theta = theta_array(k);
    [r_vec, ~] = rv_from_oe(a_final, e_final, i_final, RAAN_final, omega_final, theta, mu_N);
    r_orbit(:,k) = r_vec;
    r_mag_array(k) = norm(r_vec);
end

% === Plotting Setup ===
figure; hold on; axis equal; grid on;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Neptune Rings and Science Orbit');

% Plot Neptune
[Xs, Ys, Zs] = sphere(100);
surf(R_N*Xs, R_N*Ys, R_N*Zs, 'FaceColor', [0.3 0.5 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot Science Orbit
plot3(r_orbit(1,:), r_orbit(2,:), r_orbit(3,:), 'k', 'LineWidth', 2, 'DisplayName', 'Science Orbit');

% === Plot and Check Each Ring ===
colors = lines(length(ring_radii));
ring_intersections = {};

for i = 1:length(ring_radii)
    r_ring = double(ring_radii(i));
    width = double(ring_widths(i));
    incl_rad = deg2rad(ring_inclinations(i));

    % Ring in local XY plane
    theta_ring = linspace(0, 2*pi, 300);
    x_ring = r_ring * cos(theta_ring);
    y_ring = r_ring * sin(theta_ring);
    z_ring = zeros(size(x_ring));

    % Rotate ring by inclination (about X axis)
    R_tilt = [1 0 0;
              0 cos(incl_rad) -sin(incl_rad);
              0 sin(incl_rad)  cos(incl_rad)];
    ring_coords = R_tilt * [x_ring; y_ring; z_ring];

    % Plot ring
    plot3(ring_coords(1,:), ring_coords(2,:), ring_coords(3,:), '-', ...
          'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', ring_names{i});

    % Check for intersection
    r_min = r_ring - width/2 - tolerance_km;
    r_max = r_ring + width/2 + tolerance_km;
    in_ring = (r_mag_array >= r_min) & (r_mag_array <= r_max);

    if any(in_ring)
        plot3(r_orbit(1,in_ring), r_orbit(2,in_ring), r_orbit(3,in_ring), ...
              '.', 'Color', colors(i,:), 'MarkerSize', 12);

        ring_intersections{end+1} = struct( ...
            'Ring', ring_names{i}, ...
            'TrueAnomalyRange_deg', theta_array(in_ring([1 end])), ...
            'NumPoints', sum(in_ring));
    else
        ring_intersections{end+1} = struct( ...
            'Ring', ring_names{i}, ...
            'TrueAnomalyRange_deg', [], ...
            'NumPoints', 0);
    end
end

legend show;

% === Print Results ===
fprintf('\n--- Orbit–Ring Intersection Analysis ---\n');
for i = 1:length(ring_intersections)
    data = ring_intersections{i};
    if ~isempty(data.TrueAnomalyRange_deg)
        if numel(data.TrueAnomalyRange_deg) >= 2
            fprintf('✔ Intersects %s between θ = %.1f° and %.1f° (%d points)\n', ...
                data.Ring, data.TrueAnomalyRange_deg(1), data.TrueAnomalyRange_deg(2), data.NumPoints);
        else
            fprintf('✔ Intersects %s at θ ≈ %.1f° (%d points)\n', ...
                data.Ring, data.TrueAnomalyRange_deg(1), data.NumPoints);
        end
    else
        fprintf('✘ No intersection with %s\n', data.Ring);
    end
end


% === 3D Plot Setup ===
figure; hold on; axis equal; grid on;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Neptune Rings and Science Orbit (3D)', 'FontSize', 14);

% Plot Neptune Sphere
[Xs, Ys, Zs] = sphere(100);
surf(R_N*Xs, R_N*Ys, R_N*Zs, 'FaceColor', [0.3 0.5 1], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot Science Orbit
plot3(r_orbit(1,:), r_orbit(2,:), r_orbit(3,:), 'k', 'LineWidth', 2, ...
    'DisplayName', 'Science Orbit');

% === Plot Each Ring in 3D ===
theta_ring = linspace(0, 2*pi, 300);
colors = lines(length(ring_radii));

for i = 1:length(ring_radii)
    r_ring = double(ring_radii(i));
    incl_rad = deg2rad(ring_inclinations(i));
    
    % Ring in XY plane
    x_ring = r_ring * cos(theta_ring);
    y_ring = r_ring * sin(theta_ring);
    z_ring = zeros(size(x_ring));

    % Tilt by ring inclination (about X axis)
    R_tilt = [1 0 0;
              0 cos(incl_rad) -sin(incl_rad);
              0 sin(incl_rad)  cos(incl_rad)];
    ring_coords = R_tilt * [x_ring; y_ring; z_ring];

    % Plot ring
    plot3(ring_coords(1,:), ring_coords(2,:), ring_coords(3,:), '--', ...
          'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', ring_names{i});
end

legend show;
view(45, 25);  % Adjust viewing angle

%% == Plotting Magnetic field lines ==

num_lines = 8;  % Number of meridional lines
r_eq = 3 * R_N; % Maximum radial distance of field lines

% Dipole tilt (Neptune's magnetic axis)
tilt_deg = 47;
tilt_rad = deg2rad(tilt_deg);
R_tilt = [1 0 0;
          0 cos(tilt_rad) -sin(tilt_rad);
          0 sin(tilt_rad)  cos(tilt_rad)];

theta_vals = linspace(0, pi, 100);

% Plot field lines
figure; hold on; axis equal; grid on;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Neptune Magnetic Field Lines (Dipole Model)');

% Neptune sphere
[Xs, Ys, Zs] = sphere(60);
surf(R_N*Xs, R_N*Ys, R_N*Zs, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
    'FaceColor', [0.3 0.5 1]);

phi_vals = linspace(0, 2*pi, num_lines+1);
phi_vals(end) = [];  % Avoid duplicate

for phi = phi_vals
    for L = [1.2, 1.5, 2, 2.5, 3]  % Field lines for different L-shells
        r = L * R_N * (sin(theta_vals).^2);
        x = r .* sin(theta_vals) * cos(phi);
        y = r .* sin(theta_vals) * sin(phi);
        z = r .* cos(theta_vals);

        coords = R_tilt * [x; y; z];  % Apply magnetic tilt
        plot3(coords(1,:), coords(2,:), coords(3,:), 'r--', 'LineWidth', 1.2);
    end
end

% === Plot Science Orbit on Top ===
plot3(r_orbit(1,:), r_orbit(2,:), r_orbit(3,:), 'k', 'LineWidth', 2, ...
       'DisplayName', 'Science Orbit');

legend('Neptune', 'Magnetic Field Lines', 'Science Orbit');
view(40, 20);



%% === FUNCTION ===
function plot_ground_track_neptune_modes(r_ECI_all, jd_array, texture_file, n_orbits)
    sidereal_day_neptune = 16.11 * 3600; % seconds
    omega_N = 2*pi / sidereal_day_neptune; % rad/s

    N = size(r_ECI_all, 2);
    lat = zeros(1, N);
    lon = zeros(1, N);
    alt = zeros(1, N);
    color_orbit = zeros(1, N);
    local_time_deg = zeros(1, N);

    for k = 1:N
        t_sec = (jd_array(k) - jd_array(1)) * 86400;
        theta_rot = mod(omega_N * t_sec, 2*pi);
        R = [ cos(theta_rot),  sin(theta_rot), 0;
             -sin(theta_rot),  cos(theta_rot), 0;
              0,                0,             1 ];
        r_ECEF = R' * r_ECI_all(:,k);

        x = r_ECEF(1); y = r_ECEF(2); z = r_ECEF(3);
        lon(k) = atan2d(y, x);
        r_xy = sqrt(x^2 + y^2);
        lat(k) = atan2d(z, r_xy);
        alt(k) = norm(r_ECEF) - 24764;  % altitude above Neptune radius [km]

        color_orbit(k) = floor((t_sec)/ (2*pi*sqrt((30000)^3/6835100))) + 1;
        local_time_deg(k) = mod(rad2deg(theta_rot), 360); % rotation = local time
    end

    % Load texture
    img = imread(texture_file);
    img = flipud(img);
    Rimg = georefcells([-90 90], [-180 180], size(img));

    % === FIGURE 1: Orbit color-coded ===
    figure;
    ax = axesm('eqdcylin','MapLatLimit',[-90 90],'MapLonLimit',[-180 180],...
               'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on');
    title('Ground Track Coloured by Orbit Number');
    xlabel('Longitude [°]');
    ylabel('Latitude [°]');
    setm(ax,'MLabelLocation',30,'PLabelLocation',30,'FontSize',10);
    geoshow(img, Rimg, 'DisplayType', 'texturemap'); hold on;
    scatterm(lat, lon, 10, color_orbit, 'filled');
    % colormap(lines(n_orbits)); colorbar('Ticks', 1:n_orbits, 'TickLabels', compose('Orbit %d', 1:n_orbits));
    colormap(hsv(n_orbits)); colorbar('Ticks', 1:n_orbits, 'TickLabels', compose('Orbit %d', 1:n_orbits));
    caxis([1 n_orbits]);

    % Start/End Markers
    plotm(lat(1), lon(1), 'y^', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
    plotm(lat(end), lon(end), 'yv', 'MarkerSize', 8, 'MarkerFaceColor', 'y');

    % === FIGURE 2: Altitude colored ===
    figure;
    ax2 = axesm('eqdcylin','MapLatLimit',[-90 90],'MapLonLimit',[-180 180],...
                'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on');
    title('Ground Track Coloured by Altitude');
    xlabel('Longitude [°]');
    ylabel('Latitude [°]');
    setm(ax2,'MLabelLocation',30,'PLabelLocation',30,'FontSize',10);
    geoshow(img, Rimg, 'DisplayType', 'texturemap'); hold on;
    scatterm(lat, lon, 10, alt, 'filled');
colormap(turbo);  % Use vivid colormap to contrast Neptune's blue background
hcb = colorbar; caxis([min(alt) max(alt)]);
ticks = linspace(min(alt), max(alt), 5);
hcb.Ticks = ticks;
hcb.Ticks = ticks(3);
hcb.TickLabels = compose('%.2f km', ticks(3));
hcb.FontSize = 20;
hcb.Position(3) = 0.015;

    % Mark closest approach with triangle
    [~, idx_min] = min(alt);
    plotm(lat(idx_min), lon(idx_min), 'yx', 'MarkerSize', 25, 'MarkerFaceColor', 'y');

    % === FIGURE 3: Time-of-day (Local Longitude Angle)
    figure;
    ax3 = axesm('eqdcylin','MapLatLimit',[-90 90],'MapLonLimit',[-180 180],...
                'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on');
    title('Ground Track Coloured by Local Time of Day');
    xlabel('Longitude [°]');
    ylabel('Latitude [°]');
    setm(ax3,'MLabelLocation',30,'PLabelLocation',30,'FontSize',10);
    geoshow(img, Rimg, 'DisplayType', 'texturemap'); hold on;
    scatterm(lat, lon, 10, local_time_deg, 'filled');
    colormap('hot'); colorbar; caxis([0 360]);

    % Annotate timestamps every ~1/4 orbit
    step = round(N / 12);
    for k = 1:step:N
        time_str = datestr((jd_array(k)-2451545)*86400/86400, 'HH:MM');
        textm(lat(k), lon(k), time_str, 'FontSize', 10, 'Color', 'w');
    end
end


function [r_ECI_out, v_ECI_out, jd_array_after] = rvj_from_pc(r0, v0, i_target_deg, jd0, n_orbits)
    % Find current orbital elements , so it needs to input final r and v at
    % apoapsis at science orbit initial and its epoch there
    mu_N = 6835100;
    [a0, e0, h0, i0, Omega0, omega0, theta0] = oe_from_rv(r0, v0, mu_N);

    % Compute apoapsis radius
    [r_ap, v_ap] = rv_from_oe(a0, e0, h0, i0, Omega0, omega0, 180);
    % Time from current position to apoapsis
    E0 = 2 * atan( sqrt((1-e0)/(1+e0)) * tand(theta0/2) );
    M0 = E0 - e0 * sin(E0);
    if M0 < 0, M0 = M0 + 2*pi; end
    M_ap = pi;  % Mean anomaly at apoapsis
    dM = M_ap - M0;
    if dM < 0, dM = dM + 2*pi; end
    T_orbit = 2 * pi * sqrt(a0^3 / mu_N);
    t_to_ap = dM * T_orbit / (2*pi);
    jd_array_updated = jd0 + t_to_ap / 86400;

    % Delta inclination (radians)
    i0_rad = deg2rad(i0);
    i_target_rad = deg2rad(i_target_deg);
    delta_i = abs(i_target_rad - i0_rad);
    deltaV = sqrt(2) * norm(v_ap) * sqrt(1-cos(delta_i)); %check plane change eq 

    fprintf('Apoapsis plane change from %.2f° to %.2f° requires ΔV = %.3f km/s\n', ...
        i0, i_target_deg, deltaV);

    % Propagate orbit with updated inclination
    [r_after, v_after] = rv_from_oe(a0, e0, i_target_deg, Omega0, omega0, 180, mu_N) % initial r,v, after planechange 
    T_orbit = 2 * pi * sqrt(a0^3 / mu_N);
    t = linspace(0, n_orbits * T_orbit, 1000);
    jd_array_after = jd_array_updated + t / 86400; 

    r_ECI_out = zeros(3, length(t));
    v_ECI_out = zeros(3,length(t));
    for k = 1:length(t)
        M = pi + 2 * pi * t(k) / T_orbit;
        E = kepler_E(M, e0);
        theta = 2 * atan2(sqrt(1 + e0) * sin(E/2), sqrt(1 - e0) * cos(E/2));
        [r_vec, v_vec] = rv_from_oe(a0, e0, i_target_deg, Omega0, omega0, rad2deg(theta), mu_N);
        r_ECI_out(:, k) = r_vec; % does this include the intial one being r_after, v_after
        v_ECI_out(:,k) = v_vec;
    end
end

function E = kepler_E(M, e)
    E = M;  % Initial guess
    for i = 1:100
        E = E - (E - e*sin(E) - M) / (1 - e*cos(E));
    end
end

function compute_TID(orbits, time_per_orbit_min, dose_rate_unshielded)
    % === Shielding thicknesses (mm) and attenuation factors (approximate) ===
    shielding_mm = [0, 1, 3, 5, 10];     % mm Al
    attenuation_factors = [1.00, 0.50, 0.20, 0.10, 0.03];  % dose reduction factors

    % === Total time in belt across mission [hours] ===
    total_time_hr = (time_per_orbit_min / 60) * orbits;

    % === Header ===
    fprintf('--- TID Estimate ---\n');
    fprintf('Total orbits: %d\n', orbits);
    fprintf('Time in belt per orbit: %.1f minutes\n', time_per_orbit_min);
    fprintf('Unshielded dose rate: %.0f rad/hour\n', dose_rate_unshielded);
    fprintf('Total time in radiation zone: %.1f hours\n\n', total_time_hr);

    % === TID Calculation Loop ===
    for k = 1:length(shielding_mm)
        shield = shielding_mm(k);
        attenuation = attenuation_factors(k);
        TID = dose_rate_unshielded * total_time_hr * attenuation;
        fprintf('Shielding: %2d mm Al → TID = %.0f rad(Si)\n', shield, TID);
    end
end

% === Orbit Function ===
function r_all = get_orbit(a, e, i, RAAN, omega, mu, t_array)
    T_orbit = 2 * pi * sqrt(a^3 / mu);
    r_all = zeros(3, length(t_array));
    for k = 1:length(t_array)
        M = 2 * pi * t_array(k) / T_orbit;

        % Solve Kepler's Equation
        E = M;
        for iter = 1:100
            E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
        end
        theta = 2 * atan2(sqrt(1 + e) * sin(E/2), sqrt(1 - e) * cos(E/2));
        r = a * (1 - e^2) / (1 + e * cos(theta));
        x_p = r * cos(theta);
        y_p = r * sin(theta);
        z_p = 0;

        % Rotate from perifocal to ECI
        R = [cos(RAAN)*cos(omega)-sin(RAAN)*sin(omega)*cos(i), -cos(RAAN)*sin(omega)-sin(RAAN)*cos(omega)*cos(i), sin(RAAN)*sin(i);
             sin(RAAN)*cos(omega)+cos(RAAN)*sin(omega)*cos(i), -sin(RAAN)*sin(omega)+cos(RAAN)*cos(omega)*cos(i), -cos(RAAN)*sin(i);
             sin(omega)*sin(i),                                cos(omega)*sin(i),                               cos(i)];

        r_vec = R * [x_p; y_p; z_p];
        r_all(:,k) = r_vec;
    end
end


function compute_TID_material(instrument_name, orbits, time_per_orbit_min, dose_rate_unshielded, material)
    % === Shielding thicknesses (mm) and attenuation factors (approximate) ===
    shielding_mm = [0, 1, 3, 5, 10];     
    attenuation_factors = [1.00, 0.50, 0.20, 0.10, 0.03];  

    % === Material scaling (GaAs absorbs less ionising dose than Si) ===
    switch lower(material)
        case 'si'
            material_factor = 1.00;
        case 'gaas'
            material_factor = 0.30;  % Approximate reduction
        otherwise
            error('Unsupported material: use "Si" or "GaAs"');
    end

    % === Total time in belt across mission [hours] ===
    total_time_hr = (time_per_orbit_min / 60) * orbits;

    % === Header ===
    fprintf('\n--- TID Estimate for %s (%s) ---\n', instrument_name, material);
    fprintf('Total orbits: %d\n', orbits);
    fprintf('Time in belt per orbit: %.1f minutes\n', time_per_orbit_min);
    fprintf('Unshielded dose rate: %.0f rad/hour\n', dose_rate_unshielded);
    fprintf('Total time in radiation zone: %.1f hours\n\n', total_time_hr);

    % === TID Calculation Loop ===
    for k = 1:length(shielding_mm)
        shield = shielding_mm(k);
        attenuation = attenuation_factors(k);
        TID_si = dose_rate_unshielded * total_time_hr * attenuation;
        TID_material = TID_si * material_factor;
        fprintf('Shielding: %2d mm Al → TID = %.0f rad(%s)\n', shield, TID_material, upper(material));
    end
end
