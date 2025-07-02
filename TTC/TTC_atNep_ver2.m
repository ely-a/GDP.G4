clear;
clc;

% Constants and configuration
c = 299792.458; % Speed of light (km/s)
f0 = 8.4e9; % X-band frequency (Hz)
critical_angle = 3; % Solar conjunction critical angle (degrees)
min_elevation = 10; % Minimum elevation angle for visibility (degrees)
R_neptune = 24622; % Neptune radius (km)
mu_neptune = 6.836529e6; % Neptune gravitational parameter (km³/s²)
arrival_date = datetime('2040-06-01');
mission_end = arrival_date + calyears(5);

% DSN stations [Latitude, Longitude]
stations = [35.42, -116.89;   % Goldstone
            40.43,  -4.25;    % Madrid
           -35.40, 148.98];   % Canberra

% Coarse time grid (1-day intervals for ephemeris)
time_coarse = arrival_date:days(1):mission_end;
jd_coarse = juliandate(time_coarse);

% Planet data (so you don't need the extension)
load('E_and_N_planet_data.mat');

% Fine grids for mission-wide analysis (daily)
time_daily = arrival_date:days(1):mission_end;
jd_daily = juliandate(time_daily);

% Interpolate positions and velocities to daily grid
earth_pos_daily = interp1(jd_coarse, earth_pos_coarse, jd_daily, 'spline');
earth_vel_daily = interp1(jd_coarse, earth_vel_coarse, jd_daily, 'spline');
neptune_pos_daily = interp1(jd_coarse, neptune_pos_coarse, jd_daily, 'spline');
neptune_vel_daily = interp1(jd_coarse, neptune_vel_coarse, jd_daily, 'spline');

% Vectorized calculations for daily analyses
r_vec_daily = neptune_pos_daily - earth_pos_daily;
distances_daily = vecnorm(r_vec_daily, 2, 2);
sun_earth_daily = vecnorm(earth_pos_daily, 2, 2);
sun_neptune_daily = vecnorm(neptune_pos_daily, 2, 2);

% SEP angle calculation
sep_angle = acosd(dot(earth_pos_daily, neptune_pos_daily, 2) ./ (sun_earth_daily .* sun_neptune_daily));
conjunction_binary = double(sep_angle < critical_angle);

% Light time calculation
light_time = distances_daily / c;

% Doppler shift calculation (vectorized)
radial_velocity = dot(r_vec_daily./distances_daily, neptune_vel_daily - earth_vel_daily, 2);
frequency_shift = f0 * (1 - radial_velocity/c);

% Calculate RA and Dec for all daily positions
dec = asind(r_vec_daily(:,3) ./ distances_daily);
ra = atan2d(r_vec_daily(:,2), r_vec_daily(:,1));
ra(ra < 0) = ra(ra < 0) + 360;

% Compute GMST for all daily times (vectorized)
gmst = siderealTime(jd_daily); % In degrees

% Preallocate visibility matrix
station_visibility = false(length(time_daily), size(stations,1));

% Vectorized visibility calculation for each station
for j = 1:size(stations,1)
    lat = stations(j, 1);
    lon = stations(j, 2);
    
    % Local Sidereal Time
    lst = mod(gmst + lon, 360);
    
    % Hour Angle
    ha = mod(lst - ra + 180, 360) - 180;
    
    % Elevation calculation
    sin_el = sind(lat)*sind(dec) + cosd(lat)*cosd(dec).*cosd(ha);
    station_visibility(:,j) = asind(sin_el) > min_elevation;
end

% Overall visibility: visible if any station can see it
visibility_binary = ~any(station_visibility, 2);

% Time since arrival (days)
t_daily = days(time_daily - arrival_date);



%% Detailed DSN Visibility Analysis
% Load orbital elements
load('Aerobrakes_OE.mat', 'a_brakes', 'e_brakes', 'i_brakes', 'omega_brakes', 'RAAN_brakes');
load('ScienceOrbit.mat', 'a_final', 'e_final', 'i_final', 'omega_final', 'RAAN_final', 't_days');

% User input for target day
day_input = input('\nEnter day since arrival for detailed DSN visibility analysis (0-1825): ');
if day_input < 0 || day_input > 1825
    error('Day must be between 0 and 1825 (5 years)');
end

% Create fine grid (200 points for the selected day)
n_grid = 200;
target_date = arrival_date + days(day_input);
time_fine = target_date + seconds(linspace(0, 86400, n_grid));
jd_fine = juliandate(time_fine);

% Interpolate ephemerides to fine grid
earth_pos_fine = interp1(jd_coarse, earth_pos_coarse, jd_fine, 'spline');
neptune_pos_fine = interp1(jd_coarse, neptune_pos_coarse, jd_fine, 'spline');

% Preallocate arrays
visibility = false(n_grid, 3);  % Per station visibility
occultation = false(n_grid, 1); % Occultation flags
conjunction = false(n_grid, 1); % Solar conjunction flags

% Process each time point
for i = 1:n_grid
    % Get orbital elements at current time
    t_current = days(time_fine(i) - arrival_date);
    if t_current <= t_days
        % Linear interpolation during aerobraking
        frac = t_current / t_days;
        idx = min(floor(frac * (length(a_brakes)-1) + 1), length(a_brakes)-1);
        weight = frac*(length(a_brakes)-1) - (idx-1);
        
        a = a_brakes(idx) + weight*(a_brakes(idx+1)-a_brakes(idx));
        e = e_brakes(idx) + weight*(e_brakes(idx+1)-e_brakes(idx));
        inc = i_brakes(idx) + weight*(i_brakes(idx+1)-i_brakes(idx));
        omega = omega_brakes(idx) + weight*(omega_brakes(idx+1)-omega_brakes(idx));
        RAAN = RAAN_brakes(idx) + weight*(RAAN_brakes(idx+1)-RAAN_brakes(idx));
    else
        % Final science orbit
        a = a_final;
        e = e_final;
        inc = i_final;
        omega = omega_final;
        RAAN = RAAN_final;
    end
    
    % Compute spacecraft position relative to Neptune
    [r_sc_neptune, ~] = kepler_orbit(a, e, deg2rad(inc), deg2rad(omega), deg2rad(RAAN), mu_neptune, time_fine(i), arrival_date);
    r_sc_neptune = r_sc_neptune';  % Ensure row vector
    
    % Position in solar system
    r_sc_sun = neptune_pos_fine(i,:) + r_sc_neptune;
    r_earth_to_sc = r_sc_sun - earth_pos_fine(i,:);
    distance_earth_sc = norm(r_earth_to_sc);
    
    % Compute exact SEP angle for spacecraft
    earth_to_sun = earth_pos_fine(i,:);
    earth_to_probe = r_earth_to_sc;
    sep_angle_fine = acosd(dot(earth_to_sun, earth_to_probe) / (norm(earth_to_sun) * norm(earth_to_probe)));
    conjunction(i) = sep_angle_fine < critical_angle;
    
    % Compute RA and Dec
    dec_val = asind(r_earth_to_sc(3) / distance_earth_sc);
    ra_val = atan2d(r_earth_to_sc(2), r_earth_to_sc(1));
    if ra_val < 0
        ra_val = ra_val + 360;
    end
    
    % Occultation check
    r_earth_to_neptune = neptune_pos_fine(i,:) - earth_pos_fine(i,:);
    dist_earth_neptune = norm(r_earth_to_neptune);
    
    % Angle between vectors to Neptune and spacecraft from Earth
    cos_theta = dot(r_earth_to_neptune, r_earth_to_sc) / (dist_earth_neptune * distance_earth_sc);
    angular_separation = acosd(cos_theta);
    
    % Angular radius of Neptune
    angular_radius_neptune = asind(R_neptune / dist_earth_neptune);
    
    % Check if spacecraft is behind Neptune and within angular radius
    if distance_earth_sc > dist_earth_neptune && angular_separation < angular_radius_neptune
        occultation(i) = true;
    else
        occultation(i) = false;
    end
    
    % Check visibility for each station (only if not occulted and not in conjunction)
    gmst_val = siderealTime(jd_fine(i));
    for j = 1:3
        if ~occultation(i) && ~conjunction(i)
            lat = stations(j,1);
            lon = stations(j,2);
            lst = mod(gmst_val + lon, 360);
            ha = mod(lst - ra_val + 180, 360) - 180;
            sin_el = sind(lat)*sind(dec_val) + cosd(lat)*cosd(dec_val)*cosd(ha);
            elevation = asind(sin_el);
            visibility(i,j) = elevation > min_elevation;
        else
            visibility(i,j) = false;
        end
    end
end

%% PLOTS

% Plot mission-wide results
figure('Name', 'Doppler Analysis', 'Position', [100 100 1200 800]);
subplot(2,1,1);
plot(t_daily, radial_velocity);
xlabel('Time Since Arrival (days)');
ylabel('Radial Velocity (km/s)');
title('Radial Velocity vs Time');
grid on;

subplot(2,1,2);
plot(t_daily, frequency_shift);
xlabel('Time Since Arrival (days)');
ylabel('Frequency (Hz)');
title('Frequency Shift vs Time');
grid on;

figure('Name', 'Solar Conjunction Analysis', 'Position', [100 100 1200 800]);
subplot(2,1,1);
plot(t_daily, sep_angle);
xlabel('Time Since Arrival (days)');
ylabel('SEP Angle (degrees)');
title('Sun-Earth-Probe Angle');
grid on;

subplot(2,1,2);
stairs(t_daily, conjunction_binary, 'LineWidth', 1.5);
ylim([-0.1 1.1]);
xlabel('Time Since Arrival (days)');
ylabel('Conjunction Status (1 = in conjunction)');
title('Solar Conjunction Binary');
grid on;

figure('Name', 'DSN Visibility (Mission-wide)', 'Position', [100 100 1200 400]);
stairs(t_daily, visibility_binary, 'LineWidth', 1.5);
ylim([-0.1 1.1]);
xlabel('Time Since Arrival (days)');
ylabel('Visibility Status (1 = not visible)');
title('DSN Visibility Binary (Without Occultation)');
grid on;

figure('Name', 'Light Time Delay', 'Position', [100 100 1200 400]);
plot(t_daily, light_time);
xlabel('Time Since Arrival (days)');
ylabel('Signal Travel Time (seconds)');
title('Light Time Delay from Neptune to Earth');
grid on;

% Generate detailed visibility plot
plot_detailed_visibility(time_fine, visibility, occultation, conjunction, stations, day_input);

%% Helper functions
function [r, v] = kepler_orbit(a, e, i, omega, RAAN, mu, current_time, arrival_date)
    % Keplerian orbit propagator
    t_since_arrival = seconds(current_time - arrival_date);
    
    % Calculate orbital period
    T = 2*pi*sqrt(a^3/mu);
    
    % Mean motion
    n = 2*pi/T;
    
    % Mean anomaly
    M = n * t_since_arrival;
    M = mod(M, 2*pi);
    
    % Solve Kepler's equation
    E = M; % Initial guess
    for iter = 1:50
        E_new = M + e*sin(E);
        if abs(E_new - E) < 1e-12
            break;
        end
        E = E_new;
    end
    
    % True anomaly
    sin_theta = sqrt(1-e^2)*sin(E)/(1-e*cos(E));
    cos_theta = (cos(E)-e)/(1-e*cos(E));
    theta = atan2(sin_theta, cos_theta);
    
    % Perifocal coordinates
    r_peri = a*(1 - e^2)/(1 + e*cos(theta)) * [cos(theta); sin(theta); 0];
    v_peri = sqrt(mu/(a*(1-e^2))) * [-sin(theta); e+cos(theta); 0];
    
    % Rotation matrices
    R_omega = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
    R_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R_RAAN = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
    
    % Transform to inertial frame
    R = R_RAAN * R_i * R_omega;
    r = R * r_peri;
    v = R * v_peri;
end

function plot_detailed_visibility(time, visibility, occultation, conjunction, stations, day_input)
    % Convert to hours since start of day
    hoursT = hours(time - time(1));
    
    % Calculate communication availability
    station_available = visibility & ~occultation & ~conjunction;
    any_station_visible = any(visibility, 2);
    combined_visibility = any(visibility, 2) & ~occultation;
    overall_comm = combined_visibility & ~conjunction;
    
    figure('Position', [100, 100, 1200, 1000], 'Name', sprintf('Detailed DSN Visibility - Day %d', day_input));
    
    % 1. Individual station visibility
    subplot(5,1,1);
    hold on;
    colors = lines(3);
    for j = 1:3
        plot(hoursT, visibility(:,j), 'Color', colors(j,:), 'LineWidth', 1.5);
    end
    title('Station Visibility (Without Constraints)');
    xlabel('Hours Since 00:00 UTC');
    ylabel('Visible');
    legend('Goldstone', 'Madrid', 'Canberra', 'Location', 'best');
    grid on;
    ylim([-0.1 1.1]);
    set(gca, 'YTick', [0 1], 'YTickLabel', {'No', 'Yes'});
    
    % 2. Combined DSN visibility (without conjunction)
    subplot(5,1,2);
    stairs(hoursT, combined_visibility, 'LineWidth', 2, 'Color', [0 0.5 0]);
    title('Combined Visibility (Not Occulted)');
    xlabel('Hours Since 00:00 UTC');
    ylabel('Visible');
    grid on;
    ylim([-0.1 1.1]);
    set(gca, 'YTick', [0 1], 'YTickLabel', {'No', 'Yes'});
    
    % 3. Occultation status
    subplot(5,1,3);
    stairs(hoursT, occultation, 'LineWidth', 2, 'Color', [0.8 0 0]);
    title('Occultation by Neptune');
    xlabel('Hours Since 00:00 UTC');
    ylabel('Occulted');
    grid on;
    ylim([-0.1 1.1]);
    set(gca, 'YTick', [0 1], 'YTickLabel', {'No', 'Yes'});
    
    % 4. Solar conjunction status
    subplot(5,1,4);
    stairs(hoursT, conjunction, 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
    title('Solar Conjunction');
    xlabel('Hours Since 00:00 UTC');
    ylabel('In Conjunction');
    grid on;
    ylim([-0.1 1.1]);
    set(gca, 'YTick', [0 1], 'YTickLabel', {'No', 'Yes'});
    
    % 5. Overall communication availability
    subplot(5,1,5);
    stairs(hoursT, overall_comm, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    title('Overall Communication Availability');
    xlabel('Hours Since 00:00 UTC');
    ylabel('Available');
    grid on;
    ylim([-0.1 1.1]);
    set(gca, 'YTick', [0 1], 'YTickLabel', {'No', 'Yes'});
    
    % Add solar conjunction status for the day
    if any(conjunction)
        conjunction_status = 'IN SOLAR CONJUNCTION at some point today';
        color = [0.8 0 0];
    else
        conjunction_status = 'NOT in solar conjunction today';
        color = [0 0.5 0];
    end
    
    annotation('textbox', [0.3, 0.02, 0.4, 0.04], 'String', ...
        sprintf('Day %d Status: %s', day_input, conjunction_status), ...
        'Color', color, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    
end

