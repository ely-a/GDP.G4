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

% Planet data
load('E_and_N_planet_data.mat');

% Load orbital elements
load('Aerobrakes_OE.mat', 'a_brakes', 'e_brakes', 'i_brakes', 'omega_brakes', 'RAAN_brakes');
load('ScienceOrbit.mat', 'a_final', 'e_final', 'i_final', 'omega_final', 'RAAN_final', 't_days');

% Animation parameters
days_animate = ceil(t_days):2:1825; % Every other day from science mission
num_frames = length(days_animate);
frame_rate = 60;
animation_duration = 40; % seconds
frames_per_plot = round(frame_rate * animation_duration / num_frames);

% Set up figure for animation
f = figure('Position', [100, 100, 1200, 200]);
ax = axes('Parent', f);
hold(ax, 'on');

% Preallocate video structure
final_animation(num_frames) = struct('cdata',[],'colormap',[]);

% Generate each frame
for idx = 1:num_frames
    day_input = days_animate(idx);
    target_date = arrival_date + days(day_input);
    
    % Create fine grid (200 points)
    n_grid = 200;
    time_fine = target_date + seconds(linspace(0, 86400, n_grid));
    jd_fine = juliandate(time_fine);
    
    % Interpolate ephemerides
    earth_pos_fine = interp1(jd_coarse, earth_pos_coarse, jd_fine, 'spline');
    neptune_pos_fine = interp1(jd_coarse, neptune_pos_coarse, jd_fine, 'spline');
    
    % Preallocate arrays
    visibility = false(n_grid, 3);
    occultation = false(n_grid, 1);
    conjunction = false(n_grid, 1);
    
    % Calculate positions and visibility
    for i = 1:n_grid
        % Determine orbital elements
        t_current = days(time_fine(i) - arrival_date);
        if t_current <= t_days
            frac = t_current / t_days;
            idx_oe = min(floor(frac * (length(a_brakes)-1) + 1), length(a_brakes)-1);
            weight = frac*(length(a_brakes)-1) - (idx_oe-1);
            
            a = a_brakes(idx_oe) + weight*(a_brakes(idx_oe+1)-a_brakes(idx_oe));
            e = e_brakes(idx_oe) + weight*(e_brakes(idx_oe+1)-e_brakes(idx_oe));
            inc = i_brakes(idx_oe) + weight*(i_brakes(idx_oe+1)-i_brakes(idx_oe));
            omega = omega_brakes(idx_oe) + weight*(omega_brakes(idx_oe+1)-omega_brakes(idx_oe));
            RAAN = RAAN_brakes(idx_oe) + weight*(RAAN_brakes(idx_oe+1)-RAAN_brakes(idx_oe));
        else
            a = a_final;
            e = e_final;
            inc = i_final;
            omega = omega_final;
            RAAN = RAAN_final;
        end
        
        % Compute spacecraft position
        [r_sc_neptune, ~] = kepler_orbit(a, e, deg2rad(inc), deg2rad(omega),...
                                         deg2rad(RAAN), mu_neptune, time_fine(i), arrival_date);
        r_sc_neptune = r_sc_neptune';
        r_sc_sun = neptune_pos_fine(i,:) + r_sc_neptune;
        r_earth_to_sc = r_sc_sun - earth_pos_fine(i,:);
        distance_earth_sc = norm(r_earth_to_sc);
        
        % Check solar conjunction
        earth_to_sun = earth_pos_fine(i,:);
        earth_to_probe = r_earth_to_sc;
        sep_angle_fine = acosd(dot(earth_to_sun, earth_to_probe) /...
                          (norm(earth_to_sun) * norm(earth_to_probe)));
        conjunction(i) = sep_angle_fine < critical_angle;
        
        % Check occultation
        r_earth_to_neptune = neptune_pos_fine(i,:) - earth_pos_fine(i,:);
        dist_earth_neptune = norm(r_earth_to_neptune);
        cos_theta = dot(r_earth_to_neptune, r_earth_to_sc) /...
                   (dist_earth_neptune * distance_earth_sc);
        angular_separation = acosd(cos_theta);
        angular_radius_neptune = asind(R_neptune / dist_earth_neptune);
        
        occultation(i) = (distance_earth_sc > dist_earth_neptune) &&...
                         (angular_separation < angular_radius_neptune);
        
        % Check station visibility
        if ~occultation(i) && ~conjunction(i)
            dec_val = asind(r_earth_to_sc(3) / distance_earth_sc);
            ra_val = atan2d(r_earth_to_sc(2), r_earth_to_sc(1));
            ra_val = mod(ra_val, 360);
            
            gmst_val = siderealTime(jd_fine(i));
            for j = 1:3
                lat = stations(j,1);
                lon = stations(j,2);
                lst = mod(gmst_val + lon, 360);
                ha = mod(lst - ra_val + 180, 360) - 180;
                sin_el = sind(lat)*sind(dec_val) + cosd(lat)*cosd(dec_val)*cosd(ha);
                elevation = asind(sin_el);
                visibility(i,j) = elevation > min_elevation;
            end
        end
    end
    
    % Compute overall communication availability
    overall_comm = any(visibility, 2) & ~occultation & ~conjunction;
    hoursT = hours(time_fine - time_fine(1));
    
    % Plot communication availability
    font = 25;
    cla(ax);
    stairs(ax, hoursT, overall_comm, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
    title(ax, sprintf('Overall Communication Availability - Day %d', day_input), 'FontSize', font);
    xlabel(ax, 'Hours Since 00:00 UTC', "FontSize", font);
    ylabel(ax, 'Available', "FontSize", font);
    grid(ax, 'on');
    ylim(ax, [-0.1 1.1]);
    set(ax, 'YTick', [0 1], 'YTickLabel', {'No', 'Yes'}, 'FontSize', 16);
    
    % Capture frame
    final_animation(idx) = getframe(f);
end

% Play animation
figure;
movie(final_animation, 1, frame_rate / frames_per_plot);

% Save animation as video
v = VideoWriter('CommAvailability.mp4', 'MPEG-4');
v.FrameRate = frame_rate / frames_per_plot;
open(v);
writeVideo(v, final_animation);
close(v);

% Helper functions
function [r, v] = kepler_orbit(a, e, i, omega, RAAN, mu, current_time, arrival_date)
    t_since_arrival = seconds(current_time - arrival_date);
    T = 2*pi*sqrt(a^3/mu);
    n = 2*pi/T;
    M = n * t_since_arrival;
    M = mod(M, 2*pi);
    
    E = M;
    for iter = 1:50
        E_new = M + e*sin(E);
        if abs(E_new - E) < 1e-12
            break;
        end
        E = E_new;
    end
    
    sin_theta = sqrt(1-e^2)*sin(E)/(1-e*cos(E));
    cos_theta = (cos(E)-e)/(1-e*cos(E));
    theta = atan2(sin_theta, cos_theta);
    
    r_peri = a*(1 - e^2)/(1 + e*cos(theta)) * [cos(theta); sin(theta); 0];
    v_peri = sqrt(mu/(a*(1-e^2))) * [-sin(theta); e+cos(theta); 0];
    
    R_omega = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
    R_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R_RAAN = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
    
    R = R_RAAN * R_i * R_omega;
    r = R * r_peri;
    v = R * v_peri;
end

function gmst = siderealTime(jd)
    T = (jd - 2451545.0) / 36525.0;
    gmst = 280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*T.^2 - T.^3/38710000;
    gmst = mod(gmst, 360);
end