% =========================================================================
% SETUP & CONSTANTS
% =========================================================================
clear; clc; close all;
% Load trajectory and planetary data
sc_data = load('sc_data.mat');
planet_data = load('planet_data.mat');
orbit_oe = load('ScienceOrbit.mat');
% NOTE: The orbital phase requires ephemeris data for Earth and Neptune
% covering the 2040-2045 period. In a real mission, this would come from a source like JPL's HORIZONS.
load('E_and_N_planet_data.mat'); % This loads earth_pos_coarse and neptune_pos_coarse
% Constants
c = 299792.458;                 % Speed of light (km/s)
AU_km = 149597870.7;            % 1 AU in km
f_dl = 8.4e9;                   % Downlink frequency (Hz) - X-band
f_ul = 7.9e9;                   % Uplink frequency (Hz) - X-band
k = 1.38065e-23;                % Boltzmann constant (J/K)
mu_neptune = 6.836529e12 / 1e9;  % Neptune gravitational parameter (km^3/s^2)
% =========================================================================
% PHASE 1: INTERPLANETARY CRUISE TRAJECTORY (2032 - 2040)
% =========================================================================
% Time setup for the cruise phase
t_start_cruise = datetime('2032-03-23');
t_end_cruise = datetime('2040-06-01');
time_vec_cruise = t_start_cruise:days(1):t_end_cruise;
% Extract position vectors from loaded data
sc_pos_cruise = sc_data.r_out;
earth_pos_cruise = planet_data.earth_pos;
% =========================================================================
% PHASE 2: NEPTUNE SCIENCE ORBIT (2040 - 2045)
% =========================================================================
% Time setup for the 1825-day science orbit phase
t_start_orbit = t_end_cruise + days(1);
t_end_orbit = t_start_orbit + days(1825);
time_vec_orbit = t_start_orbit:days(1):t_end_orbit; % This has 1826 entries
% Get Earth and Neptune positions during the orbital phase
% Ensure the loaded orbital planet data matches the orbital time vector length
num_orbit_days = length(time_vec_orbit); % This will be 1826
earth_pos_orbit = earth_pos_coarse(1:num_orbit_days, :);
neptune_pos_orbit = neptune_pos_coarse(1:num_orbit_days, :);
% Pre-allocate array for spacecraft's heliocentric position
sc_pos_orbit = zeros(length(time_vec_orbit), 3); % This will have 1826 rows
% Calculate spacecraft position in orbit around Neptune for each day
a = orbit_oe.a_final;
e = orbit_oe.e_final;
inc = orbit_oe.i_final;
raan = orbit_oe.RAAN_final;
aop = orbit_oe.omega_final;
n = sqrt(mu_neptune / a^3); % Mean motion (rad/s)
T_seconds = 2*pi/n;         % Orbital period in seconds
for i = 1:length(time_vec_orbit)
    % Time since start of orbit in seconds
    t = (i-1) * 86400;
    
    % Mean Anomaly
    M = n * mod(t, T_seconds);
    
    % Solve Kepler's Equation for Eccentric Anomaly (E) using Newton's method
    E = M; % Initial guess
    for k_iter = 1:10
        E = E - (E - e*sin(E) - M) / (1 - e*cos(E));
    end
    
    % True Anomaly
    nu = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    
    % Convert orbital elements to state vector relative to Neptune
    [r_sc_wrt_neptune, ~] = coe2rv(a, e, inc, raan, aop, nu, mu_neptune);
    
    % Calculate heliocentric position of the spacecraft
    % r_sc_heliocentric = r_neptune_heliocentric + r_sc_wrt_neptune
    sc_pos_orbit(i, :) = neptune_pos_orbit(i, :) + r_sc_wrt_neptune';
end
% =========================================================================
% COMBINE TRAJECTORIES & CALCULATE LINK MARGIN
% =========================================================================
% Combine time vectors and position data from both phases
time_vec = [time_vec_cruise, time_vec_orbit]';
sc_pos = [sc_pos_cruise; sc_pos_orbit];
earth_pos = [earth_pos_cruise; earth_pos_orbit];
N = length(time_vec); % Total number of data points for the entire mission
% Compute Earth-SC distance (km) for the entire mission
r_earth_sc = sc_pos - earth_pos;
dist_earth_sc_km = vecnorm(r_earth_sc, 2, 2);
dist_AU = dist_earth_sc_km / AU_km;
% =========================================================================
% LINK BUDGET PARAMETERS (Aligned with Detailed Script where applicable)
% --- Calculations for constant parameters are performed ONCE here ---
% =========================================================================
% --- Downlink Specific Inputs (from Detailed Script) ---
% Spacecraft HGA Characteristics
dish_diameter     = 2.5;    % meters (Spacecraft HGA)
dish_efficiency   = 0.75;   % Spacecraft HGA efficiency
% Ground DSN Antenna Characteristics
ground_diameter   = 134.3;     % meters (Ground DSN antenna, as per your initial script's G_r_dl value)
ground_efficiency = 0.75;   % Ground DSN efficiency
% Transmitter Power (Orbiter)
P_tx_orbiter = 100; % Watts (Orbiter TX power to HGA for Earth)
% Link Performance Requirements (for Downlink)
requiredEbN0_down    = 0.8;  % Required Eb/N0 (dB)
implementationLoss   = 1.0;  % dB (Total system degradation due to non-ideal components/synchronization)
% Component Losses (dB) - Orbiter TX Side (Summed to 'total_tx_losses_orbiter')
tx_amp_backoff_loss_dB = 0.5;
tx_feeder_loss_dB      = 0.5;
tx_radome_loss_dB      = 0.2;
tx_other_loss_dB       = 1.0;
total_tx_losses_orbiter = tx_amp_backoff_loss_dB + tx_feeder_loss_dB + tx_radome_loss_dB + tx_other_loss_dB; % Sum of positive losses
% Component Losses (dB) - Ground RX Side (Summed to 'total_rx_losses_ground')
rx_feeder_loss_dB      = 0.5;
rx_radome_loss_dB      = 0.2;
rx_other_loss_dB       = 0.3;
total_rx_losses_ground = rx_feeder_loss_dB + rx_radome_loss_dB + rx_other_loss_dB; % Sum of positive losses
% Antenna Pointing Error (Fixed, assumed system capability)
system_pointing_error_deg = 0.20; % degrees
% Polarization Mismatch Loss
polarization_mismatch_loss_dB = 0.2;
% Atmospheric Effects Inputs (for ITU models)
elevation_angle_deg             = 10;      % degrees
rain_rate_mm_hr                 = 5;       % mm/hr
polarization_tilt_angle_deg     = 45;      % degrees (for rain model)
cloud_temp_celsius              = 0;       % Celsius
cloud_liquid_water_density_gm3  = 0.5;     % g/m^3
air_temp_celsius                = 15;      % Celsius
atmospheric_pressure_hPa        = 1013.25; % hPa
site_latitude_deg               = 35.3;    % degrees
water_vapor_density_gm3         = 7.5;     % g/m^3
% Downlink System Noise Temperature (Kelvin) - Earth Ground Station <--- ADDED DEFINITION
T_s_dl = 20.6;
% --- Derived Downlink Parameters (using helper functions - calculated once) ---
% SC HGA Gain (G_t_dl)
G_t_dl = computeAntennaGain(dish_diameter, dish_efficiency, f_dl);
% Ground DSN Gain (G_r_dl)
G_r_dl = computeAntennaGain(ground_diameter, ground_efficiency, f_dl);
% Spacecraft HGA Beamwidth (needed for pointing loss calculation)
beamwidth_hga_down = computeBeamwidth(dish_diameter, f_dl);
% Antenna Pointing Loss (Downlink)
antenna_pointing_loss_dl = computePointingLoss(system_pointing_error_deg, beamwidth_hga_down);
% Atmospheric Losses (Downlink) - calculated once as inputs are constant
rain_att_dl = computeRainAttenuationITU(f_dl, elevation_angle_deg, rain_rate_mm_hr, polarization_tilt_angle_deg, site_latitude_deg);
cloud_att_dl = computeCloudAttenuationITU(f_dl, cloud_temp_celsius, cloud_liquid_water_density_gm3);
gas_att_dl = computeGasAttenuationITU(f_dl, air_temp_celsius, atmospheric_pressure_hPa, water_vapor_density_gm3, elevation_angle_deg);
total_atmospheric_losses_dl = rain_att_dl + cloud_att_dl + gas_att_dl;
% Effective Isotropic Radiated Power (EIRP) for Downlink - calculated once
EIRP_dl = computeEIRP(P_tx_orbiter, G_t_dl, total_tx_losses_orbiter); % Ptx + Gtx - Ltx_comp
% --- Uplink Parameters (Unchanged from original - No Detailed Counterpart Provided) ---
P_t_ul = 20000;                 % Uplink Transmit Power (Watts) - Earth Ground Station
G_t_ul = 73.5;                  % Uplink Transmitter Antenna Gain (dBi) - Earth Ground Station
L_tx_ul = -2.2;                 % Uplink Transmitter Line Losses (dB) - Earth Ground Station (original lumped loss)
EIRP_ul = 10*log10(P_t_ul*1000) + L_tx_ul + G_t_ul; % Uplink EIRP (dBm)
G_r_ul = 44.47;                 % Uplink Receiver Antenna Gain (dBi) - Spacecraft
T_s_ul = 100;                   % Uplink System Noise Temperature (Kelvin) - Spacecraft
L_rx_ul = -1;                   % Uplink Receiver Line Losses (dB) - Spacecraft (original lumped loss)
L_atm_ul = -1.945;              % Uplink Atmospheric Loss (dB) - Earth's atmosphere (original lumped loss)
% Uplink Noise Power Spectral Density (N0_ul_dBm) - calculated once
N0_ul_dBm = 10*log10(k*T_s_ul*1000); 
req_EbN0_ul = 1;                % Uplink Required Eb/N0 (dB)
% --- Common Data Rate Parameter ---
R_bps = 42500;                  % Data Rate (bits per second) - applies to both links for margin calcs
R_dBHz = 10*log10(R_bps);       % Data Rate (dBHz)
% Initialize arrays for margins
downlink_margin = zeros(N,1);
uplink_margin = zeros(N,1);
% Loop through the entire mission timeline
for i = 1:N
    dist_m = dist_earth_sc_km(i) * 1000; % Distance in meters for FSPL - this changes per day
    % --- Downlink Calculations (Only distance-dependent calculations inside loop) ---
    FSPL_dl = computeFSPL(dist_m, f_dl); % Free Space Path Loss for Downlink
    % Total propagation losses (all losses that reduce signal strength)
    % This sum includes the distance-dependent FSPL and other constant losses
    total_propagation_losses_dl = FSPL_dl + total_atmospheric_losses_dl + polarization_mismatch_loss_dB + antenna_pointing_loss_dl;
    % Calculate Received Power (P_r_dl)
    % Prx = EIRP - PropagationLosses + Grx - Lrx
    P_r_dl = EIRP_dl - total_propagation_losses_dl + G_r_dl - total_rx_losses_ground; % All values are in dB/dBm
    % Calculate Energy per Bit to Noise Power Spectral Density (Eb/N0)
    EbN0_dl = computeEbN0(P_r_dl, R_bps, T_s_dl, k);
    
    % Calculate Downlink Link Margin
    downlink_margin(i) = EbN0_dl - requiredEbN0_down - implementationLoss;
    % --- Uplink Calculations (Only distance-dependent calculations inside loop) ---
    FSPL_ul = computeFSPL(dist_m, f_ul);
    P_r_ul = EIRP_ul - FSPL_ul + L_atm_ul + G_r_ul; % L_atm_ul is the original lumped atmospheric loss
    
    % EbN0_ul calculation using original lumped parameters (consistent with original script's uplink)
    EbN0_ul = P_r_ul - R_dBHz - N0_ul_dBm + L_rx_ul; 
    % Calculate Uplink Link Margin
    uplink_margin(i) = EbN0_ul - req_EbN0_ul;
end
% =========================================================================
% RESULTS TABLE (EXTENDED INTERVALS)
% =========================================================================
% Generate intervals for the full mission duration
interval_months = 0:6:162; % 13.5 years * 12 months/year
interval_dates = t_start_cruise + calmonths(interval_months);
valid_indices = interval_dates <= t_end_orbit;
interval_dates = interval_dates(valid_indices);
% Find closest trajectory indices
interval_idx = zeros(size(interval_dates));
for i = 1:length(interval_dates)
    [~, idx] = min(abs(time_vec - interval_dates(i)));
    interval_idx(i) = idx;
end
% Create results table
results_table = table(...
    interval_dates',...
    dist_AU(interval_idx),...
    downlink_margin(interval_idx),...
    uplink_margin(interval_idx),...
    'VariableNames',...
    {'Date','Distance_AU','Downlink_Margin_dB','Uplink_Margin_dB'});
disp('Link Margins at 6-Month Intervals (Full Mission):');
disp(results_table);
% =========================================================================
% PLOT FULL MISSION MARGIN
% =========================================================================
time_years = days(time_vec - t_start_cruise)/365.25;
figure('Position', [100, 100, 1100, 600]);
hold on;
% Define a global font size for plots
global_font_size = 25; % Your requested global font size
grid_line_width = 2.5; % Your requested grid line width
curve_line_width = 4;  % Your requested curve line width
% Plot margins
plot(time_years, downlink_margin, 'b-', 'LineWidth', curve_line_width, 'DisplayName', 'Downlink');
plot(time_years, uplink_margin, 'r-', 'LineWidth', curve_line_width, 'DisplayName', 'Uplink');

% Set plot limits BEFORE calculating annotation coordinates
ylim_max = max(uplink_margin) + 5;
ylim([0, ylim_max]);
xlim([0, time_years(end)]);

% Annotate key milestones
[~, jupiter_idx] = min(abs(time_vec - datetime('2033-07-01')));
jupiter_x = time_years(jupiter_idx);
jupiter_y_downlink = downlink_margin(jupiter_idx);

% % Plot a dot marker at the intersection point
% plot(jupiter_x, jupiter_y_downlink, 'o', ... % 'o' is for circle (dot)
%      'MarkerSize', global_font_size * 0.8, ... % Adjust dot size
%      'MarkerFaceColor', 'g', ...
%      'MarkerEdgeColor', 'k', ...
%      'DisplayName', '', ...
%      'HandleVisibility', 'off'); % No legend entry for this marker

% Add text with arrow pointing down
text_x_offset = -2.8; % Adjusted: smaller value for horizontal offset (e.g., 0.8 years)
text_y_offset = 8;   % Adjusted: smaller value for vertical offset (e.g., 3 dB)
text_start_x = jupiter_x - text_x_offset; % Calculate the starting X for the text

text(text_start_x, jupiter_y_downlink + text_y_offset, ... % Offset text above and to the left of the marker
     sprintf('Jupiter Flyby - %.1f dB', jupiter_y_downlink), ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'Interpreter', 'latex', ...
     'FontSize', global_font_size * 0.8, ...
     'HandleVisibility', 'off');

% Add arrow pointing from text to the dot
% Get current axes properties AFTER setting limits
ax_pos = get(gca, 'Position'); % Get axes position within the figure [left, bottom, width, height]
x_lim = get(gca, 'XLim');      % Get current X limits [xmin, xmax]
y_lim = get(gca, 'YLim');      % Get current Y limits [ymin, ymax]

% Convert data coordinates to normalized figure coordinates for the arrow start and end points
% Arrow starts from below the text
arrow_start_x_norm = ax_pos(1) + (text_start_x - x_lim(1)) / (x_lim(2) - x_lim(1)) * ax_pos(3);
arrow_start_y_norm = ax_pos(2) + (jupiter_y_downlink + text_y_offset * 0.8 - y_lim(1)) / (y_lim(2) - y_lim(1)) * ax_pos(4); % Point slightly below text

% Arrow ends slightly above the dot
arrow_end_x_norm = ax_pos(1) + (jupiter_x - x_lim(1)) / (x_lim(2) - x_lim(1)) * ax_pos(3);
arrow_end_y_norm = ax_pos(2) + (jupiter_y_downlink + 0.1 - y_lim(1)) / (y_lim(2) - y_lim(1)) * ax_pos(4); % Point slightly above dot

annotation('arrow',...
    [arrow_start_x_norm arrow_end_x_norm], ...
    [arrow_start_y_norm arrow_end_y_norm],...
    'HeadStyle','vback2','HeadSize',10,'HeadWidth',15,'Color','k', ...
    'LineWidth', 2); % Added LineWidth property for a thicker line
xline(time_years(jupiter_idx), 'b:', 'DisplayName','Jupiter Flyby Time','LineWidth',curve_line_width*0.8); % Make xline blue



[~, neptune_idx] = min(abs(time_vec - t_end_cruise));
% text(time_years(neptune_idx), downlink_margin(neptune_idx)-2, ...
%      sprintf('Neptune Arrival\\newline%.1f dB', downlink_margin(neptune_idx)), ... % Use \newline for line break in LaTeX
%      'HorizontalAlignment','right', ...
%      'Interpreter', 'latex', ... % Set interpreter to latex
%      'FontSize', global_font_size * 0.8);
% Find the minimum downlink margin and its index
[min_downlink_margin, min_downlink_idx] = min(downlink_margin);
min_downlink_time_years = time_years(min_downlink_idx);
% Add vertical line for start of orbital phase
orbit_start_years = time_years(neptune_idx + 1);
xline(orbit_start_years, 'k:', 'LineWidth', curve_line_width*0.8, 'DisplayName', 'Start of Orbital Phase');
% Format plot
xlabel('Time Since Departure (years)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('Link Margin (dB)', 'Interpreter', 'latex', 'FontSize', 10);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
ax = gca; % Get current axes handle
ax.GridLineWidth = grid_line_width;
ax.FontSize = global_font_size; % Apply global font size to axes tick labels
ax.TickLabelInterpreter = 'latex'; % Ensure tick labels also use LaTeX

% Adjust ylim to start at 0 dB, or slightly below if any data points are there
ylim_min_display = min(3, min(downlink_margin)); % Ensure it includes at least 0, or lower if margins go very negative
ylim_max = max(uplink_margin) + 5;
ylim([0, ylim_max]);
xlim([0, time_years(end)+0.5]);
xline(time_years(end), 'm:', 'DisplayName','End of mission','LineWidth',curve_line_width*0.8); % Make xline blue

% % Highlight positive margin region (green, darker)
% patch([min(xlim) max(xlim) max(xlim) min(xlim)],...
%       [0 0 ylim_max ylim_max],...
%       [0 0.6 0], 'EdgeColor','none', 'FaceAlpha',0.3, 'DisplayName','Positive Margin'); % Darker green
% 
% % Highlight negative margin region (red, darker) - assumes negative margin is anything below 0 dB
% patch([min(xlim) max(xlim) max(xlim) min(xlim)],...
%       [ylim_min_display ylim_min_display 0 0],... % From min display limit up to 0 dB
%       [0.8 0 0], 'EdgeColor','none', 'FaceAlpha',0.3, 'DisplayName','Negative Margin'); % Darker red

hold off;
% =========================================================================
% HELPER FUNCTIONS (Copied/Adapted from Detailed Script)
% =========================================================================
% computeFSPL: Calculates Free Space Path Loss in dB
function FSPL_dB = computeFSPL(distance_m, freq_Hz)
    % distance_m: distance in meters
    % freq_Hz: frequency in Hertz
    c_val = 299792458; % Speed of light in m/s (consistent with the detailed script)
    lambda_val = c_val / freq_Hz;
    FSPL_dB = 20 * log10(4 * pi * distance_m / lambda_val);
end
% computeAntennaGain: Calculates antenna gain in dBi
function G_dBi = computeAntennaGain(D, eta, freq)
    % D: Antenna Diameter in meters
    % eta: Antenna Efficiency (0 to 1)
    % freq: Frequency in Hertz
    c_val = 299792458; % Speed of light in m/s
    lambda_val = c_val / freq;
    G_linear = eta * (pi * D / lambda_val)^2;
    G_dBi = 10 * log10(G_linear);
end
% computeEIRP: Calculates Equivalent Isotropically Radiated Power in dBm
function EIRP_dBm = computeEIRP(P_tx_W, G_tx_dBi, tx_loss_dB)
    % P_tx_W: Transmit Power in Watts
    % G_tx_dBi: Transmitter Antenna Gain in dBi
    % tx_loss_dB: Total Transmitter Line Losses in dB (should be a positive value, subtracted)
    P_tx_dBm = 10 * log10(P_tx_W) + 30; % Convert Watts to dBm
    EIRP_dBm = P_tx_dBm + G_tx_dBi - tx_loss_dB;
end
% computeEbN0: Calculates Energy per Bit to Noise Power Spectral Density (Eb/N0) in dB
function EbN0_dB = computeEbN0(P_rx_dBm, R_bps, T_sys_K, k_boltzmann)
    % P_rx_dBm: Received Power in dBm
    % R_bps: Data Rate in bits per second
    % T_sys_K: System Noise Temperature in Kelvin
    % k_boltzmann: Boltzmann's constant in J/K
    
    % Noise spectral density (N0) in dBW/Hz
    No_dB_per_Hz = 10*log10(k_boltzmann) + 10*log10(T_sys_K);
    
    % Convert received power from dBm to dBW
    P_rx_dBW = P_rx_dBm - 30;
    
    % Data rate in dB-Hz
    R_dBHz = 10*log10(R_bps);
    
    % Eb/N0 in dB
    EbN0_dB = P_rx_dBW - No_dB_per_Hz - R_dBHz;
end
% computeBeamwidth: Calculates 3dB beamwidth in degrees for a parabolic antenna
function beamwidth_deg = computeBeamwidth(D, freq)
    % D: Antenna Diameter in meters
    % freq: Frequency in Hertz
    c_val = 299792458; % Speed of light in m/s
    lambda_val = c_val / freq;
    k_beamwidth = 70; % Constant for parabolic dish (approximate)
    if D <= 0
        beamwidth_deg = Inf; % Handle zero or negative diameter to avoid division by zero
    else
        beamwidth_deg = k_beamwidth * (lambda_val / D);
    end
end
% computePointingLoss: Calculates pointing loss in dB based on pointing error and beamwidth
function L_p_dB = computePointingLoss(pointing_error_deg, beamwidth_deg)
    % pointing_error_deg: Pointing error in degrees
    % beamwidth_deg: 3dB Beamwidth of the antenna in degrees
    if beamwidth_deg <= 0 || pointing_error_deg < 0
        L_p_dB = Inf; % Invalid inputs
    elseif pointing_error_deg > (beamwidth_deg * 1.5) % Beyond a certain point, loss becomes very high
        L_p_dB = 20; % Cap the loss for extreme errors (arbitrary cap, but reasonable)
    else
        % Common approximation for pointing loss for small errors relative to beamwidth
        L_p_dB = 12 * (pointing_error_deg / beamwidth_deg)^2;
    end
end
% ITU ATTENUATION FUNCTIONS (from Detailed Script)
% computeRainAttenuationITU: Calculates rain attenuation using a simplified ITU-R P.618 model.
function A_rain_dB = computeRainAttenuationITU(freq, elevation_angle_deg, rain_rate_mm_hr, polarization_tilt_angle_deg, site_latitude_deg)
    f_GHz = freq / 1e9;
    if f_GHz < 1 || f_GHz > 1000
        A_rain_dB = 0;
        warning('ITU-R P.618 model valid for 1-1000 GHz. Setting rain attenuation to 0 for %.2f GHz.', f_GHz);
        return;
    end
    if elevation_angle_deg < 0.1
        A_rain_dB = 1000; % Very high loss for extremely low angles
        warning('Elevation angle too low (<0.1 deg) for practical rain attenuation model. Setting to high loss.');
        return;
    end
    % ITU-R P.838-3 Coefficients for X-band (8.4 GHz, 20 deg C) - Approximate
    k_H = 0.0188; alpha_H = 1.050;
    k_V = 0.0173; alpha_V = 1.030;
    
    k_eff = (k_H + k_V + (k_H - k_V) * cosd(2 * polarization_tilt_angle_deg) * cosd(elevation_angle_deg)) / 2;
    alpha_eff = (k_H * alpha_H + k_V * alpha_V + (k_H * alpha_H - k_V * alpha_V) * cosd(2 * polarization_tilt_angle_deg) * cosd(elevation_angle_deg)) / (2 * k_eff);
    
    gamma_R = k_eff * (rain_rate_mm_hr^alpha_eff); % Specific attenuation (dB/km)
    
    % ITU-R P.839-4 Rain Height (0-degree isotherm height) model - Simplified
    abs_lat = abs(site_latitude_deg);
    if abs_lat <= 36
        H_R = 3.6; % km (average for tropical/temperate latitudes)
    else
        H_R = 3.6 - 0.075 * (abs_lat - 36);
        if H_R < 1.0 % Ensure minimum plausible rain height
            H_R = 1.0;
        end
    end
    L_S_km = (H_R - 0) / sind(elevation_angle_deg); % Slant path length through rain (assuming ground height = 0km)
    
    % Path reduction factor (Simplified ITU-R P.618 approach)
    if L_S_km > 10
        r_p = 1 / (1 + (L_S_km / 35));
    else
        r_p = 1;
    end
    effective_rain_path_km = r_p * L_S_km;
    A_rain_dB = gamma_R * effective_rain_path_km; % Total rain attenuation (dB)
end
% computeCloudAttenuationITU: Calculates cloud attenuation using a simplified ITU-R P.840 model.
function A_cloud_dB = computeCloudAttenuationITU(freq, cloud_temp_celsius, liquid_water_density_gm3)
    f_GHz = freq / 1e9;
    if f_GHz < 10 || f_GHz > 1000
        A_cloud_dB = 0;
        warning('ITU-R P.840 model valid for 10-1000 GHz. Setting cloud attenuation to 0 for %.2f GHz.', f_GHz);
        return;
    end
    % T_K = cloud_temp_celsius + 273.15; % Not directly used in this simplified k_L_approx
    k_L_approx = 0.008 * (f_GHz/8.4)^1.5; % Simplified specific attenuation coefficient
    if f_GHz > 100
        k_L_approx = k_L_approx * 2; % Ad-hoc adjustment from detailed script for higher frequencies
    end
    gamma_cloud_dB_per_gm3_km = k_L_approx;
    specific_atten_dB_km = gamma_cloud_dB_per_gm3_km * liquid_water_density_gm3;
    effective_cloud_path_km = 1.0; % Assuming 1 km effective path for simplicity
    A_cloud_dB = specific_atten_dB_km * effective_cloud_path_km;
end
% computeGasAttenuationITU: Calculates atmospheric gas attenuation (oxygen and water vapor)
function A_gas_dB = computeGasAttenuationITU(freq, air_temp_celsius, atmospheric_pressure_hPa, water_vapor_density_gm3, elevation_angle_deg)
    f_GHz = freq / 1e9;
    if f_GHz < 1 || f_GHz > 1000
        A_gas_dB = 0;
        warning('ITU-R P.676 model valid for 1-1000 GHz. Setting gas attenuation to 0 for %.2f GHz.', f_GHz);
        return;
    end
    if elevation_angle_deg < 0.1
        A_gas_dB = 1000; % Very high loss for extremely low angles
        warning('Elevation angle too low (<0.1 deg) for practical gas attenuation model. Setting to high loss.');
        return;
    end
    % Simplified Zenith Attenuation coefficients
    zenith_oxygen_att_dB = 0.02; % Approximate for X-band
    zenith_water_vapor_att_dB = 0.1 * (water_vapor_density_gm3 / 7.5); % Scaled based on reference density
    total_zenith_att_dB = zenith_oxygen_att_dB + zenith_water_vapor_att_dB;
    A_gas_dB = total_zenith_att_dB / sind(elevation_angle_deg); % Slant path scaling
end
% coe2rv: Converts classical orbital elements (COE) to Cartesian state vectors (r, v).
function [r, v] = coe2rv(a, e, i, RAAN, omega, nu, mu)
    % Converts classical orbital elements (COE) to Cartesian state vectors (r, v).
    % All angles should be in radians as per standard definition, but this function
    % converts degrees input to radians for calculation.
    i = deg2rad(i);
    RAAN = deg2rad(RAAN);
    omega = deg2rad(omega);
    
    p = a * (1 - e^2); % Semi-latus rectum
    r_norm = p / (1 + e * cos(nu));
    
    % Position and velocity in perifocal frame
    r_pqw = [r_norm * cos(nu); r_norm * sin(nu); 0];
    v_pqw = [-sqrt(mu/p) * sin(nu); sqrt(mu/p) * (e + cos(nu)); 0];
    
    % Rotation matrices (from perifocal to inertial frame)
    R3_RAAN = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_omega = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
    
    % Combined rotation from perifocal to inertial frame (assuming Z-X-Z sequence)
    Q_pqw2ijk = R3_RAAN' * R1_i' * R3_omega'; % Inverse rotations applied in reverse order
    
    % State vectors in inertial frame
    r = Q_pqw2ijk * r_pqw;
    v = Q_pqw2ijk * v_pqw;
end