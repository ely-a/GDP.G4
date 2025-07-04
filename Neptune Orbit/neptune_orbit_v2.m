clear
clc
close all

%% SET GRAPHING FONTS AND SIZES
set(groot,'defaultLineLineWidth',4) % Set all line widths to 2 unless stated otherwise
set(groot,'defaultAxesFontSize',34) % Set all axes font size to 16 unless stated otherwise
set(groot,'defaulttextfontsize',36) % Set all text font sizes to 16 unless stated otherwise
set(groot,'defaultLineMarkerSize',12) % Set all marker sizes to 16 unless stated otherwise
set(groot,'defaultAxesXGrid','on') % Turn x grid lines on
set(groot,'defaultAxesYGrid','on') % Turn y grid lines on

%% Orbit around Neptune 

% Neptune Orbital Trajectory Design
% Author: Team 04 - Trajectory Design
% Date: May 2025

%% === CONSTANTS ===
mu_S = 1.32712440018e11;      % Sun gravitational parameter [km^3/s^2]
mu_N = 6.8351e6;              % Neptune gravitational parameter [km^3/s^2]
R_N = 24764;                  % Neptune radius [km]
R_T = 13526;

%% === TRITON PROPERTIES ===
if exist('moon_orbits.mat', 'file') ~= 2 
    load_moon_data()
end 
load('moon_orbits.mat','orbits_mat','moon_names', 'max_len', 'times_all', 'velocities_mat');
triton_orbit = [orbits_mat(:, :, 14),velocities_mat(:,:,14)];
triton_times = times_all{14};


arrival_epoch = 14762.860319752668; % MJD2000
arrival_epoch = arrival_epoch + 2451545.0; 

interpolated_time = interp1(triton_times, triton_orbit, arrival_epoch, 'spline');

r_Triton_0 = interpolated_time(1, 1:3);
v_Triton_0 = interpolated_time(1, 4:6);

[a_triton, e_triton, h_triton, i_triton, Omega_triton, omega_triton, ~] = oe_from_rv(r_Triton_0, v_Triton_0, mu_N);

%% Get initial spacecraft orbit

% === HIGHLY ELLIPTICAL ORBIT ===
% BENS VECTORS 
% r_sc_initial = [1.9793e+04; 1.5867e+04; 0];
% v_sc_initial = [1.2467e+01; -1.2873e+01; 0];

ra_cap = 5.634e6;           % Apoapsis [km]
rp_cap = R_N + 10000;       % Periapsis [km] (initial guess - can change now - no more aerocapture - doesnt change much if 2000km ect)
a_cap = (rp_cap + ra_cap)/2;
e_cap = (ra_cap - rp_cap)/(ra_cap + rp_cap);
h_cap = sqrt(mu_N * rp_cap * (1 + e_cap));
i_cap = 180 - i_triton;
Omega_cap = Omega_triton-180;
omega_cap = 90;
[r0_sc, v0_sc] = rv_from_oe(a_cap, e_cap, i_cap, Omega_cap, omega_cap, 0, mu_N);
P_NOW = 2 * pi * a_cap^(3/2) / sqrt(mu_N);

save("vectors.mat", "r0_sc", "v0_sc")

figure 
hold on
plotOrbit(a_triton, e_triton, i_triton, Omega_triton, omega_triton, mu_N, 'Triton')
plotOrbit(a_cap, e_cap, i_cap, Omega_cap, omega_cap, mu_N, 'spacecraft')
axis equal
legend()


%% Triton at apogee

P = 2 * pi * sqrt(a_cap^3/mu_N);
tof_a = P/2; 

ap_epoch = arrival_epoch + tof_a / (3600 * 24);

statevector = interp1(triton_times, triton_orbit, arrival_epoch, 'spline');

r_Triton_ap = statevector(1, 1:3);
v_Triton_ap = statevector(1, 4:6);

%% 

rp_try = linspace(3e4, 3.6e5, 1000);

for k = 1 : length(rp_try)
    rp_int = rp_try(k);
    a_int = (rp_int + ra_cap)/2; % so this changes to ra_initial 
    e_int = (ra_cap - rp_int)/(ra_cap + rp_int);
    h_int = sqrt(mu_N * rp_int * (1 + e_int));

    theta_int = acosd((h_int^2 / (mu_N*a_triton) - 1)/e_int);
    if theta_int <=180 
        theta_int = 360-theta_int;
    end

    E = 2 * atan( sqrt((1 - e_int)/(1 + e_int)) * tand(theta_int/2) );
    M = E - e_int * sin(E);
    if M < 0
        M = M + 2*pi;
    end
    P = 2 * pi * sqrt(a_int^3 / mu_N);
    T = M * P/(2*pi) - P/2;

    k;

    r_sc_int = rv_from_oe(a_int, e_int, i_cap, Omega_cap, omega_cap, theta_int, mu_N);

    int_epoch = ap_epoch + T / (3600 * 24);

    statevector = interp1(triton_times, triton_orbit, int_epoch, 'spline');

    r_Triton_int = statevector(1, 1:3)';


    if norm(r_Triton_int - r_sc_int) <= R_T + 5000
        rp_int;
        dv_int = h_cap / ra_cap-h_int/ra_cap;
        break 
    end 

end 



%% during flyby 

mu_T = 1428.495;

[~, v_sc_int] = rv_from_oe(a_int, e_int, i_cap, Omega_cap, omega_cap, theta_int, mu_N);
v_Triton_int = statevector(1, 4:6)';

v_inf_in = v_sc_int - v_Triton_int;
v_inf_in_mag = norm(v_inf_in)

a_flyby = mu_T / v_inf_in_mag^2; 
rp_flyby = R_T + 1000;
%e_flyby = rp_flyby/a_flyby + 1;
e_flyby = 1 + ((rp_flyby * v_inf_in_mag^2) / mu_T) % tf????? 1047 ????
delta_flyby = 2*asin(1/e_flyby);
h_flyby = sqrt(rp_flyby * mu_T * (1 + e_flyby));

v_hat_in = v_inf_in / v_inf_in_mag;

% Step 2: Rotation plane (find perpendicular unit vector)
h_vec = cross(v_inf_in, [0; 0; 1]);   % arbitrary axis if motion isn't planar

n_hat = cross(h_vec, v_hat_in);
n_hat = n_hat / norm(n_hat);

% Step 3: Rotate v_inf_in by turn angle in the plane
v_inf_out = v_inf_in_mag * (cos(delta_flyby) * v_hat_in + sin(delta_flyby) * n_hat);

% Step 4: Back to inertial frame
v_sc_after = v_inf_out + v_Triton_int;

delta_v_flyby = norm(v_sc_after) - norm(v_sc_int)

angle_deg = acosd(dot(v_sc_after, v_sc_int) / (norm(v_sc_after)*norm(v_sc_int)))


% angle between v in and out needs to be 0 - 60 ?


%%

% % === MULTIPLE FLYBYS LOOP ===
% num_flybys = 5;
% r_sc_now = r_sc_int;
% v_sc_now = v_sc_after;
% epoch_now = int_epoch;
% 
% e_list = NaN(1, num_flybys);
% a_list = NaN(1, num_flybys);
% flyby_epoch_list = NaN(1, num_flybys);
% 
% e_prev = e_int;  % start with high value otherwise 1 
% 
% for n = 1:num_flybys
% 
%     [a_now, e_now, h_now, i_now, Omega_now, omega_now, theta_now] = find_OE(r_sc_now, v_sc_now, mu_N);
% 
%     % Check for divergence
%     if e_now > 0.995
%         warning("Eccentricity diverging at flyby %d (e = %.4f). Aborting loop.", n, e_now);
%         break;
%     end
% 
%     e_list(n) = e_now;
%     a_list(n) = a_now;
%     flyby_epoch_list(n) = epoch_now;
% 
%     % Skip if eccentricity increases
%     if n > 1 && e_now > e_prev
%         fprintf("❌ Flyby %d rejected (e = %.4f > %.4f)\n", n, e_now, e_prev);
%         continue;
%     end
%     e_prev = e_now;
% 
%     % Time to next apoapsis
%     P_now = 2 * pi * sqrt(a_now^3 / mu_N);
%     T_to_ap = P_now / 2;
%     epoch_ap = epoch_now + T_to_ap / (3600 * 24);
% 
%     flyby_found = false;
% 
%     % Try various periapses to find next intersection
%     for k = 1:length(rp_try)
%         rp_next = rp_try(k);
%         a_next = (rp_next + ra_cap)/2;
%         e_next = (ra_cap - rp_next)/(ra_cap + rp_next);
%         h_next = sqrt(mu_N * rp_next * (1 + e_next));
% 
%         cos_theta = (h_next^2 / (mu_N * a_triton) - 1)/e_next;
%         if abs(cos_theta) > 1, continue; end
% 
%         theta_next = acosd(cos_theta);
%         if theta_next <= 180
%             theta_next = 360 - theta_next;
%         end
% 
%         E = 2 * atan( sqrt((1 - e_next)/(1 + e_next)) * tand(theta_next/2) );
%         M = E - e_next * sin(E);
%         if M < 0
%             M = M + 2*pi;
%         end
%         %if ~isreal(M), continue; end
% 
%         P = 2 * pi * sqrt(a_next^3 / mu_N);
%         T = M * P/(2*pi) - P/2;
%         epoch_int = epoch_ap + T / (3600 * 24);
%         if ~isreal(epoch_int) || isnan(epoch_int), continue; end
% 
%         statevector = interp1(triton_times, triton_orbit, epoch_int, 'spline');
%         r_Triton = statevector(1, 1:3)';
%         v_Triton = statevector(1, 4:6)';
% 
%         r_sc = rv_from_oe(a_next, e_next, i_now, Omega_now, omega_now, theta_next, mu_N);
%         flyby_window = R_T + 5000; %was 15000
%         if norm(r_Triton - r_sc) <= flyby_window
%             flyby_found = true;
%             break;
%         end
%     end
% 
%     if ~flyby_found
%         warning("No flyby opportunity found at step %d.", n);
%         break;
%     end
% 
%     % === Perform Flyby ===
%     [~, v_sc_before] = rv_from_oe(a_next, e_next, i_now, Omega_now, omega_now, theta_next, mu_N);
%     v_inf = v_sc_before - v_Triton;
%     v_inf_mag = norm(v_inf);
% 
%     % Reduce rp over time to increase energy exchange
%     rp_flyby = max(R_T + 500, R_T + 5000 - n * 150);
%     a_flyby = mu_T / v_inf_mag^2;
%     e_flyby = 1 + ((rp_flyby * v_inf_mag^2) / mu_T);
%     delta = 2 * asin(1 / e_flyby);
% 
%     v_hat_in = v_inf / v_inf_mag;
%     h_vec = cross(v_inf, [0; 0; 1]);
%     if abs(dot(v_hat_in, [1;0;0])) < 0.95
%     ref_axis = [1;0;0];
%     else
%         ref_axis = [0;1;0];
%     end
%     h_vec = cross(v_inf, ref_axis);
%     n_hat = cross(h_vec, v_inf);
%     n_hat = n_hat / norm(n_hat);
% 
%     %n_hat = cross(h_vec, v_hat_in); 
%     % n_hat = n_hat / norm(n_hat);
% 
%     % REVERSED sign for velocity deflection
%     v_inf_out = v_inf_mag * (cos(delta) * v_hat_in - sin(delta) * n_hat);
% 
%     v_sc_now = v_inf_out + v_Triton;
%     r_sc_now = r_sc;
%     epoch_now = epoch_int;
% 
%     fprintf("✅ Flyby %d complete: e = %.4f, rp_flyby = %.0f km\n", n, e_now, rp_flyby - R_T);
% end

% % === MULTIPLE FLYBYS LOOP (ENERGY-TOLERANT, LOGGING ENHANCED) ===
% num_flybys = 5;
% r_sc_now = r_sc_int;
% v_sc_now = v_sc_after;
% epoch_now = int_epoch;
% 
% e_list = NaN(1, num_flybys);
% a_list = NaN(1, num_flybys);
% flyby_epoch_list = NaN(1, num_flybys);
% delta_list = NaN(1, num_flybys);
% 
% % Initial energy
% E_prev = norm(v_sc_now)^2 / 2 - mu_N / norm(r_sc_now);
% 
% for n = 1:num_flybys
%     % Get current orbital elements
%     [a_now, e_now, h_now, i_now, Omega_now, omega_now, theta_now] = find_OE(r_sc_now, v_sc_now, mu_N);
% 
%     % Check for divergence
%     if e_now > 0.995
%         warning("Eccentricity diverging at flyby %d (e = %.4f). Aborting loop.", n, e_now);
%         break;
%     end
% 
%     e_list(n) = e_now;
%     a_list(n) = a_now;
%     flyby_epoch_list(n) = epoch_now;
% 
%     % Time to apoapsis
%     P_now = 2 * pi * sqrt(a_now^3 / mu_N);
%     T_to_ap = P_now / 2;
%     epoch_ap = epoch_now + T_to_ap / (3600 * 24);
% 
%     flyby_found = false;
% 
%     for k = 1:length(rp_try)
%         rp_next = rp_try(k);
%         a_next = (rp_next + ra_cap)/2;
%         e_next = (ra_cap - rp_next)/(ra_cap + rp_next);
%         h_next = sqrt(mu_N * rp_next * (1 + e_next));
% 
%         cos_theta = (h_next^2 / (mu_N * a_triton) - 1)/e_next;
%         if abs(cos_theta) > 1, continue; end
% 
%         theta_next = acosd(cos_theta);
%         if theta_next <= 180
%             theta_next = 360 - theta_next;
%         end
% 
%         E = 2 * atan( sqrt((1 - e_next)/(1 + e_next)) * tand(theta_next/2) );
%         M = E - e_next * sin(E);
%         if M < 0, M = M + 2*pi; end
% 
%         P = 2 * pi * sqrt(a_next^3 / mu_N);
%         T = M * P/(2*pi) - P/2;
%         epoch_int = epoch_ap + T / (3600 * 24);
%         if ~isreal(epoch_int) || isnan(epoch_int), continue; end
% 
%         statevector = interp1(triton_times, triton_orbit, epoch_int, 'spline');
%         r_Triton = statevector(1, 1:3)';
%         v_Triton = statevector(1, 4:6)';
% 
%         r_sc = rv_from_oe(a_next, e_next, i_now, Omega_now, omega_now, theta_next, mu_N);
%         flyby_window = R_T + 5000;
%         if norm(r_Triton - r_sc) <= flyby_window
%             flyby_found = true;
%             break;
%         end
%     end
% 
%     if ~flyby_found
%         warning("No flyby opportunity found at step %d.", n);
%         break;
%     end
% 
%     % === Perform Flyby ===
%     [~, v_sc_before] = rv_from_oe(a_next, e_next, i_now, Omega_now, omega_now, theta_next, mu_N);
%     v_inf = v_sc_before - v_Triton;
%     v_inf_mag = norm(v_inf);
% 
%     % Refine flyby periapsis and angle
%     rp_flyby = max(R_T + 500, R_T + 5000 - n * 150);
%     a_flyby = mu_T / v_inf_mag^2;
%     e_flyby = 1 + ((rp_flyby * v_inf_mag^2) / mu_T);
%     delta = 2 * asin(1 / e_flyby);
%     delta_list(n) = delta;
% 
%     % Build turn frame
%     if abs(dot(v_inf, [1; 0; 0]) / v_inf_mag) < 0.95
%         ref_axis = [1; 0; 0];
%     else
%         ref_axis = [0; 1; 0];
%     end
%     h_vec = cross(v_inf, ref_axis);
%     n_hat = cross(h_vec, v_inf); n_hat = n_hat / norm(n_hat);
%     v_hat_in = v_inf / v_inf_mag;
% 
%     % Rotate
%     v_inf_out = v_inf_mag * (cos(delta) * v_hat_in - sin(delta) * n_hat);
%     v_sc_now = v_inf_out + v_Triton;
%     r_sc_now = rv_from_oe(a_next, e_next, i_now, Omega_now, omega_now, theta_next, mu_N);
%     epoch_now = epoch_int;
% 
%     % Energy check
%     E_now = norm(v_sc_now)^2 / 2 - mu_N / norm(r_sc_now);
%     dE = E_now - E_prev;
% 
%     if E_now > E_prev - 1e-3
%         fprintf("❌ Flyby %d rejected (ΔE = %.4f, E = %.4f > %.4f)\n", ...
%                 n, dE, E_now, E_prev);
%         continue;
%     else
%         [~, e_now_new, ~, ~, ~, ~, ~] = find_OE(r_sc_now, v_sc_now, mu_N);
%         fprintf("✅ Flyby %d complete: e = %.4f, rp_flyby = %.0f km, ΔE = %.4f\n", ...
%                 n, e_now_new, rp_flyby - R_T, dE);
%         E_prev = E_now;
%     end
% end


% === MULTIPLE TRITON FLYBYS LOOP ===
num_flybys = 5;
r_sc_now = r_sc_int;
v_sc_now = v_sc_after;
epoch_now = int_epoch;

mu_T = 1428.495;

flyby_angles = zeros(1, num_flybys);
flyby_epochs = zeros(1, num_flybys);

for flyby_num = 1:num_flybys

    % Get orbital elements from current state
    [a_now, e_now, h_now, i_now, Omega_now, omega_now, theta_now] = oe_from_rv(r_sc_now, v_sc_now, mu_N);
    P_now = 2 * pi * sqrt(a_now^3 / mu_N);
    
    % Estimate time of next apoapsis
    tof_next_ap = P_now / 2; % half-period
    epoch_next = epoch_now + tof_next_ap / (3600 * 24);

    % Interpolate Triton state at this time
    statevector = interp1(triton_times, triton_orbit, epoch_next, 'spline');
    r_Triton = statevector(1,1:3)';
    v_Triton = statevector(1,4:6)';

    % Propagate SC to this point
    theta_ap = 180;
    r_sc_next = rv_from_oe(a_now, e_now, i_now, Omega_now, omega_now, theta_ap, mu_N);
    [~, v_sc_next] = rv_from_oe(a_now, e_now, i_now, Omega_now, omega_now, theta_ap, mu_N);

    % Check for encounter
    if norm(r_sc_next - r_Triton) > R_T + 5000
        warning("No encounter at flyby #%d", flyby_num);
        break;
    end

    % Compute v-infinity and flyby deflection
    v_inf_in = v_sc_next - v_Triton;
    v_inf_mag = norm(v_inf_in);
    a_flyby = mu_T / v_inf_mag^2;
    rp_flyby = R_T + 1000;
    e_flyby = 1 + (rp_flyby * v_inf_mag^2 / mu_T);
    delta = 2 * asin(1 / e_flyby);
    
    % Rotate velocity vector
    v_hat_in = v_inf_in / v_inf_mag;
    h_vec = cross(v_inf_in, [1; 0; 0]);
    if norm(h_vec) < 1e-6
        h_vec = cross(v_inf_in, [0; 1; 0]);
    end
    n_hat = cross(h_vec, v_hat_in); n_hat = n_hat / norm(n_hat);
    v_inf_out = v_inf_mag * (cos(delta) * v_hat_in + sin(delta) * n_hat);
    
    v_sc_out = v_inf_out + v_Triton;

    % Store and update
    flyby_angles(flyby_num) = acosd(dot(v_sc_next, v_sc_out) / (norm(v_sc_next) * norm(v_sc_out)));
    flyby_epochs(flyby_num) = epoch_next;
    
    % Prepare for next iteration
    r_sc_now = r_sc_next;
    v_sc_now = v_sc_out;
    epoch_now = epoch_next;

    fprintf("Flyby %d complete. Turn angle: %.2f deg\n", flyby_num, flyby_angles(flyby_num));

end




%% Post flyby - start of 5 years 

% imma explain what i did here before i forget, so worked out orbital
% elements after flyby, then worked out r,v at apoapsis (180deg) of this orbit then
% used that r to set the new orbits apoapsis and calculated the velocity at that same point (180deg)
% hence found dv at apoapsis to set it into new trajectory of rp = 2000 km
% + r_n to meet instruments 
%then from here i can repeat same thing and do a dv at periapsis maybe this
%time to change into new trajectory and then do plane changes after this?

[a_return, e_return, h_return, i_return, Omega_return, omega_return, theta_return] = oe_from_rv(r_sc_int, v_sc_after, mu_N);
[rp_return, vp_return] = rv_from_oe(a_return, e_return, i_return, Omega_return, omega_return, 180 , mu_N);  % at apoapsis 

ra_main1 = norm(rp_return);    % then set this ra_main1 as rp_return at 180 deg  - then change to 1.2e5    
rp_main1 = R_N + 2000;      % adjust  
a_main1 = (rp_main1 + ra_main1)/2;
e_main1 = (ra_main1 - rp_main1)/(ra_main1 + rp_main1);
h_main1 = sqrt(mu_N * rp_main1 * (1 + e_main1));
i_main1 = i_return;
Omega_main1 = Omega_return;
omega_main1 = omega_return;
[r0_main1, v0_main1] = rv_from_oe(a_main1, e_main1, i_main1, Omega_main1, omega_main1, 180, mu_N);

dv_return = norm(v0_main1 - vp_return); % dv at apoapsis 

P_return = 2 * pi * a_main1^(3/2) / (sqrt(mu_N));
P_return_days = P_return / (3600 * 24);


time_flyby_inst = find_time_flyby(R_T, h_flyby, e_flyby, a_flyby, mu_T, 50000);

time_flyby_inst * 60;

%% instgrumetns around neptune 


ra_test = norm(rp_return);    % then set this ra_main1 as rp_return at 180 deg  - then change to 1.2e5    
rp_test = R_N + 10000;      % adjust  
a_test = (rp_test + ra_test)/2;
e_test = (ra_test - rp_test)/(ra_test + rp_test);
h_test = sqrt(mu_N * rp_test * (1 + e_test));
time_new = find_time_neptune(R_N, h_test, e_test, a_test, mu_N, 50000);



%find theta for the perigeemain position then for the current flyby
%position we have? then do time for those and take away then add onto int
%epoch 

%theta where we are and then at apog 

%theta at flyby = theta_return  - finding time from periapsis to apoapsis 
[a_backtoapo, e_backtoapo, h_backtoapo, i_backtoapo, ~, ~, theta_backtoapo] = oe_from_rv(rp_return, vp_return, mu_N);


P_backtoapo = 2 * pi * sqrt(a_return^3 / mu_N);
T_backtoapo = P_backtoapo/2;
%T_backtoapo = M_backtoapo * P_backtoapo/(2*pi);  

%finding time from periapsis for where we are in flyby 
theta_return = 360 - theta_return; % so this will be time from perigee to down to initial flyby pos to reach 
E_return = 2 * atan( sqrt((1 - e_return)/(1 + e_return)) * tand(theta_return/2) );
M_return = E_return - e_return * sin(E_return);
    if M_return < 0
        M_return = M_return + 2*pi;
    end
P_return = 2 * pi * sqrt(a_return^3 / mu_N);
T_return = M_return * P_return/(2*pi);


time_taken = T_backtoapo + T_return; % so adding the 180deg from perig to apog + flyby to perig 

epoch_after = int_epoch + time_taken / (3600 * 24);

%% After returning to to neptune epoch 

%epoch_main1 = epoch_after + time_taken2 / (3600 * 24);


%% final orbit - hopefully?

ra_final = R_N + 4e4;   
rp_final = R_N + 1500; % can be less but should be a bit less than 2000 km for maximum like max 1990 ish so atleast 2 mins for smallest instruemnt
a_final = (rp_final + ra_final)/2;
e_final = (ra_final - rp_final)/(ra_final + rp_final);
h_final = sqrt(mu_N * rp_final * (1 + e_final));
i_final = i_main1;
Omega_final = Omega_main1;
omega_final = omega_main1;
%[r0_final, v0_final] = rv_from_oe(a_final, e_final, i_final, Omega_final, omega_final, 180, mu_N);

%dv_return = norm(v0_main1 - vp_return) % dv at apoapsis 

P_return = 2 * pi * a_final^(3/2) / (sqrt(mu_N));

timememem = find_time_neptune(R_N, h_final, e_final, a_final, mu_N, 2000);

save('ScienceOrbit_0','ra_final','rp_final','a_final','e_final','P_return')


%
%% Plotting interception

epochs = linspace(arrival_epoch, epoch_after, 1000);
dt = (epochs(2) - epochs(1)) * 24 * 3600;
rs_Triton = zeros(3, length(epochs));
rs_sc = zeros(3, length(epochs));
theta_sc = 0;

for j = 1:length(epochs)
    statevector = interp1(triton_times, triton_orbit, epochs(j), 'spline');
    rs_Triton(:,j) = statevector(1, 1:3)';

    if epochs(j) < ap_epoch
        rs_sc(:,j) = rv_from_oe(a_cap, e_cap, i_cap, Omega_cap, omega_cap, theta_sc, mu_N);
        E0 = 2 * atan( sqrt((1 - e_cap)/(1 + e_cap)) * tand(theta_sc/2) );
        M0 = E0 - e_cap * sin(E0);
        M1 = M0 + sqrt(mu_N/a_cap^3) * dt;
        E1 = solve_kepler(M1, e_cap, 1e-8, 1000);
        theta_sc = 2 * atand(sqrt((1 + e_cap)/(1 - e_cap)) * tan(E1/2));
    elseif epochs(j) > ap_epoch && epochs(j) < int_epoch
        rs_sc(:,j) = rv_from_oe(a_int, e_int, i_cap, Omega_cap, omega_cap, theta_sc, mu_N);
        E0 = 2 * atan( sqrt((1 - e_int)/(1 + e_int)) * tand(theta_sc/2) );
        M0 = E0 - e_int * sin(E0);
        M1 = M0 + sqrt(mu_N/a_int^3) * dt;
        E1 = solve_kepler(M1, e_int, 1e-8, 1000);
        theta_sc = 2 * atand(sqrt((1 + e_int)/(1 - e_int)) * tan(E1/2));
    else %epochs(j) > int_epoch && epochs(j) < epoch_after
        rs_sc(:,j) = rv_from_oe(a_return, e_return, i_cap, Omega_cap, omega_cap, theta_sc, mu_N);
        E0 = 2 * atan( sqrt((1 - e_return)/(1 + e_return)) * tand(theta_sc/2) );
        M0 = E0 - e_return * sin(E0);
        M1 = M0 + sqrt(mu_N/a_return^3) * dt;
        E1 = solve_kepler(M1, e_return, 1e-8, 1000);
        theta_sc = 2 * atand(sqrt((1 + e_return)/(1 - e_return)) * tan(E1/2));
    end

end

visual_radius = 5 * R_N;  % Scaled-up radius for visibility
[xx, yy, zz] = sphere(50);
surf(visual_radius*xx, visual_radius*yy, visual_radius*zz, ...
     'FaceColor', [0.2 0.4 0.85], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.5, ...
     'DisplayName', 'Neptune');


%% Animate

%animate_two_trajectories(rs_sc, rs_Triton)


%%

function plotOrbit(a, e, i, Omega, omega, mu, name)
    thetas = 0:0.5:360;
    r = zeros(3, length(thetas));
    for j = 1:length(thetas)
        [r(:, j) , ~] = rv_from_oe(a, e, i, Omega, omega, thetas(j), mu);
    end

    plot3(r(1,:), r(2,:), r(3,:), DisplayName=name);  % unperturbed
end

function animate_two_trajectories(r_sc, r_triton)
    % === PARAMETERS ===
    R_N = 24764; % Neptune radius [km]
    visual_radius = 2.5 * R_N; % make Neptune appear larger
    trail_len = 80;            % number of trailing points for spacecraft

    % === SETUP FIGURE ===
    figure('Position', [100, 100, 1600, 900])
    hold on;
    grid on;
    axis equal;
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    title('Neptune System Trajectory');
    view(3);

    % === COLORS ===
    col_sc = [0.6, 0.2, 0.8];    % purple spacecraft
    col_triton = [0.9, 0.1, 0.1];   % orange Triton
    col_neptune = [0.2, 0.4, 0.85];   % blue Neptune


    % === PLOT NEPTUNE ===
    visual_radius = 3 * R_N;
    [xx, yy, zz] = sphere(50);
    surf(visual_radius*xx, visual_radius*yy, visual_radius*zz, ...
        'FaceColor', col_neptune, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.5, ...
        'DisplayName', 'Neptune');

    % === PLOT STATIC ORBITS ===
    plot3(r_triton(1,:), r_triton(2,:), r_triton(3,:), '-', ...
        'Color', col_triton, 'LineWidth', 2, 'DisplayName', 'Triton Orbit');
    plot3(r_sc(1,:), r_sc(2,:), r_sc(3,:), '--', ...
        'Color', col_sc, 'LineWidth', 2, 'DisplayName', 'Spacecraft Trajectory');
    h_tail = plot3(NaN, NaN, NaN, '-', ...
        'Color', col_sc * 0.8, 'LineWidth', 1.5, 'DisplayName', 'Trajectory Path');

    % === CREATE ANIMATED OBJECTS ===
    h_sc = plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', col_sc, 'MarkerEdgeColor', 'k', 'DisplayName', 'Spacecraft');
    h_triton = plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', col_triton, 'MarkerEdgeColor', 'k', 'DisplayName', 'Triton');
    

    legend('Location', 'bestoutside');

    % === ANIMATION LOOP ===
    N = size(r_sc, 2);
    % === CREATE VIDEO WRITER ===
v = VideoWriter('neptune_trajectory.mp4', 'MPEG-4');
v.Quality = 100;         % Max quality
v.FrameRate = 30;        % 30 frames per second
open(v);

% === ANIMATION LOOP ===
for k = 1:N
    % Animate spacecraft
    set(h_sc, 'XData', r_sc(1,k), 'YData', r_sc(2,k), 'ZData', r_sc(3,k));

    % Tail
    tail_range = max(1, k - trail_len):k;
    set(h_tail, 'XData', r_sc(1,tail_range), ...
                'YData', r_sc(2,tail_range), ...
                'ZData', r_sc(3,tail_range));

    % Triton
    set(h_triton, 'XData', r_triton(1,k), 'YData', r_triton(2,k), 'ZData', r_triton(3,k));

    drawnow;

    % Capture current frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% === CLOSE VIDEO FILE ===
close(v);
disp('🎥 Video saved as "neptune_trajectory.mp4".');

end




% function animate_two_trajectories(r1, r2)
%     % r1, r2 are 3xN matrices; columns are position vectors over time
% 
%     % Setup figure
%     figure;
%     hold on;
%     grid on;
%     axis equal;
%     xlabel('X [km]');
%     ylabel('Y [km]');
%     zlabel('Z [km]');
%     title('Trajectory Animation');
%     view(3);
% 
%     % Plot full trajectories as lines for context
%     plot3(r1(1,:), r1(2,:), r1(3,:), 'b--', 'DisplayName', 'Object 1 Path');
%     plot3(r2(1,:), r2(2,:), r2(3,:), 'r--', 'DisplayName', 'Object 2 Path');
% 
%     % Create animated points
%     h1 = plot3(r1(1,1), r1(2,1), r1(3,1), 'bo', 'MarkerSize', 8, 'DisplayName', 'Object 1');
%     h2 = plot3(r2(1,1), r2(2,1), r2(3,1), 'ro', 'MarkerSize', 8, 'DisplayName', 'Object 2');
%     legend;
% 
%     % Animate
%     N = size(r1, 2);
%     for k = 1:N
%         set(h1, 'XData', r1(1,k), 'YData', r1(2,k), 'ZData', r1(3,k));
%         set(h2, 'XData', r2(1,k), 'YData', r2(2,k), 'ZData', r2(3,k));
%         drawnow;
%         pause(0.001); % Adjust speed here
%     end
% end


% === Triton Aerobraking Simulation ===

% === PARAMETERS ===
Cd = 2.0;               % Drag coefficient
A = 20;                 % Area in m² (e.g. inflatable LOFTID-type)
m = 500;                % Spacecraft mass in kg
v_entry = 1500;         % Relative velocity at periapsis [m/s]
mu_Triton = 1428.495e9; % Triton GM [m³/s²]
R_Triton = 1353.4e3;    % Radius of Triton [m]
rp = R_Triton + 100e3;  % Periapsis altitude [m] (e.g., 100 km)
ra = R_Triton + 1000e3; % Apoapsis altitude [m]

rho0 = 1.2e-9;          % Surface density [kg/m³]
H = 8000;               % Scale height [m]

% === BALLISTIC COEFFICIENT ===
BC = m / (Cd * A);  % [kg/m²]

% === DENSITY AT PERIAPSIS ===
alt = linspace(0, 1000e3, 500);     % Altitude range [m]
rho = rho0 * exp(-alt / H);         % Atmospheric profile
rho_p = rho0 * exp(-(rp - R_Triton)/ H);  % Density at periapsis

% === DECELERATION ===
a_drag = 0.5 * rho * v_entry^2 * Cd * A / m;

% === CUMULATIVE ΔV PER PASS ===
% Assume drag acts only near periapsis: ΔV = a * Δt
% Simplify: use density at periapsis and short time duration (e.g., 60 s)
duration = 60;  % seconds spent in atmosphere near periapsis
a_peak = 0.5 * rho_p * v_entry^2 * Cd * A / m;
deltaV_per_pass = a_peak * duration;

% === ESTIMATE NUMBER OF PASSES ===
a_elliptical = (rp + ra) / 2;
v_periapsis = sqrt(mu_Triton * (2/rp - 1/a_elliptical));
v_circular = sqrt(mu_Triton / rp);
dv_required = v_periapsis - v_circular;
num_passes = ceil(dv_required / deltaV_per_pass);

% === PLOTS ===
figure;
subplot(2,1,1);
plot(alt/1e3, rho);
xlabel('Altitude [km]'); ylabel('Density [kg/m³]');
title('Triton Atmospheric Density Profile');
grid on;

subplot(2,1,2);
plot(alt/1e3, a_drag);
xlabel('Altitude [km]'); ylabel('Deceleration [m/s²]');
title('Aerobraking Deceleration vs Altitude');
grid on;

% === OUTPUT SUMMARY ===
fprintf('Ballistic Coefficient: %.2f kg/m²\n', BC);
fprintf('Density at Periapsis: %.2e kg/m³\n', rho_p);
fprintf('ΔV per Aerobraking Pass: %.4f m/s\n', deltaV_per_pass);
fprintf('ΔV required for Circularization: %.2f m/s\n', dv_required);
fprintf('Estimated Number of Passes: %d\n', num_passes);

% === Atmospheric Density Profile Comparison ===

% Altitude range (0 to 1000 km)
alt_km = linspace(0, 1000, 500);       % Altitude in km
alt_m = alt_km * 1e3;                  % Convert to meters

% Exponential density model: rho = rho0 * exp(-h/H)
rho_exp = @(h, rho0, H) rho0 * exp(-h / H);

% Parameters for each body
params.Triton.rho0 = 1.2e-9;  params.Triton.H = 8000;     % [kg/m^3, m]
params.Mars.rho0   = 0.02;    params.Mars.H   = 11100;
params.Pluto.rho0  = 1.0e-7;  params.Pluto.H  = 60000;

% Calculate density profiles
rho_Triton = rho_exp(alt_m, params.Triton.rho0, params.Triton.H);
rho_Mars   = rho_exp(alt_m, params.Mars.rho0, params.Mars.H);
rho_Pluto  = rho_exp(alt_m, params.Pluto.rho0, params.Pluto.H);

% Plotting
figure;
semilogy(alt_km, rho_Triton, 'b', 'DisplayName', 'Triton'); hold on;
semilogy(alt_km, rho_Mars, 'r', 'DisplayName', 'Mars');
semilogy(alt_km, rho_Pluto, 'g', 'DisplayName', 'Pluto');

xlabel('Altitude [km]');
ylabel('Atmospheric Density [kg/m^3]');
title('Atmospheric Density Profiles of Triton, Mars, and Pluto');
legend;
grid on;





%% 

function [t] = find_time_flyby(R, h, e, a, mu, altitude)

% need to work out when r = altitude + r_T , find theta from where we are
% then the time before and after , also check if under sphere of influence
% of triton ? got all a_flyby, e_flyby ect 

    r_need = R + altitude;  
    theta_need = acosd((h^2 / (mu*r_need) - 1)/e)
    
    F_need = 2 * atan( sqrt((e - 1)/(e + 1)) * tand(theta_need/2) );
    M_need = e * sin(F_need) - F_need;
    if M_need < 0
        M_need = M_need + 2*pi;
    end
    P_need = 2 * pi * sqrt(a^3 / mu);
    T_need = M_need * P_need/(2*pi);
    t = T_need * 2;
    t = t /(3600);


end 


function [t] = find_time_neptune(R, h, e, a, mu, altitude)

    r_need = R + altitude;  
    % theta_need = acosd((h^2 / (mu*r_need) - 1)/e)
    % 
    % E_need = 2 * atan( sqrt((1-e)/(1+e)) * tand(theta_need/2) );

    % Clamp cos(theta) to avoid NaNs due to floating point precision
    cos_theta = (h^2 / (mu * r_need) - 1) / e;
    cos_theta = max(-1, min(1, cos_theta));
    theta_rad = acos(cos_theta);  % radians

    % Eccentric anomaly (elliptical orbit)
    E_need = 2 * atan( sqrt((1 - e)/(1 + e)) * tan(theta_rad / 2) );

    M_need = E_need - e * sin(E_need);
    if M_need < 0
        M_need = M_need + 2*pi;
    end
    P_need = 2 * pi * sqrt(a^3 / mu);
    T_need = M_need * P_need/(2*pi);
    t = T_need * 2;
    %t = t /(3600);


end 



    








    
    






