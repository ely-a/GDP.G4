clear all
clc

% Preallocate arrays for speed
theta_array = linspace(0, 2*pi, 100);
num_angles = length(theta_array);
oe_matrix = zeros(num_angles, 7); % Columns: h, e, RA, i, w, TA, a

for k = 1:num_angles
    thetaB_deg = theta_array(k) * 180/pi;

    % approach conditions to neptune (heliocentric)
    r_sc_neptune = [3.7386e+09	2.2776e+09	8.3916e+08] - [3.7410e+09	2.2792e+09	8.3974e+08];
    v_sc_neptune = [12.3155	11.5601	4.4616] - [-2.3148    3.4722    1.8519];
    mu = 6.8365e6;
    v_inf_mag = 17.14;
    result = bplane(r_sc_neptune, v_sc_neptune, mu, thetaB_deg);
    vinf_in = v_sc_neptune; % km/s
    
    % position and velocity of periapsis
    % vinf_mag = norm(vinf_in);  % Magnitude of incoming v_infinity [km/s]
    rp = result.rp_km;               % Periapsis radius [km]
    rp = 34764;
    delta = result.turn_angle_deg * pi/180;        % Turning angle [rad]
    
    % B-plane unit vectors
    S = result.S_hat; % Unit vector of incoming V_inf
    R = result.R_hat;
    T = result.T_hat;
    
    % B-vector (in the B-plane)
    B = result.B_vec;
    B_hat = B / norm(B);       % Normalize
    
    % Angular momentum direction
    h = cross(B, vinf_in, 2); 
    h_hat = h / norm(h);
    
    % Periapsis direction (unit vector)
    p_hat = cross(h_hat, S);
    
    % Periapsis position vector
    r_p_vec = rp * p_hat;
    
    % Velocity magnitude at periapsis
    vp_mag = sqrt(v_inf_mag^2 + 2 * mu / rp);
    
    % Periapsis velocity vector (perpendicular to r_p_vec in orbital plane)
    v_p_hat = cross(h_hat, p_hat);
    v_p_vec = vp_mag * v_p_hat;
    
    % --- Output ---
    % fprintf('Periapsis position vector [km]: [%f, %f, %f]\n', r_p_vec);
    % fprintf('Periapsis velocity vector [km/s]: [%f, %f, %f]\n', v_p_vec);
    
    dv_mag = 6.5;
    dv = dv_mag * v_p_vec / norm(v_p_vec);
    v_p_vec = v_p_vec - dv;
    oe = sv_2_oe(r_p_vec, v_p_vec, mu);
    
    % fprintf('\nClassical Orbital Elements:\n');
    % fprintf('----------------------------\n');
    % fprintf('Specific angular momentum (h) : %.6f km^2/s\n', oe(1));
    % fprintf('Eccentricity (e)              : %.6f\n', oe(2));
    % fprintf('RA of ascending node (RA)     : %.6f deg\n', oe(3));
    % fprintf('Inclination (i)               : %.6f deg\n', oe(4));
    % fprintf('Argument of periapsis (w)     : %.6f deg\n', oe(5));
    % fprintf('True anomaly (TA)             : %.6f deg\n', oe(6));
    % fprintf('Semi-major axis (a)           : %.6f km\n', oe(7));

    oe_matrix(k, :) = oe;
end


figure;

theta_array = rad2deg(theta_array);

subplot(3,1,1)
plot(theta_array, oe_matrix(:,4)) % Inclination
xlabel('B-plane rotation angle [deg]')
ylabel('Inclination [deg]')
title('Inclination vs B-plane Angle')
grid on

subplot(3,1,2)
plot(theta_array, oe_matrix(:,5)) % Eccentricity
xlabel('B-plane rotation angle [deg]')
ylabel('RAAN [deg]')
title('Argument of periapsis vs B-plane Angle')
grid on

subplot(3,1,3)
plot(theta_array, oe_matrix(:,3)) % RAAN
xlabel('B-plane rotation angle [deg]')
ylabel('RAAN [deg]')
title('RAAN vs B-plane Angle')
grid on

function result = bplane(r_sc_neptune, v_sc_neptune, mu_neptune, thetaB_deg)
    % Inputs:
    % r_sc_jup     - Position vector of spacecraft relative to Jupiter [km]
    % v_sc_jup     - Velocity vector of spacecraft relative to Jupiter [km/s]
    % mu_jup       - Jupiter gravitational parameter [km^3/s^2]
    % thetaB_deg   - Angle of B-vector in B-plane (from T̂ toward R̂) [degrees]

    % Compute v_inf vector (incoming relative velocity)
    v_inf = v_sc_neptune;
    v_inf_mag = 17.14;

    % Define Ŝ = v_inf direction
    S_hat = v_inf / v_inf_mag;

    % Define R̂ = direction of h = r x v_inf
    h_vec = cross(r_sc_neptune, v_sc_neptune);
    R_hat = h_vec / norm(h_vec);

    % Define T̂ = R̂ × Ŝ
    T_hat = cross(R_hat, S_hat);

    % Compute orbit parameters
    h = norm(h_vec);
    a_hyp = -mu_neptune / v_inf_mag^2;
    e = 1 + h^2 / (a_hyp * mu_neptune);
    rp = a_hyp * (e - 1);  % Periapsis radius [km]
    delta = 2 * asin(1 / e);  % Turn angle [rad]
    B_mag = rp * sqrt(e^2 - 1);  % Impact parameter [km]

    % B-vector in B-plane using thetaB
    thetaB = deg2rad(thetaB_deg);
    B_vector = B_mag * (cos(thetaB) * T_hat + sin(thetaB) * R_hat);

    % Output structure
    result.S_hat = S_hat;
    result.R_hat = R_hat;
    result.T_hat = T_hat;
    result.B_vec = B_vector;
    result.B_mag = B_mag;
    result.turn_angle_deg = rad2deg(delta);
    result.rp_km = rp;
end