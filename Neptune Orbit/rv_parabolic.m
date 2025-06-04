function [r, v] = rv_parabolic(r_p, v_p, i, Omega, omega)
    % Input:
    % r_p     : periapsis radius [km]
    % i       : inclination [deg]
    % Omega   : RAAN [deg]
    % omega   : argument of periapsis [deg]
    % mu      : gravitational parameter [km^3/s^2]

    % Position at periapsis in perifocal frame
    r_perif = [r_p; 0; 0];
    
    % Velocity at periapsis in perifocal frame (for parabolic orbit)
    v_perif = [0; v_p; 0];

    % Rotation matrices
    R3_Omega = [cosd(Omega), -sind(Omega), 0;
                sind(Omega),  cosd(Omega), 0;
                0,            0,           1];

    R1_i = [1, 0, 0;
            0, cosd(i), -sind(i);
            0, sind(i),  cosd(i)];

    R3_omega = [cosd(omega), -sind(omega), 0;
                sind(omega),  cosd(omega), 0;
                0,            0,           1];

    % Final rotation
    R = R3_Omega * R1_i * R3_omega;

    % Convert to inertial
    r = R * r_perif;
    v = R * v_perif;
end
