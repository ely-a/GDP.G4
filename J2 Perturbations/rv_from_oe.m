function [r, v] = rv_from_oe(a, e, i, Omega, omega, theta_deg, mu)
    theta = deg2rad(theta_deg);

    % Perifocal coordinates
    p = a * (1 - e^2);
    r_perif = (p / (1 + e * cos(theta))) * [cos(theta); sin(theta); 0];
    v_perif = sqrt(mu/p) * [-sin(theta); e + cos(theta); 0];

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

    Q = R3_Omega * R1_i * R3_omega;

    r = Q * r_perif;
    v = Q * v_perif;
end
