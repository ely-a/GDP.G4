clc 
clear
close all

% Constants
mu_S = 1.327e11;


% 467 days to get to jupiter
trajectories = load("trajectory.mat");
velocities = load("velocity.mat");
r = trajectories.r_out(1:467,:);
v = velocities.v_out(1:467,:);

% Postion and velocity wrt Jupiter
r_jupiter = trajectories.r(:,10:12);
r_J = r_jupiter(467,:);
v_J = [6.97159704,10.70124642,4.41714965];


% Define B-Plane
v_inf_vec = v(465,:);                     % Incoming hyp. velocity vector
S_vec = v_inf_vec/norm(v_inf_vec);        % Unit vector parallel to incoming hyp. velocity vector (normal to b_plane)
h_vec = cross(r(465,:), v(465,:));        % Specific angular momentum vector
R_vec = h_vec / norm(h_vec);              % R vector
T_vec = cross(R_vec, S_vec);              % T vector


% Find point of intersection with the B-Plane
dist = 27*71492;        % 27 radii
phi = -3.05;            % B-plane rotation angle in º (-ve for clockwise rotation)
% Rodrigues' formula
B_vec_unit = T_vec * cosd(phi) + cross(S_vec, T_vec) * sind(phi) + S_vec * dot(S_vec, T_vec) * (1 - cosd(phi));
B_vec = dist*B_vec_unit;
tp = B_vec;  % Target point wrt to Neptune


% Lambert solver
days = length(r);               % Days to get to Jupiter from Earth
days_sec = days*24*3600;        % Days to get to Jupiter from Earth in seconds
delta_V = zeros(1,length(r));   % Initialise delta V vector
for i = 1:length(r)
    %[v1, v2] = lambert(r(i,:), r(days,:), days_sec - i*(24*3600), "pro",mu_S);
    [v1, v2] = lambert(r(i,:), B_vec, days_sec - i*(24*3600), "pro",mu_S);
    dV_vec = v2-v1;
    dV_mag = norm(dV_vec);
    delta_V(i) = dV_mag;
end


figure
plot(1:days, delta_V)
grid on
xlabel("Days from launch");
ylabel("\DeltaV (km/s)");






