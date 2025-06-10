clc 
clear
close all

% Constants
mu_S = 1.327e11;

% Load data
trajectories = load("trajectory.mat");
velocities = load("velocity.mat");


% Postion and velocity of Jupiter
r_neptune = trajectories.r(:,19:21);
r_N = r_neptune(end,:);
v_N = [-2.99130437	4.21070026	1.79770374];


% 467 days to get to jupiter
r_helio = trajectories.r_out(467:end,:);
v_helio = velocities.v_out(467:end,:);


% Position and velocity wrt Jupiter
r_sc_j = r_helio - r_N;
v_sc_j = v_helio - v_N;


% Define B-Plane
v_inf_vec = v_sc_j(end-2,:);                     % Incoming hyp. velocity vector
S_vec = v_inf_vec/norm(v_inf_vec);        % Unit vector parallel to incoming hyp. velocity vector (normal to b_plane)
h_vec = cross(r_sc_j(end-2,:), v_sc_j(end-2,:));        % Specific angular momentum vector
R_vec = h_vec / norm(h_vec);              % R vector
T_vec = cross(R_vec, S_vec);              % T vector


% Find point of intersection with the B-Plane
dist = 10000 + 24764;        
phi = -181.3;            % B-plane rotation angle in º (-ve for clockwise rotation)
% Rodrigues' formula
B_vec_unit = T_vec * cosd(phi) + cross(S_vec, T_vec) * sind(phi) + S_vec * dot(S_vec, T_vec) * (1 - cosd(phi));
B_vec = dist*B_vec_unit;
tp = B_vec;                 % Target point wrt to Jupiter
tp_helio = tp + r_N;        % Target point wrt to the Sun


% Lambert solver
days = length(r_helio);         % Days to get to Neptune from Jupiter
days_sec = days*24*3600;        % Days to get to Neptune from Jupiter in seconds
delta_V = zeros(1,days);        % Initialise delta V vector
for i = 1:days
    [v1_helio, v2_helio] = lambert(r_helio(i,:), r_N, days_sec - i*(24*3600), "pro",mu_S);
    [v1, v2] = lambert(r_helio(i,:), tp_helio, days_sec - i*(24*3600), "pro",mu_S);
    dV_vec = v1 - v1_helio;
    dV_mag = norm(dV_vec);
    delta_V(i) = dV_mag*1000;
end


figure
semilogy(467:days+466-1, delta_V(1:end-1),"r",LineWidth=1.5)
grid on
xlabel("Days from Earth departure");
ylabel("\DeltaV (m/s)");
xlim([0 3050]);
ylim([1e-1 2e4]);



% from the graph: minimum at 467, dv=0.185433230751888 m/s