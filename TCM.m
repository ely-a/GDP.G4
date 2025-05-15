% make sure r_plot is in the workspace first
% of the order [Venus_x Venus_y Venus_z Earth_x ... S/C_x S/C_y S/C_z
v_plot = diff(r_plot, 1, 1);
tcm_interval = 30; % days
