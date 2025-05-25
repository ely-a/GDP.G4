function [r, v] = propogate_lagrangian(r0, v0, mu, delta_t)
    
    [~, e, h, ~, ~, ~, theta] = find_OE(r0, v0, mu);
    delta_theta = kepler_eqn(r0, v0, delta_t, mu);

    p = h^2/mu;
    theta2 = theta + delta_theta;

    r0_mag = vecnorm(r0);
    r_mag = p/(1 + e * cosd(theta2));
    
    f = 1 - r_mag / p * (1 - cosd(delta_theta));
    g = r_mag * r0_mag * sind(delta_theta) / h;

    r = f * r0 + g * v0;
    
    g_dot = 1 - r0_mag / p * (1 - cosd(delta_theta));
    f_dot = (f * g_dot - 1) / g;
    
    v = f_dot * r0 + g_dot * v0;
end