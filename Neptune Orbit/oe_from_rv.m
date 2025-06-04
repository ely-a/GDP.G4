%% Find Orbital elements from r and v

function [a, e_mag, h_mag, i, Omega, omega, theta] = oe_from_rv(r, v, mu)
    r_mag = vecnorm(r);
    v_mag = vecnorm(v);
    
    h = cross(r, v);
    h_mag = vecnorm(h);
    
    a = -mu / (v_mag^2 - 2*mu/r_mag);
    
    e = 1/mu * (cross(v,h) - mu*r/r_mag);
    e_mag = vecnorm(e);
    
    i = acosd(h(3)/h_mag);
    
    N = cross([0;0;1], h);
    N_mag = vecnorm(N);
    Omega = acosd(N(1)/N_mag);
    if N(2) < 0
        Omega = 360 - Omega;
    end
    
    omega = acosd(dot(N, e)/(N_mag * e_mag));
    if e(3) < 0
        omega = 360 - omega;
    end
    
    cos_theta = dot(r, e)/(r_mag * e_mag);
    cos_theta = min(max(cos_theta, -1), 1); % Clamp to [-1, 1]
    theta = acosd(cos_theta);
    if dot(r, v) < 0
        theta = 360 - theta;
    end
end