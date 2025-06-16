function delta_theta = kepler_eqn(r, v, delta_t, mu)

    [a, e, ~, ~, ~, ~, theta] = oe_from_rv(r, v, mu);
    n = sqrt(mu/a^3);
    
    E0 = 2 * atan(sqrt((1 - e)/(1 + e)) * tand(theta/2)); % radians
    M0 = E0 - e * sin(E0);

    M1 = M0 + n * delta_t;
    E1 = solve_kepler(M1, e, 1e-8, 1000);
    
    theta2 = 2 * atand(sqrt((1 + e)/(1 - e)) * tan(E1/2));
    delta_theta = theta2 - theta;
end