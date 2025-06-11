%tvishi tps thickness code 
function [m_TPS, t_TPS_cm] = HeatShieldMass(material, alt_km, vel_km_s, t_vals, temps, density)
    % Universal constants
    gamma = 1.3846;
    Rg    = 3.61491e3;      % J/(kg·K)
    sigma = 5.67e-8;        % W/(m^2·K^4)
    emis  = 0.864;          % emissivity
    Rn    = 3;              % nose radius (m)
    T0    = 288;            % ambient temperature (K)

    % --- Pull in material properties ---
    [k, rhoTPS, cpTPS] = materialProperties(material);

    % --- Aerodynamic profiles ---
    Tprof = temps(alt_km);
    a     = sqrt(gamma * Rg .* Tprof);
    v     = vel_km_s * 1e3;    % convert km/s to m/s
    v_fake= vel_km_s * 1e2;
    Mach  = v ./ a;

    % --- Convective heat flux (Sutton-Graves) ---
    K2     = 1/((0.667/0.0395)+(0.333/0.0797));
    density_profile = density(alt_km) / 1e9;
    q_conv = (K2/2) * sqrt(density_profile./Rn) .* v_fake.^3;

    % --- Wall temperature (radiative equilibrium) ---
    T_wall = (q_conv ./ (sigma * emis)).^0.25;

    % --- Radiative back-flux (match original Python) ---
    q_rad  = emis * (T_wall - T0);

    % --- Net heat flux & integration ---
    q_total = q_conv - q_rad;
    Q_total = trapz(t_vals, q_total);

    % --- TPS sizing ---
    A_nose   = pi * Rn^2;
    deltaT   = max(T_wall) - T0;
    m_areal  = Q_total / (cpTPS * deltaT);  % kg/m^2
    t_TPS    = m_areal / rhoTPS;           % m
    m_TPS    = m_areal * A_nose;           % kg
    t_TPS_cm = 100 * t_TPS;

    % --- Pack results ---
    % res.k         = k;
    % res.rho       = rhoTPS;
    % res.cp        = cpTPS;
    % res.Mach      = Mach;
    % res.velocity  = vel_km_s;
    % res.T_wall    = T_wall;
    % res.q_total   = q_total;
    % res.Q_total   = Q_total;
    % res.t_TPS     = t_TPS;
    % res.t_TPS_cm  = t_TPS * 100;
    % res.m_TPS     = m_TPS;
end

function [k, rho, cp] = materialProperties(material)
    % materialProperties  Lookup TPS material constants
    switch lower(material)
        case 'tiles'
            k   = 0.88;    % W/(m·K)
            rho = 961;    % kg/m^3
            cp  = 1250;    % J/(kg·K) 
        otherwise
            error('Unsupported material: %s', material);
    end
end