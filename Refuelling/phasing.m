clear all
clc

%% delta v calculations based on masses + launch capacity

% injection
m_injection = 1500; % kg
dv_injection = 6390; % m/s
isp_injection = 310; % s (MMH/NTO)
mass_frac_cruise = 0.85; % proportion of cruise stage that is fuel
m_cruise = m_injection * (1 + exp(dv_injection ...
    / (isp_injection * 9.81)) / mass_frac_cruise);

% earth ejection (direct)
dv_ejection = 7780; % m/s
isp_ejection = 397; % s (LOx/RP-1)
mass_frac_kick = 0.9; % proportion of kick stage that is fuel
m_direct = m_cruise * (1 + exp(dv_ejection ...
    / (isp_ejection * 9.81)) / mass_frac_kick);

% earth ejection (with refuelling)
FH1_list = [];
FH2_list = [];
apo_list = [];
mu = 398600; % km^3/s^2
Re = 6378; % km
rp = Re + 300; % km
v_LEO = sqrt(mu / rp); % km/s
v_esc = v_LEO * sqrt(2); % km/s

% vehicle parameters
fuel_distribution = 0; % proportion of fuel carried in FH1
mass_frac_kick = 0.9; % proportion of kick stage that is fuel
mass_frac_tanker = 1; % proportion of tanker that is fuel
isp_FHUS = 397; % s, FH upper stage isp

for dv_1 = linspace(0, (v_esc - v_LEO) * 1000, 100) % LEO to escape

    % mass calculations
    dv_2 = dv_ejection - dv_1; % m/s, performed after docking
    m_docked = m_cruise * (1 + exp(dv_2 ...
        / (isp_ejection * 9.81)) / mass_frac_kick); % kg
    m_FH1 = m_cruise + (m_docked - m_cruise) * (1 - mass_frac_kick); % kg
    m_FH2 = (m_docked - m_FH1) / mass_frac_tanker; % kg

    % move fuel around
    m_FH1 = m_FH1 + (m_docked - m_FH1) * fuel_distribution;
    m_FH2 = m_FH2 - (m_docked - m_FH1) * fuel_distribution;

    % transfer to docking orbit
    m_LEO_FH1 = m_FH1 * (1 + exp(dv_1 ...
        / (isp_FHUS * 9.81))); % kg
    m_LEO_FH2 = m_FH2 * (1 + exp(dv_1 ...
        / (isp_FHUS * 9.81))); % kg
    FH1_list = [FH1_list, m_LEO_FH1 / 1000];
    FH2_list = [FH2_list, m_LEO_FH2 / 1000];

    % orbit calculations
    vp = v_LEO + dv_1 / 1000; % km/s
    h = vp * rp;
    e = h^2 / (rp * mu) - 1;
    ra = h^2 / mu * (1 / (1 - e));
    apo_list = [apo_list, ra];
end

% plotting
figure;
semilogx(apo_list, FH1_list, LineWidth=2, color="red");
hold on
semilogx(apo_list, FH2_list, LineWidth=2, color="blue");
yline(63.8, LineStyle="--")
yline(63.8 * 0.9, LineStyle="--", color="red")
grid on
xlabel("Apoapsis (km)")
ylabel("Mass to LEO (tonnes)")
legend(["Cruise Vehicle", "Fuel in tanker", "Falcon Heavy Capacity", "10% Margin"])
ylim([0, max(max(FH1_list, FH2_list))])

%% phasing calculations

% new apoapsis
rp = 300 + Re; % km
ra = 48788; % km
sma = 0.5 * (ra + rp); % km
orbit_num = 4; % number of phasing orbits
e = (ra - rp) / (ra + rp);
T_orbit = 2 * pi * sqrt(sma^3 / mu); % s
dt_phasing = 3600 / orbit_num; % number of seconds we are off by
T_new = T_orbit + dt_phasing; % s
sma_new = ((T_new / (2 * pi)) ^ 2 * mu) ^ (1/3); % km
ra_new = 2 * sma_new - rp; % km
e_new = (ra_new - rp) / (ra_new + rp);

% delta v
h = sqrt(rp * mu * (1 + e)); % km^2/s
h_new = sqrt(rp * mu * (1 + e_new)); % km^2/s
vp = h / rp; % km/s
vp_new = h_new / rp; % km/s
dv_total = 2 * abs(vp - vp_new); % km/s

%% finite burn losses

% for r_p = 42078 km
% dv_1 = 2.4244 km/s before docking
% dv_2 = 5.2756 km/s after docking

% initial orbit
e = 0;
ra = 48878; % km
rp = 6678; % km'
sma = rp;
dv = 2.4244; % km/s
acc = 10; % burn acceleration, m/s^2
t = dv * 1000 / acc;
grav_loss_1 = 1/24 * (mu * (1 + e)) / (sma^3 * (1 - e)^3) * t^2 * dv; % km/s

% second orbit
e = (ra - rp) / (ra + rp);
sma = (6678 + 42078) / 2; % km
dv = 5.2756; % km/s
acc = 10; % burn acceleration, m/s^2
t = dv * 1000 / acc;
grav_loss_2 = 1/24 * (mu * (1 + e)) / (sma^3 * (1 - e)^3) * t^2 * dv; % km/s