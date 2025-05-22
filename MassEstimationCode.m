% Housekeeping
clear
clc
close all

% Enter the rocket boost capabilities (kg)
launchers = struct('Ariane6', 21650, 'FalconHeavy', 50000, 'Starship', 100000, 'Atlas5', 18850);

% Terminology: GA = gravity assist, dir = direct transfer
% init = initial trajectory burn, inj = Neptune injection

% =========================================================================
% CHANGE ΔV HERE (in km/s)
dV_init_GA = 9.5;
dV_init_dir = 6.5;
dV_inj_GA = 9.5;
dV_inj_dir = 8;
% =========================================================================

% =========================================================================
% CHANGE STRUCTURAL MASS FRACTION HERE
epsilon_s = 0.15;
% =========================================================================

% Specific impulses (s)
Isp_H2O2 = 450; % for initial burn
Isp_N2O4 = 290; % for injection burn (chemical)
Isp_EP = 1350; % for injection burn (electric)

g0 = 9.81; % m/s^2

% Pre-allocate structures for results
propellant = struct();

launcherNames = fieldnames(launchers);

for i = 1:length(launcherNames)
    launcher = launcherNames{i};
    m0 = launchers.(launcher);

    % Structural mass 
    m_struct = epsilon_s * m0;

    % Usable mass for propellant + payload
    m_avail = m0 - m_struct;

    % Convert ΔV to m/s
    dV_init_GA_m = dV_init_GA * 1000;
    dV_init_dir_m = dV_init_dir * 1000;
    dV_inj_GA_m = dV_inj_GA * 1000;
    dV_inj_dir_m = dV_inj_dir * 1000;

    % Initial burn propellant mass (always H2O2)
    propellant.(launcher).init_GA = m_avail * (1 - exp(-dV_init_GA_m / (Isp_H2O2 * g0)));
    propellant.(launcher).init_dir = m_avail * (1 - exp(-dV_init_dir_m / (Isp_H2O2 * g0)));

    % Injection burn propellant mass with N2O4
    propellant.(launcher).inj_GA_N2O4 = (m_avail - propellant.(launcher).init_GA) * (1 - exp(-dV_inj_GA_m / (Isp_N2O4 * g0)));
    propellant.(launcher).inj_dir_N2O4 = (m_avail - propellant.(launcher).init_dir) * (1 - exp(-dV_inj_dir_m / (Isp_N2O4 * g0)));

    % Injection burn propellant mass with Electric Propulsion
    propellant.(launcher).inj_GA_EP = (m_avail - propellant.(launcher).init_GA) * (1 - exp(-dV_inj_GA_m / (Isp_EP * g0)));
    propellant.(launcher).inj_dir_EP = (m_avail - propellant.(launcher).init_dir) * (1 - exp(-dV_inj_dir_m / (Isp_EP * g0)));

    % Final payload masses delivered into Neptunian orbit (excluding structure and propellant)
    propellant.(launcher).nep_orb_GA_N2O4 = m_avail - propellant.(launcher).init_GA - propellant.(launcher).inj_GA_N2O4;
    propellant.(launcher).nep_orb_GA_EP = m_avail - propellant.(launcher).init_GA - propellant.(launcher).inj_GA_EP;
    propellant.(launcher).nep_orb_dir_N2O4 = m_avail - propellant.(launcher).init_dir - propellant.(launcher).inj_dir_N2O4;
    propellant.(launcher).nep_orb_dir_EP = m_avail - propellant.(launcher).init_dir - propellant.(launcher).inj_dir_EP;

    % Store structural mass for reporting
    propellant.(launcher).struct_mass = m_struct;
end

% Write mission results to a text file
fileID = fopen('mission_results.txt', 'w');

fprintf(fileID, '\nMission Results to Neptunian Orbit (Gravity Assist + N₂O₄):\n');
fprintf(fileID, '------------------------------------------------------------\n');
for i = 1:length(launcherNames)
    launcher = launcherNames{i};
    total_prop_mass = propellant.(launcher).init_GA + propellant.(launcher).inj_GA_N2O4;
    delivered_payload = propellant.(launcher).nep_orb_GA_N2O4;

    fprintf(fileID, '%s\n', launcher);
    fprintf(fileID, '  Structural Mass (kg): %.2f\n', propellant.(launcher).struct_mass);
    fprintf(fileID, '  Propellant Mass (kg): %.2f\n', total_prop_mass);
    fprintf(fileID, '  Delivered Payload (kg): %.2f\n\n', delivered_payload);
end

fprintf(fileID, '\nMission Results to Neptunian Orbit (Gravity Assist + EP):\n');
fprintf(fileID, '------------------------------------------------------------\n');
for i = 1:length(launcherNames)
    launcher = launcherNames{i};
    total_prop_mass = propellant.(launcher).init_GA + propellant.(launcher).inj_GA_EP;
    delivered_payload = propellant.(launcher).nep_orb_GA_EP;

    fprintf(fileID, '%s\n', launcher);
    fprintf(fileID, '  Structural Mass (kg): %.2f\n', propellant.(launcher).struct_mass);
    fprintf(fileID, '  Propellant Mass (kg): %.2f\n', total_prop_mass);
    fprintf(fileID, '  Delivered Payload (kg): %.2f\n\n', delivered_payload);
end

fprintf(fileID, '\nMission Results to Neptunian Orbit (Direct + N₂O₄):\n');
fprintf(fileID, '------------------------------------------------------------\n');
for i = 1:length(launcherNames)
    launcher = launcherNames{i};
    total_prop_mass = propellant.(launcher).init_dir + propellant.(launcher).inj_dir_N2O4;
    delivered_payload = propellant.(launcher).nep_orb_dir_N2O4;

    fprintf(fileID, '%s\n', launcher);
    fprintf(fileID, '  Structural Mass (kg): %.2f\n', propellant.(launcher).struct_mass);
    fprintf(fileID, '  Propellant Mass (kg): %.2f\n', total_prop_mass);
    fprintf(fileID, '  Delivered Payload (kg): %.2f\n\n', delivered_payload);
end

fprintf(fileID, '\nMission Results to Neptunian Orbit (Direct + EP):\n');
fprintf(fileID, '------------------------------------------------------------\n');
for i = 1:length(launcherNames)
    launcher = launcherNames{i};
    total_prop_mass = propellant.(launcher).init_dir + propellant.(launcher).inj_dir_EP;
    delivered_payload = propellant.(launcher).nep_orb_dir_EP;

    fprintf(fileID, '%s\n', launcher);
    fprintf(fileID, '  Structural Mass (kg): %.2f\n', propellant.(launcher).struct_mass);
    fprintf(fileID, '  Propellant Mass (kg): %.2f\n', total_prop_mass);
    fprintf(fileID, '  Delivered Payload (kg): %.2f\n\n', delivered_payload);
end

fclose(fileID);