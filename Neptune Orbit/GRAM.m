clear 
clc
close all

winds = readmatrix("Winds.txt");
density = readmatrix("Density.txt");
temperatures = readmatrix("Temp.txt");

alts = winds(:, 1);
winds_EW = winds(:, 4);
winds_NS = winds(:, 7);
rhos = density(:, 5);
temps = temperatures(:, 2);

alt_steps = 0:10:4000;

MCs = length(alts)/length(alt_steps);

% Determine number of altitude steps
num_altitudes = length(alt_steps);

% Reshape data: rows = altitudes, cols = MC runs
rhos_mat      = reshape(rhos,      num_altitudes, MCs);
winds_EW_mat  = reshape(winds_EW,  num_altitudes, MCs);
winds_NS_mat  = reshape(winds_NS,  num_altitudes, MCs);
temps_mat = reshape(temps, num_altitudes, MCs);

% Compute mean and standard deviation for each row (altitude)
rho_mean   = mean(rhos_mat, 2);
rho_std    = std(rhos_mat,  0, 2);

ew_mean    = mean(winds_EW_mat, 2);
ew_std     = std(winds_EW_mat,  0, 2);

ns_mean    = mean(winds_NS_mat, 2);
ns_std     = std(winds_NS_mat,  0, 2);

temps_mean = mean(temps_mat, 2);
temps_std     = std(temps_mat,  0, 2);

% old_rhos = readmatrix("densitymodel.txt");
% old_alt = old_rhos(:, 2)/1000;
% old_rho = old_rhos(:,1);

% Optional: plot
figure
semilogx(rho_mean, alt_steps); hold on;
semilogx(rho_mean + 3*rho_std, alt_steps, 'r--');
semilogx(rho_mean - 3*rho_std, alt_steps, 'r--');
% semilogx(old_rho, old_alt, 'b--');
xlabel('Density [kg/m^3]'); ylabel('Altitude [km]');
title('Density ± 3σ (log scale)'); grid on;

figure
subplot(1,2,1)
plot(ew_mean, alt_steps); hold on;
plot(ew_mean + 3*ew_std, alt_steps, 'r--');
plot(ew_mean - 3*ew_std, alt_steps, 'r--');
xlabel('Eastward Wind [m/s]');
title('Eastward Wind ± 3σ'); grid on;

subplot(1,2,2)
plot(ns_mean, alt_steps); hold on;
plot(ns_mean + 3*ns_std, alt_steps, 'r--');
plot(ns_mean - 3*ns_std, alt_steps, 'r--');
xlabel('Northward Wind [m/s]');
title('Northward Wind ± 3σ'); grid on;

figure
plot(temps_mean, alt_steps); hold on;
plot(temps_mean + 3*temps_std, alt_steps, 'r--');
plot(temps_mean - 3*temps_std, alt_steps, 'r--');
xlabel('Temperature [K]');
title('Temperature ± 3σ'); grid on;

save('gram_profiles.mat', 'alt_steps', ...
     'rho_mean', 'rho_std', ...
     'ew_mean', 'ew_std', ...
     'ns_mean', 'ns_std', ...
     'temps_mean', 'temps_std' ...
     );