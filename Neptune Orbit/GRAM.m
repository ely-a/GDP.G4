clear 
clc
close all

set(groot,'defaultLineLineWidth',2) % Set all line widths to 2 unless stated otherwise
set(groot,'defaultAxesFontSize',16) % Set all axes font size to 16 unless stated otherwise
set(groot,'defaulttextfontsize',16) % Set all text font sizes to 16 unless stated otherwise
set(groot,'defaultLineMarkerSize',12) % Set all marker sizes to 16 unless stated otherwise
set(groot,'defaultAxesXGrid','on') % Turn x grid lines on
set(groot,'defaultAxesYGrid','on') % Turn y grid lines on

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

old_rhos = readmatrix("densitymodel.txt");
old_alt = old_rhos(:, 2)/1000;
old_rho = old_rhos(:,1);

% Optional: plot
figure
subplot(1, 3, 1)
semilogx(rho_mean, alt_steps, 'b'); hold on;
semilogx(rho_mean + 3*rho_std, alt_steps, 'r--');
semilogx(rho_mean - 3*rho_std, alt_steps, 'r--');
% semilogx(old_rho, old_alt, 'k-.');
xlabel('Density [kg/m^3]'); ylabel('Altitude [km]');
legend('Mean', 'Mean + 3σ', 'Mean - 3σ', 'Location', 'Best');
% xticks([1e-12 1e-10 1e-8 1e-15 1e-10 1e-5 1e0])
xticks([1e-15 1e-10 1e-5 1e0])

subplot(1,3,2)
plot(ew_mean, alt_steps); hold on;
plot(ew_mean + 3*ew_std, alt_steps, 'r--');
plot(ew_mean - 3*ew_std, alt_steps, 'r--');
xlabel('Eastward Wind [m/s]');
xticks(-500 : 250 : 500)


subplot(1,3,3)
plot(ns_mean, alt_steps); hold on;
plot(ns_mean + 3*ns_std, alt_steps, 'r--');
plot(ns_mean - 3*ns_std, alt_steps, 'r--');
xlabel('Northward Wind [m/s]');


% figure
% plot(temps_mean, alt_steps); hold on;
% plot(temps_mean + 3*temps_std, alt_steps, 'r--');
% plot(temps_mean - 3*temps_std, alt_steps, 'r--');
% xlabel('Temperature [K]');
% title('Temperature ± 3σ'); grid on;

save('gram_profiles.mat', 'alt_steps', ...
     'rho_mean', 'rho_std', ...
     'ew_mean', 'ew_std', ...
     'ns_mean', 'ns_std', ...
     'temps_mean', 'temps_std' ...
     );