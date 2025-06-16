clc
clear 

days_after = 0:0.25:40; % days after launch window

data = load("earth_orbit.mat");
r = data.r_perturbed_big;
t = data.t_perturbed_big;
days_in_orbit = ceil(t(end) / 86400); % days in parking orbit
t_start = juliandate(2032,3,23 - days_in_orbit);
mu_moon = 0.00490e6; 
dv_list = zeros(size(days_after));

for day = 1:length(days_after)
    disp(day)
    r_moon = [];
    delay = days_after(day);
    t_query = t_start + t / 86400 + delay;
    [r_moon(:, 1:3), ~] = planetEphemeris(t_query', "Earth", "Moon", "430");
    r_moon = r_moon';
    
    dr = r - r_moon; % vector between s/c and moon
    norm_r = vecnorm(dr);
    a = -mu_moon * dr ./ norm_r.^3; % perturbing acceleration
    norm_a = vecnorm(a);
    dv = trapz(t, norm_a) * 1000; % convert to m/s
    dv_list(day) = dv;
    
    % plot3(r(1, :), r(2, :), r(3, :), color="black", LineWidth=1)
    % hold on
    % plot3(r_moon(1, :), r_moon(2, :), r_moon(3, :), color="red", LineWidth=1)
    % axis equal
    % legend(["Spacecraft", "Moon"])
    % xlabel("X (km)")
    % ylabel("Y (km)")
    % zlabel("Z (km)")
end

plot(days_after, dv_list, color="black", LineWidth=2)
grid on
xlabel("Launch delay (days)")
ylabel("N-body perturbing Δv (m/s)")
set(gca, "FontSize", 16)
