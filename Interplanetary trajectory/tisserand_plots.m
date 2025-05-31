radii = [0.72334, 1.00000, 1.52366, 5.20260, 9.55491, 19.2184, 30.0583] * 1.496e8; % km
min_flyby = [300, 300, 200, 1e5, 5e4, 5e4, 5e4];
planet_names = {'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};
mu_sun = 1.3271244e11; % km^3/s^2

v_inf_values = 2:12; % km/s
alpha_values = linspace(0, pi, 1e3); % fine angle resolution

colors = [0.9 0.6 0;
          0 0.6 0.8;
          0.9 0.2 0;
          0.8 0.5 0.3;
          0.6 0.6 0.3;
          0.5 0.7 0.7;
          0 0 0.8];


figure; hold on;

for planet = 1:7
    r = radii(planet);
    v_circ = sqrt(mu_sun / r);
    
    periapsis_all = [];
    energy_all = [];
    v_inf_all = [];
    
    for v_inf = v_inf_values
        for alpha = alpha_values
            v_tan = v_circ + v_inf * cos(alpha);
            v_rad = v_inf * sin(alpha);
            final_v = [v_tan, v_rad];
            
            h = r * v_tan;
            energy = 0.5 * norm(final_v)^2 - mu_sun / r;
            a = -mu_sun / (2 * energy);
            if a <= 0, continue; end % skip unbound
            
            e = sqrt(1 - (h^2 / (mu_sun * a)));
            periapsis = a * (1 - e);
            
            periapsis_all(end+1) = periapsis;
            energy_all(end+1) = energy;
            v_inf_all(end+1) = v_inf;
        end

        % max deflection code
        ...

    end

    % Plot lines for each v_inf
    for i = 1:length(v_inf_values)
        vi = v_inf_values(i);
        mask = (v_inf_all == vi);
        x = periapsis_all(mask);
        y = energy_all(mask);
        if isempty(x), continue; end

        % Sort by periapsis for clean plotting
        [x_sorted, idx] = sort(x);
        y_sorted = y(idx);

        plot(x_sorted / 1e6, y_sorted, ...
            'Color', colors(planet, :), ...
            'LineWidth', 0.5);

        if i == 1 || i == length(v_inf_values)
            % Label the line with v_inf value at a chosen point
            label_idx = round(length(x_sorted) * 0.9);  % 90% along the line
            text_x = x_sorted(label_idx) / 1e6;
            text_y = y_sorted(label_idx);
            text_str = sprintf('v_\\infty = %d km/s', vi);
            text(text_x, text_y, text_str, ...
                'Color', "black", ...
                'FontSize', 10, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom');
        end
    end
end

set(gca, 'XScale', 'log');
xlim([20 1000])
ylim([-900 0])
xlabel('Periapsis Radius [10^6 km]');
ylabel('Specific Energy [km^2/s^2]');
grid on;