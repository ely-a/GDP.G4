% Constants
mu_sun = 132712440000;        % Solar gravitational parameter (km^3/s^2)
AU_to_km = 149597870.7;       % 1 AU in km

% Planetary data [Venus to Neptune]
planets = {'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};
r_planet_AU = [0.72, 1.0, 1.52, 5.20, 9.58, 19.22, 30.05];
planet_radii = [6051, 6378, 3389, 69911, 58232, 25362, 24622]; % km
planet_mu = [324858, 398600, 42828, 126686534, 37931207, 5793964, 6835100]; % km^3/s^2
min_flyby = [300, 300, 200, 1e5, 5e4, 5e4, 5e4]; % km

% Convert planetary distances to km
r_planet_km = r_planet_AU * AU_to_km;

% Velocity parameters
V_infs = 2:2:12;             % km/s
theta = 0:0.5:180;           % Fine angular resolution

% Initialize figure
figure;
hold on;
colors = lines(length(planets));
legend_handles = gobjects(length(planets), 1);

% Main plotting loop
for i = 1:length(planets)
    % Planetary parameters
    r_p_planet = r_planet_km(i);
    V_p = sqrt(mu_sun/r_p_planet);
    r_flyby = planet_radii(i) + min_flyby(i);
    mu_planet = planet_mu(i);
    
    % Legend entry
    legend_handles(i) = plot(NaN, NaN, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Contour generation
    for j = 1:length(V_infs)
        V_inf = V_infs(j);
        [E_vals, r_p_vals] = deal([]);
        
        % Calculate maximum bending angle per flyby
        delta_max = 2*asind(1/(1 + (r_flyby*V_inf^2)/mu_planet));
        
        % Calculate tick mark spacing (max angular change per flyby)
        tick_spacing = delta_max;
        
        % Generate tick angles
        tick_angles = 0:tick_spacing:180;
        if mod(180, tick_spacing) ~= 0
            tick_angles = [tick_angles 180];
        end
        
        % Store previous valid tick position
        prev_tick = [];
        
        for k = 1:length(theta)
            th = theta(k);
            V_helio = sqrt(V_p^2 + V_inf^2 + 2*V_p*V_inf*cosd(th));
            E = V_helio^2/2 - mu_sun/r_p_planet;
            
            if E > 0, continue; end % Skip hyperbolic
            
            % Calculate orbital parameters
            h = r_p_planet*(V_p + V_inf*cosd(th));
            a = -mu_sun/(2*E);
            e = sqrt(1 + (2*E*h^2)/mu_sun^2);
            r_p_peri = a*(1-e)/1e6; % Convert to million km
            
            % Store values
            E_vals(end+1) = E;
            r_p_vals(end+1) = r_p_peri;
            
            % Add tick marks at specified intervals
            if any(abs(th - tick_angles) < 0.5)
                if ~isempty(prev_tick) && (abs(r_p_peri - prev_tick(1)) < 5 && abs(E - prev_tick(2)) < 5)
                    continue; % Prevent overlapping ticks
                end
                scatter(r_p_peri, E, 40, 'filled', ...
                       'MarkerEdgeColor', colors(i,:), ...
                       'MarkerFaceColor', colors(i,:));
                prev_tick = [r_p_peri, E];
            end
        end
        
        % Plot contour line
        if ~isempty(r_p_vals)
            plot(r_p_vals, E_vals, 'Color', colors(i,:), 'LineWidth', 1.5);
        end
    end
end

% Plot configuration
set(gca, 'XScale', 'log', 'XLim', [20 1000], 'YLim', [-900 0]);
xlabel('Periapsis Distance (Mkm)');
ylabel('Specific Energy (km²/s²)');
legend(legend_handles, planets, 'Location', 'northeastoutside');
title('Tisserand Plot with Flyby Tick Marks');
grid on;
hold off;
