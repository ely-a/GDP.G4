clear all
close all
clc

% inputs
mu = 6.836529e6;
RN = 24622; % neptune radius, km
rp = 10000 + RN; % periapsis radius, km
ra_array = (1e4:500:1.5e6) + RN; % apoapsis radius, km

% instrument ranges
instrument_ranges = {
    'CCD/CMOS Imager',      [500 50000];
    'Dual Fluxgate/Scalar Magnetometer',  [37000 738000];
    'IR Spectrometer',            [500 30000];
    'VIS/UV Spectrometer',        [500 50000];
    'Laser Altimeter',    [500 2000];
    'Microwave Radiometer',     [1000 6000];
    'Bolometer',     [0 2000];
    'Plasma Spectrometer',        [25 1e6];
    'Radar Array',         [1000 4000];
};

time_list = [];
j = 1;
n_inst = size(instrument_ranges,1);
n_apo = length(ra_array);
proportions = zeros(n_inst, n_apo);

for apo = ra_array
    a = (rp + apo)/2;
    e = (apo - rp)/(apo + rp);
    T = 2*pi*sqrt(a^3/mu);  % orbital period

    fprintf('\n--- Apoapsis: %.0f km ---\n', apo);
    
    for i = 1:size(instrument_ranges,1)
        name = instrument_ranges{i,1};
        range = instrument_ranges{i,2} + RN;  % convert to radius from Neptune center
        rmin = max(range(1), rp);
        rmax = min(range(2), apo);

        if rmin >= rmax
            proportion = 0;
        else
            % Orbit equation: r = a(1 - e^2) / (1 + e*cos(theta))
            cos_th1 = (a*(1 - e^2)/rmax - 1)/e;
            cos_th2 = (a*(1 - e^2)/rmin - 1)/e;

            if abs(cos_th1) > 1.01 || abs(cos_th2) > 1.01
                proportion = 0;
            else
                th1 = acos(cos_th1);
                th2 = acos(cos_th2);

                % Eccentric anomalies
                E1 = 2*atan( sqrt((1 - e)/(1 + e)) * tan(th1/2) );
                E2 = 2*atan( sqrt((1 - e)/(1 + e)) * tan(th2/2) );
                if E1 < 0, E1 = E1 + 2*pi; end
                if E2 < 0, E2 = E2 + 2*pi; end

                t1 = (E1 - e*sin(E1)) * sqrt(a^3/mu);
                t2 = (E2 - e*sin(E2)) * sqrt(a^3/mu);

                % Time in range is symmetric around periapsis
                time_in_range = 2*(t1 - t2);
                proportion = time_in_range / T;
            end
        end

        fprintf('%s: %.2f%% of orbit in range\n', name, proportion*100);
        time_list(i, j) = proportion*100;
    end
    j = j + 1;
end

% Plotting
figure;
hold on;
colors = lines(n_inst);  % distinct colors
lineStyles = {'-', '--', ':', '-.', '-', '--', ':', '-.', '-'};  % at least 9 styles

for i = 1:n_inst
    style = lineStyles{mod(i-1, length(lineStyles)) + 1};  % cycle through line styles
    plot(ra_array - RN, time_list(i,:), ...
         'LineWidth', 2, ...
         'DisplayName', instrument_ranges{i,1}, ...
         'Color', colors(i,:), ...
         'LineStyle', style);
end

xlabel('Apoapsis Altitude (km)');
ylabel('Time in Range (% of Orbit)');
title('Instrument Access Time vs. Apoapsis');
legend('Location', 'eastoutside');
grid on;
xlim([min(ra_array - RN), max(ra_array - RN)]);
ylim([0 100]);

% calculations
e = (ra_array - rp) ./ (ra_array + rp);
h = sqrt(rp*mu*(1+e));
vp = h / rp;
dv = 26.234 - vp;
a = (rp + ra_array)/2;
T = 2*pi*sqrt(a.^3/mu);  % orbital period

% dv plot
figure
plot(ra_array - RN, dv, color="black", LineWidth=2);
xlabel('Apoapsis Altitude (km)');
ylabel('ΔV from cruise trajectory (km/s)');
title('Capture ΔV vs. Apoapsis');
grid on;
xlim([min(ra_array - RN), max(ra_array - RN)]);

% period plot
figure
plot(ra_array - RN, T/3600, color="black", LineWidth=2);
xlabel('Apoapsis Altitude (km)');
ylabel('Orbital period (hours)');
title('Orbital period vs. Apoapsis');
grid on;
xlim([min(ra_array - RN), max(ra_array - RN)]);