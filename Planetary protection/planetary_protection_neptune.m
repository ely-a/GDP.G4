names = {"Naiad", "Thalassa", "Despina", "Galatea", "Larissa", "Proteus"};
sma = [48224, 50075, 52526, 61953, 73548, 117646]; % km
diameter = [66, 80, 150, 174, 194, 420]; % km

% calculations
moon_area = pi * diameter .^ 2 / 4; % assume sphere as we dont have rotation
orbit_area = 4 * pi * sma .^ 2; % area swept out by orbit
probability = moon_area ./ orbit_area * 100; % in percent

% output
for i = 1:length(names)
    fprintf('Probability of impacting %s: %.6f%%\n', names{i}, probability(i));
end