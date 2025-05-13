function plotInterplanetaryTransfer(mu, aunit, ip1, ip2, jdate1, jdate2, ri, vito)
    % Plot Lambert transfer and planetary orbits
    deltat = 10; % days

    % Compute periods
    [r1, v1] = planet(ip1, jdate1);
    oe1 = eci2orb1(mu, r1, v1);
    period1 = 2 * pi * oe1(1) * sqrt(oe1(1) / mu) / 86400;

    [r2, v2] = planet(ip2, jdate1);
    oe2 = eci2orb1(mu, r2, v2);
    period2 = 2 * pi * oe2(1) * sqrt(oe2(1) / mu) / 86400;

    npts1 = floor(period1 / deltat);
    npts2 = floor(period2 / deltat);
    npts3 = floor((jdate2 - jdate1) / deltat);

    x1 = zeros(npts1 + 1, 1); y1 = zeros(npts1 + 1, 1);
    x2 = zeros(npts2 + 1, 1); y2 = zeros(npts2 + 1, 1);
    x3 = zeros(npts3 + 1, 1); y3 = zeros(npts3 + 1, 1);

    for i = 0:npts1
        [r, ~] = planet(ip1, jdate1 + i * deltat);
        x1(i+1) = r(1) / aunit;
        y1(i+1) = r(2) / aunit;
    end

    for i = 0:npts2
        [r, ~] = planet(ip2, jdate1 + i * deltat);
        x2(i+1) = r(1) / aunit;
        y2(i+1) = r(2) / aunit;
    end

    for i = 0:npts3
        tau = i * deltat * 86400;
        [rft, ~] = twobody2(mu, tau, ri, vito);
        x3(i+1) = rft(1) / aunit;
        y3(i+1) = rft(2) / aunit;
    end

    figure;
    hold on;
    plot(x1, y1, '-b', 'DisplayName', 'Departure Orbit');
    plot(x2, y2, '-g', 'DisplayName', 'Arrival Orbit');
    plot(x3, y3, '-r', 'DisplayName', 'Transfer Trajectory');
    plot(0, 0, 'yo', 'MarkerSize', 10, 'DisplayName', 'Sun');
    legend;
    axis equal;
    grid on;
    xlabel('X [AU]');
    ylabel('Y [AU]');
    title('Interplanetary Lambert Transfer');
end