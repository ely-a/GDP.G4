% lambert2.m     January 26, 2024

% solution of the interplanetary Lambert problem

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global mu

% define Astronomical Unit (kilometers)

aunit = 149597870.691;

% load planet name array

pdata = ['Mercury'; 'Venus  '; 'Earth  '; 'Mars   '; ...
    'Jupiter'; 'Saturn '; 'Uranus '; 'Neptune'; 'Pluto  '];

pname = cellstr(pdata);

% begin simulation

clc; home;

fprintf('\n            program lambert2\n');

fprintf('\n   < interplanetary lambert problem >\n');

fprintf('\ndeparture conditions\n');

% request departure calendar date

[month, day, year] = getdate;

% request departure ut

[uthr, utmin, utsec] = gettime;

jdate1 = julian(month, day, year);

jdate1 = jdate1 + (uthr + utmin / 60 + utsec / 3600) / 24;

fprintf('\narrival conditions \n');

% request arrival calendar date

[month, day, year] = getdate;

jdate2 = julian(month, day, year);

% request arrival ut

[uthr, utmin, utsec] = gettime;

jdate2 = jdate2 + (uthr + utmin / 60 + utsec / 3600) / 24;

% transfer time

taud = jdate2 - jdate1;

tof = taud * 86400;

% request departure and arrival planets

for i = 1:1:2

    fprintf('\n planetary menu\n');

    fprintf('\n  <1> Mercury');
    fprintf('\n  <2> Venus');
    fprintf('\n  <3> Earth');
    fprintf('\n  <4> Mars');
    fprintf('\n  <5> Jupiter');
    fprintf('\n  <6> Saturn');
    fprintf('\n  <7> Uranus');
    fprintf('\n  <8> Neptune');
    fprintf('\n  <9> Pluto');

    if (i == 1)

        while(1)

            fprintf('\n\nplease select the departure planet\n');

            ip1 = input('? ');

            if (ip1 >= 1 && ip1 <= 9)

                break;

            end

        end

    end

    if (i == 2)

        while(1)

            fprintf('\n\nplease select the arrival planet\n');

            ip2 = input('? ');

            if (ip2 >= 1 && ip2 <= 9)

                break;

            end

        end

    end

end

% compute classical orbital elements of departure and arrival planets

[ri, vi] = planet(ip1, jdate1);

oevi = eci2orb1(mu, ri, vi);

[rf, vf] = planet(ip2, jdate2);

oevf = eci2orb1(mu, rf, vf);

% solve Lambert's problem using Gooding's method

sv1 = zeros(6, 1);

sv2 = zeros(6, 1);

sv1(1:3) = ri(1:3);

sv1(4:6) = vi(1:3);

sv2(1:3) = rf(1:3);

sv2(4:6) = vf(1:3);

[vito, vfto] = lambert_gooding(mu, sv1, sv2, tof, 0);

% initial orbital elements of the transfer orbit

oevtoi = eci2orb1(mu, ri, vito);

% final orbital elements of the transfer orbit

oevtof = eci2orb1(mu, rf, vfto);

% initial delta-v vector (kilometers/second)

dvi(1) = vito(1) - vi(1);
dvi(2) = vito(2) - vi(2);
dvi(3) = vito(3) - vi(3);

% final delta-v vector (kilometers/second)

dvf(1) = vf(1) - vfto(1);
dvf(2) = vf(2) - vfto(2);
dvf(3) = vf(3) - vfto(3);

% convert time and dates to strings

[cdstr1, utstr1] = jd2str(jdate1);

[cdstr2, utstr2] = jd2str(jdate2);

% display results

fprintf('\n\n         program lambert2\n');

fprintf('\n< interplanetary lambert problem >\n');

fprintf('\ndeparture planet             ');

disp(string(pname(ip1)));

fprintf('departure calendar date      ');

disp(cdstr1);

fprintf('departure universal time     ');

disp(utstr1);

fprintf('\ndeparture julian date        %12.6f', jdate1);

fprintf('\n\narrival planet               ');

disp(string(pname(ip2)));

fprintf('arrival calendar date        ');

disp(cdstr2);

fprintf('arrival universal time       ');

disp(utstr2);

fprintf('\narrival julian date          %12.6f', jdate2');

fprintf('\n\ntransfer time              %12.6f  days \n ', taud);

fprintf('\nheliocentric ecliptic orbital elements of the departure planet');
fprintf('\n--------------------------------------------------------------\n');

oeprint1(mu, oevi);

fprintf('\nheliocentric ecliptic orbital elements of the transfer orbit after the initial delta-v');
fprintf('\n--------------------------------------------------------------------------------------\n');

oeprint1(mu, oevtoi);

fprintf('\nheliocentric ecliptic orbital elements of the transfer orbit prior to the final delta-v');
fprintf('\n---------------------------------------------------------------------------------------\n');

oeprint1(mu, oevtof);

fprintf('\nheliocentric ecliptic orbital elements of the arrival planet');
fprintf('\n------------------------------------------------------------\n');

oeprint1(mu, oevf);

fprintf('\ninitial delta-v vector and magnitude\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dvi(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dvi(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dvi(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n', 1000.0 * norm(dvi));

fprintf('\nenergy                      %12.6f  kilometers^2/seconds^2\n', norm(dvi) * norm(dvi));

fprintf('\nfinal delta-v vector and magnitude\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dvf(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dvf(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dvf(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n', 1000.0 * norm(dvf));

fprintf('\nspecific orbital energy     %12.6f  kilometers2/seconds^2\n', norm(dvf) * norm(dvf));

%%%%%%%%%%%%%%%%%%%%%%%
% trajectory graphics %
%%%%%%%%%%%%%%%%%%%%%%%

while(1)

    fprintf('\nwould you like to display this trajectory (y = yes, n = no)\n');

    slct = input('? ', 's');

    if (slct == 'y' || slct == 'n')

        break;

    end

end

if (slct == 'y')

    while(1)

        fprintf('\nplease input the graphics step size (days)\n');

        deltat = input('? ');

        if (deltat > 0.0)

            break;

        end

    end

    % planet periods

    [r1, v1] = planet(ip1, jdate1);

    oev1 = eci2orb1(mu, r1, v1);

    period1 = 2.0 * pi * oev1(1) * sqrt(oev1(1) / mu) / 86400;

    [r2, v2] = planet(ip2, jdate1);

    oev2 = eci2orb1(mu, r2, v2);

    period2 = 2.0 * pi * oev2(1) * sqrt(oev2(1) / mu) / 86400;

    xve = oev1(1) / aunit;

    if (oev2(1) > oev1(1))

        xve = oev2(1) / aunit;

    end

    % determine number of data points to plot

    npts1 = fix(period1 / deltat);

    npts2 = fix(period2 / deltat);

    npts3 = fix((jdate2 - jdate1) / deltat);

    [rti, vti] = orb2eci(mu, oevtoi);

    % create departure planet orbit data points

    x1 = zeros(npts1 + 2, 1);

    y1 = zeros(npts1 + 2, 1);

    for i = 0:1:npts1

        jdate = jdate1 + i * deltat;

        [r1, ~] = planet(ip1, jdate);

        x1(i + 1) = r1(1) / aunit;

        y1(i + 1) = r1(2) / aunit;

    end

    % compute last data point of departure planet

    [r1, v1] = planet(ip1, jdate1 + period1);

    x1(npts1 + 2) = r1(1) / aunit;

    y1(npts1 + 2) = r1(2) / aunit;

    % create arrival planet orbit data points

    x2 = zeros(npts2 + 2, 1);

    y2 = zeros(npts2 + 2, 1);

    for i = 0:1:npts2

        jdate = jdate1 + i * deltat;

        [r2, ~] = planet(ip2, jdate);

        x2(i + 1) = r2(1) / aunit;

        y2(i + 1) = r2(2) / aunit;

    end

    % compute last data point of arrival planet

    [r2, ~] = planet(ip2, jdate1 + period2);

    x2(npts2 + 2) = r2(1) / aunit;

    y2(npts2 + 2) = r2(2) / aunit;

    % create transfer orbit data points

    x3 = zeros(npts3 + 2, 1);

    y3 = zeros(npts3 + 2, 1);

    for i = 0:1:npts3

        tau = 86400.0 * i * deltat;

        [rft, ~] = twobody2(mu, tau, ri, vito');

        x3(i + 1) = rft(1) / aunit;

        y3(i + 1) = rft(2) / aunit;

    end

    % compute last data point

    tau = 86400.0 * (jdate2 - jdate1);

    [rft, vft] = twobody2(mu, tau, rti, vti);

    x3(npts3 + 2) = rft(1) / aunit;

    y3(npts3 + 2) = rft(2) / aunit;

    % plot planet orbits and transfer trajectory

    hold on;

    plot(x1, y1, '.b');

    plot(x1, y1, '-b');

    plot(x2, y2, '.g');

    plot(x2, y2, '-g');

    plot(x3, y3, '.r');

    plot(x3, y3, '-r');

    % plot sun

    plot(0, 0, 'hy', 'MarkerSize', 10, 'MarkerFaceColor', 'y');

    axis equal;

    % plot and label vernal equinox direction

    line ([0.05, xve], [0, 0], 'Color', 'black');

    text(1.05 * xve, 0, '\Upsilon');

    % label departure and arrival locations

    [r2, ~] = planet(ip2, jdate1);

    plot(r2(1) / aunit, r2(2) / aunit, '*r');

    text(r2(1) / aunit + 0.05, r2(2) / aunit + 0.05, 'D');

    [r2, v2] = planet(ip1, jdate2);

    plot(r2(1) / aunit, r2(2) / aunit, '*r');

    text(r2(1) / aunit + 0.05, r2(2) / aunit + 0.05, 'A');

    plot(x3(1), y3(1), '*r');

    text(x3(1) + 0.05, y3(1) + 0.05, 'D');

    plot(x3(npts3 + 2), y3(npts3 + 2), '*r');

    text(x3(npts3 + 2) + 0.05, y3(npts3 + 2) + 0.05, 'A');

    % label launch and arrival dates

    text(0.85 * xve, -xve + 0.8, 'Departure', 'FontSize', 8);

    text(0.875 * xve, -xve + 0.7, pname(ip1), 'FontSize', 8);

    text(0.875 * xve, -xve + 0.6, cdstr1, 'FontSize', 8);

    text(0.85 * xve, -xve + 0.4, 'Arrival', 'FontSize', 8);

    text(0.875 * xve, -xve + 0.3, pname(ip2), 'FontSize', 8);

    text(0.875 * xve, -xve + 0.2, cdstr2, 'FontSize', 8);

    % label plot and axes

    xlabel('X coordinate (AU)', 'FontSize', 12);

    ylabel('Y coordinate (AU)', 'FontSize', 12);

    title('Interplanetary Lambert Problem', 'FontSize', 16);

    zoom on;

    % create trajectory graphics disk files
    
    print('lambert2.eps', '-depsc');

    print('lambert2.jpg', '-djpeg');

    print('lambert2.tif', '-dtiff');

end


