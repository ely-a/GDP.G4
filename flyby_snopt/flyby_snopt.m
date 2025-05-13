% flyby_snopt.m           April 6, 2025

% patched-conic, single gravity assist interplanetary
% trajectory design, analysis and optimization

% JPL ephemeris, MICE routines and SNOPT

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global rtd emu smu xmu aunit ip1 ip2 ip3 jdtdb0 jdtdbi1

global eq2000 ec2000 pmu req rsoi flyby_alt fbalt_lwr fbalt_upr

global dv_departure dvm_departure dv_arrival dvm_arrival vinfm_in vinfm_out

global c3departure c3arrival jdtdb1 jdtdb2 jdtdb3

global otype rito1 vito1 rito2 vito2 dvh

global vinf_in vinf_out jdtdb_soi_in rp2sc vp2sc

vf4 = zeros(3, 1);

% J2000 ecliptic-to-equatorial transformation matrix

eq2000(1, 1) =  1.000000000000000;
eq2000(1, 2) =  4.403600000000000e-007;
eq2000(1, 3) = -1.909190000000000e-007;

eq2000(2, 1) = -4.799660000000000e-007;
eq2000(2, 2) =  0.917482137087000;
eq2000(2, 3) = -0.397776982902000;

eq2000(3, 1) = 0.000000000000000;
eq2000(3, 2) = 0.397776982902000;
eq2000(3, 3) = 0.917482137087000;

% create J2000 equatorial-to-ecliptic transformation matrix

ec2000 = eq2000';

% Astronomical Unit (kilometers)

aunit = 149597870.691;

% radians to degrees conversion factor

rtd = 180.0 / pi;

% define "reference" tdb julian date (January 1, 2000)

jdtdb0 = julian(1, 1, 2000);

% furnish jpl ephemeris

cspice_furnsh('de430.bsp');

% define planet name vector

pdata = ['Mercury       '; 'Venus         '; 'Earth         '; ...
         'Mars          '; 'Jupiter       '; 'Saturn        '; ...
         'Uranus        '; 'Neptune       '; 'Pluto         ';];

pname = cellstr(pdata);

% planet gravitational constants (kilometer^3/second^2)

pmu(1) = 22032.08;
pmu(2) = 324858.599;
pmu(3) = 398600.436;
pmu(4) = 42828.314;
pmu(5) = 126712767.863;
pmu(6) = 37940626.063;
pmu(7) = 5794549.007;
pmu(8) = 6836534.064;
pmu(9) = 981.601;

% sun gravitational constant (kilometer^3/second^2)

xmu = 0.295912208285591149e-03;

smu = xmu * aunit^3 / 86400.0^2;

% earth gravitational constant (kilometer^3/second^2)

xmu = 0.899701152970881167e-09;

emu = xmu * aunit^3 / 86400.0^2;

% planet equatorial radius (kilometers)

req(1) = 2439.7;
req(2) = 6051.9;
req(3) = 6378.14;
req(4) = 3397.0;
req(5) = 71492.0;
req(6) = 60268.0;
req(7) = 25559.0;
req(8) = 24764.0;
req(9) = 1151.0;

% sphere-of-influence for each planet (kilometers)

rsoi(1) = 112409.906759936;
rsoi(2) = 616277.129297850;
rsoi(3) = 924647.107795632;
rsoi(4) = 577231.618115568;
rsoi(5) = 48204698.6886979;
rsoi(6) = 54553723.6086354;
rsoi(7) = 51771106.3741412;
rsoi(8) = 86696170.2674129;
rsoi(9) = 15110628.1503479;

% begin simulation

clc; home;

% read simulation definition data file

[filename, pathname] = uigetfile('*.in', 'Please select the input file to read');

[fid, jdtdbi1, ddays1, jdtdbi2, ddays2, jdtdbi3, ddays3, ip1, ip2, ip3, ...
    otype] = flyby_readdata(filename);

% ----------------------------------------------
% solve the patched-conic gravity assist problem
% ----------------------------------------------

% number of control variables

ncv = 3;

xg = zeros(ncv, 1);

% number of mission constraints

nmc = 2;

% control variables initial guesses
% (departure, flyby, and arrival tdb julian dates)

xg(1) = jdtdbi1 - jdtdb0;

xg(2) = jdtdbi2 - jdtdb0;

xg(3) = jdtdbi3 - jdtdb0;

% bounds on control variables

xlwr = zeros(ncv, 1);
xupr = zeros(ncv, 1);

xlwr(1) = xg(1) - ddays1;
xupr(1) = xg(1) + ddays1;

xlwr(2) = xg(2) - ddays2;
xupr(2) = xg(2) + ddays2;

xlwr(3) = xg(3) - ddays3;
xupr(3) = xg(3) + ddays3;

% bounds on objective function

flow = zeros(nmc + 1, 1);
fupp = zeros(nmc + 1, 1);

flow(1) = 0.0;
fupp(1) = +inf;

% "normalized" bounds on flyby altitude inequality constraint

flow(2) = fbalt_lwr / req(ip2);
fupp(2) = fbalt_upr / req(ip2);

% bounds on v-infinity matching (equality) constraint

flow(3) = 0.0;
fupp(3) = 0.0;

xmul = zeros(ncv, 1);

xstate = zeros(ncv, 1);

fmul = zeros(nmc + 1, 1);

fstate = zeros(nmc + 1, 1);

% use snopt to find optimal solution

[x, ~, ~, ~, ~] = snopt(xg, xlwr, xupr, xmul, xstate, ...
    flow, fupp, fmul, fstate, 'flyby_func');

% evaluate optimal solution

[~, g] = flyby_func(x);

% compute orientation of the departure hyperbola
% (Earth mean equator and equinox of J2000)

dveq1 = eq2000 * dv_departure;

dveqm1 = norm(dveq1);

decl1 = 90.0 - rtd * acos(dveq1(3) / dveqm1);

rasc1 = rtd * atan3(dveq1(2), dveq1(1));

% compute orientation of the arrival hyperbola
% (Earth mean equator and equinox of J2000)

dveq2 = eq2000 * dv_departure;

dveqm2 = norm(dveq2);

decl2 = 90.0 - rtd * acos(dveq2(3) / dveqm2);

rasc2 = rtd * atan3(dveq2(2), dveq2(1));

% compute tdb julian dates of optimal transfer

jdtdb1 =  jdtdb0 + x(1);

jdtdb2 =  jdtdb0 + x(2);

jdtdb3 =  jdtdb0 + x(3);

% convert solution julian dates to calendar dates and tdb times

[cdstr1, tdbstr1] = jd2str(jdtdb1);

[cdstr2, tdbstr2] = jd2str(jdtdb2);

[cdstr3, tdbstr3] = jd2str(jdtdb3);

% display results

fprintf('\nflyby_snopt - patched-conic gravity assist - SNOPT version');
fprintf('\n==========================================================\n');

fprintf('\ninput data file ==> ');

disp(filename);

if (otype == 1)
    
    fprintf('\nminimize departure delta-v\n');
    
end

if (otype == 2)
    
    fprintf('\nminimize arrival delta-v\n');
    
end

if (otype == 3)
    
    fprintf('\nminimize total delta-v\n');
    
end

fprintf('\ndeparture planet   ');

disp(string(pname(ip1)));

fprintf('flyby planet       ');

disp(string(pname(ip2)));

fprintf('arrival planet     ');

disp(string(pname(ip3)));

% departure conditions

fprintf('\nDEPARTURE CONDITIONS');
fprintf('\n====================\n');

display_epoch(jdtdb1, 1);

fprintf('\nheliocentric departure delta-v vector and magnitude');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

fprintf('\ndeparture delta-vx         %12.6f meters/second', 1000.0 * dv_departure(1));
fprintf('\ndeparture delta-vy         %12.6f meters/second', 1000.0 * dv_departure(2));
fprintf('\ndeparture delta-vz         %12.6f meters/second', 1000.0 * dv_departure(3));

fprintf('\n\ndeparture delta-v          %12.6f meters/second\n', 1000.0 * dvm_departure);

fprintf('\ndeparture hyperbola characteristics');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nasymptote right ascension  %12.6f degrees\n', rasc1);

fprintf('\nasymptote declination      %12.6f degrees\n', decl1);

fprintf('\ndeparture energy           %12.6f kilometer^2/second^2\n', c3departure);

fprintf('\nspacecraft post-impulse heliocentric orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev1 = eci2orb1(smu, rito1, vito1);

oeprint1(smu, oev1, 3);

svprint(rito1, vito1);

fprintf('departure planet heliocentric orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

[rp1, vp1] = p2000_ecl(ip1, jdtdb1);

oev1 = eci2orb1(smu, rp1, vp1);

oeprint1(smu, oev1, 3);

svprint(rp1, vp1);

% flyby conditions

fprintf('PATCHED-CONIC FLYBY CONDITIONS');
fprintf('\n==============================\n');

display_epoch(jdtdb2, 1);

fprintf('\ndeparture-to-flyby time    %12.6f days\n', jdtdb2 - jdtdb1);

fprintf('\nv-infinity in              %12.6f meters/second', 1000.0 * vinfm_in);

fprintf('\nv-infinity out             %12.6f meters/second\n', 1000.0 * vinfm_out);

fprintf('\nflyby altitude             %12.6f kilometers\n', flyby_alt);

% turn angles

tmp = req(ip2) * vinfm_in * vinfm_in / pmu(ip2);

tangle1 = 2.0 * asin(1.0 / (1.0 + tmp));

fprintf('\nmaximum turn angle         %12.6f degrees', tangle1 * rtd);

tmp = (flyby_alt + req(ip2)) * vinfm_in * vinfm_in / pmu(ip2);

tangle2 = 2.0 * asin(1.0 / (1.0 + tmp));

fprintf('\nactual turn angle          %12.6f degrees\n', tangle2 * rtd);

fprintf('\nheliocentric delta-v       %12.6f meters/second\n', 1000.0 * dvh);

% heliocentric deltavs

dvh_max = sqrt(pmu(ip2) / req(ip2));

fprintf('\nmax heliocentric delta-v   %12.6f meters/second\n', 1000.0 * dvh_max);

% planet-centered orbital elements at periapsis

flyby_rp = req(ip2) + flyby_alt;

oev_periapsis = flyby_hyper(pmu(ip2), vinf_in, vinf_out, flyby_rp);

[rp2sc_periapsis, vp2sc_periapsis] = orb2eci(pmu(ip2), oev_periapsis);

fprintf('\nspacecraft planet-centered orbital elements and state vector at periapsis');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oeprint1(pmu(ip2), oev_periapsis, 1);

svprint(rp2sc_periapsis, vp2sc_periapsis);

[bplane, ~, ~, ~] = rv2bplane(pmu(ip2), rp2sc_periapsis, vp2sc_periapsis);

fprintf('spacecraft b-plane coordinates at periapsis');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

bpprint(rp2sc_periapsis, vp2sc_periapsis, bplane);

cdecl_asy = cos(bplane(2));
crasc_asy = cos(bplane(3));

srasc_asy = sin(bplane(3));
sdecl_asy = sin(bplane(2));

fprintf('\nspacecraft heliocentric orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

tau = 86400.0 * (jdtdb2 - jdtdb1);

[rper, vper] = twobody2(smu, tau, rito1, vito1);

oev = eci2orb1(smu, rper, vper);

oeprint1(smu, oev, 3);

svprint(rper, vper);

fprintf('flyby planet heliocentric orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

[rp, vp] = p2000_ecl(ip2, jdtdb2);

oev = eci2orb1(smu, rp, vp);

oeprint1(smu, oev, 3);

svprint(rp, vp);

% arrival conditions

fprintf('ARRIVAL CONDITIONS');
fprintf('\n==================\n');

display_epoch(jdtdb3, 1);

fprintf('\nflyby-to-arrival time      %12.6f days\n', jdtdb3 - jdtdb2);

fprintf('\nheliocentric arrival delta-v vector and magnitude');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

fprintf('\narrival delta-vx            %12.6f meters/second', 1000.0 * dv_arrival(1));
fprintf('\narrival delta-vy            %12.6f meters/second', 1000.0 * dv_arrival(2));
fprintf('\narrival delta-vz            %12.6f meters/second', 1000.0 * dv_arrival(3));

fprintf('\n\narrival delta-v             %12.6f meters/second\n', 1000.0 * dvm_arrival);

fprintf('\narrival hyperbola characteristics');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nasymptote right ascension  %12.6f degrees\n', rasc2);

fprintf('\nasymptote declination      %12.6f degrees\n', decl2);

fprintf('\narrival energy             %12.6f kilometer^2/second^2\n', c3arrival);

% two-body propagation of the second leg of the transfer orbit

[rf3, vf3] = twobody2(smu, 86400 * (jdtdb3 - jdtdb2), rito2, vito2);

fprintf('\nspacecraft pre-impulse heliocentric orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev3 = eci2orb1(smu, rf3, vf3);

oeprint1(smu, oev3, 3);

svprint(rf3, vf3);

fprintf('spacecraft post-impulse heliocentric orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

rf4 = rf3;

for i = 1:1:3

    vf4(i) = vf3(i) + dv_arrival(i);

end

oev4 = eci2orb1(smu, rf4, vf4);

oeprint1(smu, oev4, 3);

svprint(rf4, vf4);

fprintf('arrival planet heliocentric orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

[rp, vp] = p2000_ecl(ip3, jdtdb3);

oev = eci2orb1(smu, rp, vp);

oeprint1(smu, oev, 3);

svprint(rp, vp);

fprintf('MISSION SUMMARY');
fprintf('\n===============\n');

fprintf('\ntotal delta-v              %12.6f meters/second\n', ...
    1000.0 * (dvm_departure + dvm_arrival));

fprintf('\ntotal energy               %12.6f kilometer^2/second^2\n', ...
    dvm_departure^2 + dvm_arrival^2);

fprintf('\ntotal mission duration     %12.6f days\n\n', jdtdb3 - jdtdb1);

% ---------------------------------------------
% two-body integration of the trajectory from
% first impulse to SOI entrance at flyby planet
% ---------------------------------------------

% set up for ode45

options = odeset('RelTol', 1.0e-11, 'AbsTol', 1.0e-12, 'Events', @soi_entrance);

% define maximum search time (seconds)

tof = 86400.0 * (jdtdb2 - jdtdb1);

% --------------------------------------------------
% solve for flyby planet SOI entrance conditions
% two-body forward integration of the trajectory
% from first impulse to SOI entrance at flyby planet
% --------------------------------------------------

[~, ~, tevent, ~, ~] = ode45(@two_body_heqm, [0 tof], [rito1 vito1], options);

jdtdb_soi_in = jdtdb1 + tevent / 86400.0;

fprintf('PATCHED-CONIC SOI ENTRANCE CONDITIONS');
fprintf('\n=====================================\n');

display_epoch(jdtdb_soi_in, 1);

fprintf('\nspacecraft planet-centered orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev_soi = eci2orb1(pmu(ip2), rp2sc, vp2sc);

oeprint1(pmu(ip2), oev_soi, 1);

svprint(rp2sc, vp2sc);

[bplane, ~, ~, ~] = rv2bplane(pmu(ip2), rp2sc, vp2sc);

fprintf('spacecraft b-plane coordinates at sphere-of-influence');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

bpprint(rp2sc, vp2sc, bplane);

rp2sc_in = rp2sc;

vp2sc_in = vp2sc;

% --------------------------------------------------
% solve for flyby planet SOI exit conditions
% two-body "backwards" integration of the trajectory
% from second impulse to SOI exit at flyby planet
% --------------------------------------------------

% set up for ode45

options = odeset('RelTol', 1.0e-11, 'AbsTol', 1.0e-12, 'Events', @soi_exit);

% define maximum search time (seconds)

tof = 86400.0 * (jdtdb2 - jdtdb3);

[~, ~, tevent, ~, ~] = ode45(@two_body_heqm, [0 tof], [rf3 vf3], options);

jdtdb_soi_out = jdtdb3 + tevent / 86400.0;

fprintf('\nPATCHED-CONIC SOI EXIT CONDITIONS');
fprintf('\n=================================\n');

display_epoch(jdtdb_soi_out, 1);

fprintf('\nspacecraft planet-centered orbital elements and state vector');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev_soi = eci2orb1(pmu(ip2), rp2sc, vp2sc);

oeprint1(pmu(ip2), oev_soi, 1);

svprint(rp2sc, vp2sc);

% -------------------------------
% heliocentric and flyby graphics
% -------------------------------

slct = 'm';

while(slct ~= 'y' && slct ~= 'n')
    
    fprintf('\nwould you like to display trajectory graphics for this mission (y = yes, n = no)\n');
    
    slct = input('? ', 's');
    
end

if (slct == 'y')
    
    while(1)
        
        fprintf('\nplease input the heliocentric plot step size (days)\n');
        
        deltat = input('? ');
        if (deltat > 0.0)

            break;

        end
                
    end
    
    while(1)
        
        fprintf('\nplease input the flyby plot span (days)\n');
        
        dt_flyby = input('? ');

        if (dt_flyby > 0.0)

            break;

        end
                
    end
    
    % heliocentric planetary orbits and transfer trajectory graphics
    % (mean ecliptic and equinox of j2000)
        
    % departure planet semimajor axis and period
    
    [r1, v1] = p2000_ecl(ip1, jdtdb1);
    
    oev1 = eci2orb1(smu, r1, v1);
    
    period1 = 2.0 * pi * oev1(1) * sqrt(oev1(1) / smu) / 86400.0;
    
    % flyby planet semimajor axis and period
    
    [r2, v2] = p2000_ecl(ip2, jdtdb1);
    
    oev2 = eci2orb1(smu, r2, v2);
    
    period2 = 2.0 * pi * oev2(1) * sqrt(oev2(1) / smu) / 86400.0;
    
    % arrival planet semimajor axis and period
    
    [r3, v3] = p2000_ecl(ip3, jdtdb1);
    
    oev3 = eci2orb1(smu, r3, v3);
    
    period3 = 2.0 * pi * oev3(1) * sqrt(oev3(1) / smu) / 86400.0;
    
    % number of "planet" points to plot
    
    npts1 = fix(period1 / deltat);
    
    npts2 = fix(period2 / deltat);
    
    npts3 = fix(period3 / deltat);
    
    % departure planet orbit
    
    x1 = zeros(1, npts1 + 1);
    
    y1 = zeros(1, npts1 + 1);
    
    for i = 0:1:npts1
        
        jdtdb = jdtdb1 + i * deltat;
        
        [r1, ~] = p2000_ecl(ip1, jdtdb);
        
        x1(i + 1) = r1(1) / aunit;
        
        y1(i + 1) = r1(2) / aunit;
        
    end
    
    % compute last data point
    
    [r1, v1] = p2000_ecl(ip1, jdtdb1 + period1);
    
    x1(npts1 + 1) = r1(1) / aunit;
    
    y1(npts1 + 1) = r1(2) / aunit;
    
    % flyby planet heliocentric orbit
    
    x2 = zeros(1, npts2 + 1);
    
    y2 = zeros(1, npts2 + 1);
    
    for i = 0:1:npts2
        
        jdtdb = jdtdb1 + i * deltat;
        
        [r2, ~] = p2000_ecl(ip2, jdtdb);
        
        x2(i + 1) = r2(1) / aunit;
        
        y2(i + 1) = r2(2) / aunit;
        
    end
    
    % compute last data point
    
    [r2, v2] = p2000_ecl(ip2, jdtdb1 + period2);
    
    x2(npts2 + 1) = r2(1) / aunit;
    
    y2(npts2 + 1) = r2(2) / aunit;
    
    % arrival planet heliocentric orbit
    
    x3 = zeros(1, npts3 + 1);
    
    y3 = zeros(1, npts3 + 1);
    
    for i = 0:1:npts3
        
        jdtdb = jdtdb1 + i * deltat;
        
        [r3, ~] = p2000_ecl(ip3, jdtdb);
        
        x3(i + 1) = r3(1) / aunit;
        
        y3(i + 1) = r3(2) / aunit;
        
    end
    
    % compute last data point
    
    [r3, v3] = p2000_ecl(ip3, jdtdb1 + period3);
    
    x3(npts3 + 1) = r3(1) / aunit;
    
    y3(npts3 + 1) = r3(2) / aunit;
    
    % first leg of transfer orbit
    
    npts4 = fix((jdtdb2 - jdtdb1) / deltat);
    
    x4 = zeros(1, npts4 + 1);
    
    y4 = zeros(1, npts4 + 1);
    
    for i = 0:1:npts4
        
        tau = 86400.0 * i * deltat;
        
        [rft1, ~] = twobody2(smu, tau, rito1, vito1);
        
        x4(i + 1) = rft1(1) / aunit;
        
        y4(i + 1) = rft1(2) / aunit;
        
    end
    
    % compute last data point of first leg
    
    tau = 86400.0 * (jdtdb2 - jdtdb1);
    
    [rft1, vft1] = twobody2(smu, tau, rito1, vito1);
    
    x4(npts4 + 1) = rft1(1) / aunit;
    
    y4(npts4 + 1) = rft1(2) / aunit;
    
    % second leg of heliocentric transfer orbit
    
    npts5 = fix((jdtdb3 - jdtdb2) / deltat);
    
    x5 = zeros(1, npts5 + 1);
    
    y5 = zeros(1, npts5 + 1);
    
    for i = 0:1:npts5
        
        tau = 86400.0 * i * deltat;
        
        [rft2, vft2] = twobody2(smu, tau, rito2, vito2);
        
        x5(i + 1) = rft2(1) / aunit;
        
        y5(i + 1) = rft2(2) / aunit;
        
    end
    
    % compute last data point of second leg
    
    tau = 86400.0 * (jdtdb3 - jdtdb2);
    
    [rfto2, vfto2] = twobody2(smu, tau, rito2, vito2);
    
    x5(npts5 + 1) = rfto2(1) / aunit;
    
    y5(npts5 + 1) = rfto2(2) / aunit;
    
    % determine vernal equinox "size"
    
    xve = oev1(1) / aunit;
    
    if (oev2(1) > oev1(1))
        
        xve = oev2(1) / aunit;
        
    end
    
    if (oev3(1) > oev2(1))
        
        xve = oev3(1) / aunit;
        
    end
    
    % plot heliocentric planet orbits and patched-conic transfer trajectory
    
    figure (1);
     
    hold on;
    
    plot(x1, y1, '.b');
    
    plot(x1, y1, '-b', 'LineWidth', 1.5);
    
    plot(x2, y2, '.g');
    
    plot(x2, y2, '-g', 'LineWidth', 1.5);
    
    plot(x3, y3, '.r');
    
    plot(x3, y3, '-r', 'LineWidth', 1.5);
    
    plot(x4, y4, '.k');
    
    plot(x4, y4, '-k', 'LineWidth', 1.5);
    
    plot(x5, y5, '.k');
    
    plot(x5, y5, '-k', 'LineWidth', 1.5);
    
    % plot sun

    plot(0, 0, 'hy', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    
    % plot and label vernal equinox direction
    
    line ([0.05, xve], [0, 0], 'Color', 'black');
    
    text(1.05 * xve, 0, '\Upsilon');
    
    % label departure, flyby and arrival locations
    
    plot(x4(1), y4(1), '*k');
    
    text(x4(1) + 0.06, y4(1) - 0.01, 'D');
    
    plot(x5(1), y5(1), '*k');
    
    text(x5(1) + 0.06, y5(1) + 0.01, 'F');
    
    plot(x5(npts5 + 1), y5(npts5 + 1), '*k');
    
    text(x5(npts5 + 1) + 0.06, y5(npts5 + 1) - 0.01, 'A');
    
    axis equal;
    
    % label plot and axes in astronomical units (AU)
    
    xlabel('X coordinate (AU)', 'FontSize', 14);
    
    ylabel('Y coordinate (AU)', 'FontSize', 14);
    
    title('Patched-conic Gravity Assist Trajectory', ...
        '(heliocentric planet and transfer trajectories)', 'FontSize', 16);
    
    zoom on;
    
    print('-dtiff', 'flyby_snopt1.tif');
    
    % ---------------------------------------------------------
    % create three-dimensional graphics of the flyby trajectory
    % ---------------------------------------------------------
    
    % unit pointing vector from the flyby planet to the sun
    
    [rp, vp] = p2000_ecl(ip2, jdtdb2);
    
    up2s(1) = -rp(1) / norm(rp);
    up2s(2) = -rp(2) / norm(rp);
    up2s(3) = -rp(3) / norm(rp);
    
    sun_axisx = [0, 5 * up2s(1)];
    sun_axisy = [0, 5 * up2s(2)];
    sun_axisz = [0, 5 * up2s(3)];
    
    % unit velocity vector of the flyby planet
    
    uv(1) = vp(1) / norm(vp);
    uv(2) = vp(2) / norm(vp);
    uv(3) = vp(3) / norm(vp);
    
    uv_axisx = [0, 5 * uv(1)];
    uv_axisy = [0, 5 * uv(2)];
    uv_axisz = [0, 5 * uv(3)];
    
    dtof = 43200.0 * dt_flyby;
    
    % backward two-body propagation
    
    [ri_ho, vi_ho] = twobody2(pmu(ip2), -dtof, rp2sc_periapsis, vp2sc_periapsis);
    
    deltat1 = 2.0 * dtof / 300;
    
    simtime1 = -deltat1;
    
    rp1_x = zeros(1, 301);
    rp1_y = zeros(1, 301);
    rp1_z = zeros(1, 301);
    
    for i = 1:1:301
        
        % forward two-body propagation
        
        simtime1 = simtime1 + deltat1;
        
        % create hyperbolic orbit "normalized" position vector
        
        [rf, vf] = twobody2(pmu(ip2), simtime1, ri_ho, vi_ho);
        
        rp1_x(i) = rf(1) / req(ip2);
        
        rp1_y(i) = rf(2) / req(ip2);
        
        rp1_z(i) = rf(3) / req(ip2);
        
    end
    
    % create axes vectors
    
    xaxisx = [1, 1.75];
    xaxisy = [0, 0];
    xaxisz = [0, 0];
    
    yaxisx = [0, 0];
    yaxisy = [1, 1.75];
    yaxisz = [0, 0];
    
    zaxisx = [0, 0];
    zaxisy = [0, 0];
    zaxisz = [1, 1.75];
    
    figure(2);
    
    hold on;
    
    grid on;
    
    % plot flyby planet
    
    [x, y, z] = sphere(24);
    
    h = surf(x, y, z);
    
    colormap([127/255 1 222/255]);
    
    set (h, 'edgecolor', [1 1 1]);
    
    % plot coordinate system axes
    
    plot3(xaxisx, xaxisy, xaxisz, '-r', 'LineWidth', 1.25);
    
    plot3(yaxisx, yaxisy, yaxisz, '-g', 'LineWidth', 1.25);
    
    plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1.25);
    
    % plot unit pointing vector to the sun
    
    plot3(sun_axisx, sun_axisy, sun_axisz, '-y', 'LineWidth', 2.0);
    
    % plot unit position vector of the flyby planet

    plot3(uv_axisx, uv_axisy, uv_axisz, '-k', 'LineWidth', 2.0);
    
    % plot planet-centered flyby hyperbolic orbit
    
    plot3(rp1_x, rp1_y, rp1_z, '-m', 'LineWidth', 1.25);
    
    % label incoming leg of flyby hyperbola
    
    plot3(rp1_x(1), rp1_y(1), rp1_z(1), '*m');
    
    % label periapsis of flyby hyperbola
    
    plot3(rp2sc_periapsis(1) / req(ip2), rp2sc_periapsis(2) / req(ip2), rp2sc_periapsis(3) / req(ip2), ...
        'om', 'MarkerSize', 4, 'MarkerFaceColor', 'm');
    
    % label plot in flyby planet radii (PR)
    
    xlabel('X coordinate (PR)', 'FontSize', 14);
    
    ylabel('Y coordinate (PR)', 'FontSize', 14);
    
    zlabel('Z coordinate (PR)', 'FontSize', 14);
    
    title('Patched-conic Gravity Assist Trajectory', '(planet-centered flyby trajectory)', 'FontSize', 16);
    
    axis equal;
    
    view(38, 30);
    
    rotate3d on;
    
    print('-dtiff', 'flyby_snopt2.tif');
    
end

slct = 'm';

while(slct ~= 'y' && slct ~= 'n')
    
    fprintf('\nwould you like to display flyby graphics for this mission (y = yes, n = no)\n');
    
    slct = input('? ', 's');
    
end

if (slct == 'y')

    while(1)
        
        fprintf('\nplease input the flyby plot step size (minutes)\n');
        
        deltat = input('? ');

        if (deltat > 0.0)

            break;

        end
                
    end

    % -----------------------------------------------------
    % display planet-centered orbital elements during flyby
    % -----------------------------------------------------

    % propagate backwards for 12 hours

    ti = 0.0;

    tof = -43200.0;

    npts = 1;

    % graphics step size (seconds)

    dtstep = 60.0 * deltat;

    neq = 6;

    tetol = 1.0e-12;

    % planet-centered state vector at periapsis of flyby planet
    % (kilometers and kilometers/second)

    yi(1:3) = rp2sc_periapsis;

    yi(4:6) = vp2sc_periapsis;

    tplot(npts) = 0.0;

    oev_flyby = eci2orb1(pmu(ip2), yi(1:3), yi(4:6));

    oev_sma(npts) = oev_flyby(1);

    oev_ecc(npts) = oev_flyby(2);

    oev_inc(npts) = rtd * oev_flyby(3);

    oev_argper(npts) = rtd * oev_flyby(4);

    oev_raan(npts) = rtd * oev_flyby(5);

    oev_tanom(npts) = rtd * oev_flyby(6);

    oev_speed(npts) = norm(yi(4:6));

    while(1)

        % step size guess (seconds)

        h = 30.0;

        % increment final time (seconds)

        tf = ti - dtstep;

        % check for last step size

        if (tf < tof)

            tf = tof;

        end 

        % integrate from ti to tf

        yfinal = rkf78('flyby_eqm', neq, ti, tf, h, tetol, yi);

        % compute heliocentric classical orbital elements

        oev_flyby = eci2orb1(pmu(ip2), yfinal(1:3), yfinal(4:6));

        npts = npts + 1;

        % tplot(npts) = tf / 3600.0 - 24.0 * (jdtdb2 - jdtdb_soi_in);

        tplot(npts) = tf / 3600.0;

        oev_sma(npts) = oev_flyby(1);

        oev_ecc(npts) = oev_flyby(2);

        oev_inc(npts) = rtd * oev_flyby(3);

        oev_argper(npts) = rtd * oev_flyby(4);

        oev_raan(npts) = rtd * oev_flyby(5);

        oev_tanom(npts) = rtd * oev_flyby(6);

        oev_speed(npts) = norm(yfinal(4:6));

        % reset initial time and state

        ti = tf;

        yi = yfinal;

        % check for end of current trajectory phase

        if (tf <= tof)

            break;

        end

    end

    % propagate forward for 12 hours

    ti = 0.0;

    tof = 43200.0;

    % planet-centered state vector at periapsis of flyby
    % (kilometers and kilometers/second)

    yi(1:3) = rp2sc_periapsis;

    yi(4:6) = vp2sc_periapsis;

    tplot(npts) = 0.0;

    oev_flyby = eci2orb1(pmu(ip2), yi(1:3), yi(4:6));

    oev_sma(npts) = oev_flyby(1);

    oev_ecc(npts) = oev_flyby(2);

    oev_inc(npts) = rtd * oev_flyby(3);

    oev_argper(npts) = rtd * oev_flyby(4);

    oev_raan(npts) = rtd * oev_flyby(5);

    oev_tanom(npts) = rtd * oev_flyby(6);

    oev_speed(npts) = norm(yi(4:6));

    while(1)

        % step size guess (seconds)

        h = 30.0;

        % increment final time (seconds)

        tf = ti + dtstep;

        % check for last step size

        if (tf > tof)

            tf = tof;

        end 

        % integrate from ti to tf

        yfinal = rkf78('flyby_eqm', neq, ti, tf, h, tetol, yi);

        % compute heliocentric classical orbital elements

        oev_flyby = eci2orb1(pmu(ip2), yfinal(1:3), yfinal(4:6));

        npts = npts + 1;

        tplot(npts) = tf / 3600.0;

        oev_sma(npts) = oev_flyby(1);

        oev_ecc(npts) = oev_flyby(2);

        oev_inc(npts) = rtd * oev_flyby(3);

        oev_argper(npts) = rtd * oev_flyby(4);

        oev_raan(npts) = rtd * oev_flyby(5);

        oev_tanom(npts) = rtd * oev_flyby(6);

        oev_speed(npts) = norm(yfinal(4:6));

        % reset initial time and state

        ti = tf;

        yi = yfinal;

        % check for end of current trajectory phase

        if (tf >= tof)

            break;

        end

    end

    %-------------------------------------------------
    % plot planet-centered semimajor axis during flyby
    %-------------------------------------------------

    figure(3);

    plot(tplot, oev_sma, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'b');

    xlabel('time since closest approach (hours)', 'FontSize', 14);

    ylabel('semimajor axis (kilometers)', 'FontSize', 14);

    title('Patched-conic Gravity-assist Trajectory', ...
        '(planet-centered semimajor axis during flyby)', 'FontSize', 16);

    grid on;

    print ('flyby_sma.tif', '-dtiff');

    %-----------------------------------------------
    % plot planet-centered eccentricity during flyby
    %-----------------------------------------------

    figure(4);

    plot(tplot, oev_ecc, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'b');

    xlabel('time since closest approach (hours)', 'FontSize', 14);

    ylabel('orbital eccentricity', 'FontSize', 14);

    title('Patched-conic Gravity-assist Trajectory', ...
        '(planet-centered orbital eccentricity during flyby)', 'FontSize', 16);

    grid on;

    print ('flyby_ecc.tif', '-dtiff');

    %----------------------------------------------
    % plot planet-centered inclination during flyby
    %----------------------------------------------

    figure(5);

    plot(tplot, oev_inc, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'b');

    xlabel('time since closest approach (hours)', 'FontSize', 14);

    ylabel('orbital inclination (degrees)', 'FontSize', 14);

    title('Patched-conic Gravity-assist Trajectory', ...
        '(planet-centered orbital inclination during flyby)', 'FontSize', 16);

    grid on;

    print ('flyby_inc.tif', '-dtiff');

    %------------------------------------------------------
    % plot heliocentric argument of perihelion during flyby
    %------------------------------------------------------

    figure(6);

    plot(tplot, oev_argper, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'b');

    xlabel('time since closest approach (hours)', 'FontSize', 14);

    ylabel('argument of periapsis (degrees)', 'FontSize', 14);

    title('Patched-conic Gravity-assist Trajectory', ...
        '(planet-centered argument of periapsis during flyby)', 'FontSize', 16);

    grid on;

    print ('flyby_argper.tif', '-dtiff');

    %------------------------------------------------------------------------
    % plot planet-centered right ascension of the ascending node during flyby
    %------------------------------------------------------------------------

    figure(7);

    plot(tplot, oev_raan, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'b');

    xlabel('time since closest approach (hours)', 'FontSize', 14);

    ylabel('RAAN (degrees)', 'FontSize', 14);

    title('Patched-conic Gravity-assist Trajectory', ...
        '(planet-centered RAAN during flyby)', 'FontSize', 16);

    grid on;

    print ('flyby_raan.tif', '-dtiff');

    %--------------------------------------------
    % plot planet-centered true anomaly during flyby
    %--------------------------------------------

    figure(8);

    plot(tplot, oev_tanom, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'b');

    xlabel('time since closest approach (hours)', 'FontSize', 14);

    ylabel('true anomaly (degrees)', 'FontSize', 14);

    title('Patched-conic Gravity-assist Trajectory', ...
        '(planet-centered true anomaly during flyby)', 'FontSize', 16);

    grid on;

    print ('flyby_tanom.tif', '-dtiff');

    %----------------------------------------
    % plot planet-centered speed during flyby
    %----------------------------------------

    figure(9);

    plot(tplot, oev_speed, 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'b');

    xlabel('time since closest approach (hours)', 'FontSize', 14);

    ylabel('speed (kilometers/second)', 'FontSize', 14);

    title('Patched-conic Gravity-assist Trajectory', ...
        '(planet-centered speed during flyby)', 'FontSize', 16);

    grid on;

    print ('flyby_speed.tif', '-dtiff');

end

fprintf('\n\n');

% CSPICE_KCLEAR clears the KEEPER system: unload all kernels, clears
% the kernel pool, and re-initialize the system.

cspice_kclear;

% close all non-cspice files

fclose('all');
