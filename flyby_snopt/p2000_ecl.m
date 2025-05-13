function [r, v] = p2000_ecl(ntarg, jdate)

% heliocentric planet/asteroid/comet state vector

% mean ecliptic and equinox of j2000 coordinate system

% input

%  jdate = TDB julian day
%  ntarg = "target" body

% output

%  r = heliocentric position vector (kilometers)
%  v = heliocentric velocity vector (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global smu aunit boev

if (ntarg == 10)

    % -----------------------------
    % asteroid/comet celestial body
    % -----------------------------

    rper = boev(1);

    ecc = boev(2);

    xinc = boev(3);

    argper = boev(4);

    raan = boev(5);

    jdpp = boev(6);

    % semimajor axis (au)

    sma = rper / (1.0 - ecc);

    % time since perhelion passage (seconds)

    tspp = 86400.0 * (jdate - jdpp);

    % compute mean anomaly (radians)

    manom = sqrt(smu / abs(aunit * sma)^3) * tspp;

    % solve Kepler's equation for true anomaly (radians)

    [~, tanom] = kepler1 (manom, ecc);

    % load orbital elements array

    oev(1) = aunit * sma;
    oev(2) = ecc;
    oev(3) = xinc;
    oev(4) = argper;
    oev(5) = mod(raan, 2.0 * pi);
    oev(6) = tanom;

    % determine heliocentric ecliptic state vector
    % (kilometers and kilometers/second)

    [r, v] = orb2eci(smu, oev);

    return

end

% compute heliocentric ecliptic state vector of planet
% (kilometers and kilometers/second)

ncent = 11;

[r, v] = jpleph_mice(jdate, ntarg, ncent, 'ECLIPJ2000');




