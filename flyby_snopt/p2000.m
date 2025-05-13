function [r, v] = p2000(ntarg, jdtdb)

% heliocentric celestial body state vector

% source ephemeris ==> jpl DE421

% input

%  jdtdb = tdb julian date

%  ntarg = celestial "target" body

%  NOTE: the numbering convention for 'ntarg' is:

%        1 = mercury           8 = neptune
%        2 = venus             9 = pluto
%        3 = earth            10 = moon
%        4 = mars             11 = asteroid/comet
%        5 = jupiter
%        6 = saturn
%        7 = uranus

% output

%  r = position vector (kilometers)

%  v = velocity vector (kilometers/second)

%      earth ecliptic j2000 coordinate system 

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global au smu eq2000 boev

if (ntarg == 10)
    
   % asteroid/comet
   
   rper = boev(1);

   ecc = boev(2);

   xinc = boev(3);

   argper = boev(4);

   raan = boev(5);

   jdpp = boev(6);

   % semimajor axis (au)

   sma = rper / (1.0 - ecc);

   % time since perhelion passage (seconds)

   tspp = 86400.0 * (jdtdb - jdpp);

   % compute mean anomaly (radians)

   manom = sqrt(smu / abs(au * sma)^3) * tspp;

   % solve Kepler's equation for true anomaly

   [~, tanom] = kepler1 (manom, ecc);

   % load orbital elements array

   oev(1) = au * sma;
   oev(2) = ecc;
   oev(3) = xinc;
   oev(4) = argper;
   oev(5) = mod(raan, 2.0 * pi);
   oev(6) = tanom;

   % determine heliocentric state vector

   [r, v] = orb2eci(smu, oev);

   return
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% celestial body is a planet or pluto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call jpl ephemeris NOTE: all jpl output vectors are
% referenced to the earth mean equator and equinox of j2000

ncent = 11;

result = jplephem (jdtdb, ntarg, ncent);

% load position and velocity vectors

rtmp = result(1: 3);
   
vtmp = result(4: 6);

% convert to Earth ecliptic and equinox j2000

r = eq2000' * rtmp;

v = eq2000' * vtmp;
   


